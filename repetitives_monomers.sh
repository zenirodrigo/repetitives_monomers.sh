#!/bin/bash

# Ativar autocompletar para nomes de arquivos
shopt -s progcomp

# Função para autocompletar
_autocomplete() {
    local cur=${COMP_WORDS[COMP_CWORD]}
    COMPREPLY=( $(compgen -f -- "$cur") )
}

# Atribuir a função de autocompletar ao comando read
complete -F _autocomplete read

# Pedir o nome do arquivo de biblioteca
read -e -p "Digite o nome do arquivo de biblioteca: " input_biblio

# Pedir a referência a ser usada
read -e -p "Digite a referência a ser usada: " referencia

# Pedir o nome de saída
read -e -p "Digite o nome do arquivo de saída: " saida

# Pedir o nome do arquivo de saída
read -e -p "Digite o nome do arquivo de saída, após retiradas as sequências superiores a 2 mil pares de base: " output_file
cobertura_minima=80

# Preparar biblioteca para blast
makeblastdb -in "$input_biblio" -dbtype nucl

# Executar o blastn
blastn -task blastn -outfmt 6 -db "$input_biblio" -query "$referencia" -out "$saida" -evalue 1e-50 -qcov_hsp_perc "$cobertura_minima"

# Extrair reads do output:
echo "Extraindo reads do output..."
cat "$saida" | awk '{print $2}' | sort -u > output_reads_extraidos

# Retirada dos reads:
echo "Retirando reads..."
seqtk subseq "$input_biblio" output_reads_extraidos > reads_extraidos.fasta

# Extrair coordenadas do output do BLAST:
echo "Extraindo coordenadas do output do BLAST..."
awk '{if ($2 > $3) print $1, $3, $2; else print $1, $2, $3}' "$saida" > coordenadas_temp.bed

# Ordenar coordenadas por cromossomo e posição inicial
sort -k1,1 -k2,2n coordenadas_temp.bed > coordenadas_sorted.bed

# Criar arquivo BED com coordenadas alinhadas de 2 em 2
echo "Criando arquivo BED..."
awk -v OFS='\t' '{
    cromossomo_atual = $1;
    posicao_inicial_atual = $2;

    if (NR % 2 == 0) {
        print cromossomo_anterior, posicao_inicial_anterior, posicao_inicial_atual;
    }

    cromossomo_anterior = cromossomo_atual;
    posicao_inicial_anterior = posicao_inicial_atual;
}' coordenadas_sorted.bed > coordenadas.bed

# Remover arquivos temporários
rm coordenadas_temp.bed coordenadas_sorted.bed

# Criar um arquivo FASTA final
output_final="reads_extraidos_final.fasta"
echo "Criando arquivo final..."
bedtools getfasta -fi "$input_biblio" -bed coordenadas.bed >> "$output_final"

# Adicionar cabeçalhos ao arquivo final
awk '{print ">"$2"\n"$0}' reads_extraidos.fasta >> "$output_final"

# Tamanho máximo desejado (em pares de bases)
tamanho_maximo=10000

# Inicializar variáveis
cabecalho=""
sequencia=""
contador=0

# Ler o arquivo de entrada linha por linha
while IFS= read -r linha; do
    if [ "${linha:0:1}" == ">" ]; then
        # Se a linha começar com ">", é um cabeçalho de sequência
        # Verificar se a sequência anterior é menor ou igual ao tamanho máximo
        if [ $contador -le $tamanho_maximo ] && [ -n "$cabecalho" ]; then
            echo -e "$cabecalho\n$sequencia" >> "$output_file"
        fi

        # Inicializar para a nova sequência
        cabecalho="$linha"
        sequencia=""
        contador=0
    else
        # Se não for um cabeçalho, é parte da sequência
        sequencia="$sequencia$linha"
        contador=$((contador + ${#linha}))
    fi
done < "$output_final"

# Verificar a última sequência no final do arquivo
if [ $contador -le $tamanho_maximo ] && [ -n "$cabecalho" ]; then
    echo -e "$cabecalho\n$sequencia" >> "$output_file"
fi

echo "Sequências com até 10000 pares de bases foram salvas em $output_file"
