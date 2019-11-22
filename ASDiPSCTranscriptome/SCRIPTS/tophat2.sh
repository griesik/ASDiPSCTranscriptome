#!/bin/bash

### MAPEAMENTO (Tophat2, hg19) ###

mkdir ~/results_tophat2/
cd ~/Fastqs/
nome1="a"
nome_old="b"

for i in $(ls); do
nome1=$(echo "$i" | cut -d "." -f 1)

if [ "$nome1" == "$nome_old" ]
then
tophat2 -p 40 -o ~/results_tophat2/$nome1\_out --transcriptome-index=~/databases/gencode.v19/gencode.v19.annotation ~/Homo_sapiens_UCSC_hg19/UCSC/hg19/Sequence/Bowtie2Index/genome  $arquivo_old $i 2>> $nome1\_tophat2.out 
fi
nome_old=$nome1
arquivo_old=$i 

done








