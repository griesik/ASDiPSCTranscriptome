#!/bin/bash

cd ~/results_tophat2
for i in $(ls -d */); do
cd $i
samtools sort -n accepted_hits.bam accepted_hits.sorted &
cd ..
done
wait

### Estimate abundance (HTSeq) ###

#Iniciando a estimação da abundância dos reads (running HTSeq)..."
cd ~/results_tophat2
for l in $(ls -d */); do
cd $l
nome=$(echo "$l" | cut -d "_" -f 1)
samtools view accepted_hits.sorted.bam | htseq-count -q - ~/databases/gencode.v19/gencode.v19.annotation.gtf > ~/htseq/$nome\_htseq.gencode19.out &
cd ..
done
wait
