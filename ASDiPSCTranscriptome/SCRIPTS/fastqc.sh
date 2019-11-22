#!/bin/bash

#### Qulity control (FastQC) ###

## access folder with raw data
cd ~/Fastqs/
## create a folder to save the results
mkdir ../fastqc/fastqc_results
for i in $(ls); do
if [[ $i = *.gz ]]
then
fastqc $i -o ~/fastqc/fastqc_results &
fi
done
wait

