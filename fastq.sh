#!/usr/bin/env bash

#this is Jose's script
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err
#SBATCH -p normal,hns,dpetrov
#SBATCH -n 1
#SBATCH -t 2-00:00
#SBATCH --mem=4G
#SBATCH --requeue
#SBATCH --mail-user=alyulina@stanford.edu
#SBATCH --mail-type=ALL

#go to where the zipped reads are
cd /scratch/users/alyulina/novogene_01.09.2023_X202SC22123847-Z01-F001/01.RawData/fastq

#load modules
module load biology
module load fastqc

#unzip fastq files
for i in *.fq.gz; do gzip -cd "$i" > "${i%.*}"; done

#run fastqc
fastqc *.fq

#remove unzipped fastq files
rm *.fq

#move the output files from /scratch to /home
mv *fastqc* /home/users/alyulina/pdac/qc
