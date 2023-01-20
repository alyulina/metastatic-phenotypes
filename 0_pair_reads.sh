#!/bin/bash

#SBATCH --job-name=pair_reads
#SBATCH --output=/home/users/alyulina/pdac/bc_counts/novogene_01.09.2023_X202SC22123847-Z01-F001/slurm_pair_reads_%A-%a.out
#SBATCH --partition hns,dpetrov,normal
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --array=1-44
#SBATCH --no-requeue
#SBATCH --mail-user=alyulina@stanford.edu
#SBATCH --mail-type=ALL

echo $SLURM_ARRAY_TASK_ID

# read sample id
sample = $(head -$SLURM_ARRAY_TASK_ID input_samples.txt | tail -1)

# unzip .fastq files
for i in *.fq.gz; do gzip -cd "$i" > "${i%.*}"; done

forward_reads = "/scratch/users/alyulina/novogene_01.09.2023_X202SC22123847-Z01-F001/01.RawData/$sample/*_1.fastq"
reverse_reads = "/scratch/users/alyulina/novogene_01.09.2023_X202SC22123847-Z01-F001/01.RawData/$sample/*_2.fastq"

merged_reads = "/scratch/users/alyulina/novogene_01.09.2023_X202SC22123847-Z01-F001/01.RawData/$sample/$sample_merged"

# using pear to pair reads w/ -j 16 threads (same as ntasks) and -b 30 phred scores of at least q30
/home/groups/dpetrov/SOFTWARE/pear_0.9.11/bin/pear -f forward_reads -r reverse_reads -o merged_reads -j 16 -b 30

# this should produce four output files; we will need $sample_merged.assembled.fastq
