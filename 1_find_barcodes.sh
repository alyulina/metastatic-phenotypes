#!/bin/bash

#SBATCH --job-name=find_barcodes
#SBATCH --output=/home/users/alyulina/pdac/bc_counts/novogene_04.26.2023_X202SC23040899-Z01-F001/slurm/slurm_find_barcodes_%A-%a.out
#SBATCH --partition hns,dpetrov,normal
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
#SBATCH --array=1-4
#SBATCH --no-requeue
#SBATCH --mail-user=alyulina@stanford.edu
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=1

module load python/3.6.1
module load py-numpy/1.18.1_py36

echo $SLURM_ARRAY_TASK_ID

# read sample id
sample=$(head -$SLURM_ARRAY_TASK_ID samples.txt | tail -1)

# go to where the data is
cd /scratch/users/alyulina/novogene_04.26.2023_X202SC23040899-Z01-F001/01.RawData/"$sample"/

# unzip .fq files
for i in *.fq.gz; do gzip -cd "$i" > "${i%.*}"; done

# rename unzipped .fq files
for i in *.fq; do mv $i $(echo $i | awk '{split($1,a,/_/); print a[1]"_"a[2]"_R"a[6]}'); done

srun --cpu_bind=verbose python3 /home/users/alyulina/pdac/find_barcodes.py -s $(head -$SLURM_ARRAY_TASK_ID samples.txt | tail -1) -p /scratch/users/alyulina/novogene_04.26.2023_X202SC23040899-Z01-F001/01.RawData/
