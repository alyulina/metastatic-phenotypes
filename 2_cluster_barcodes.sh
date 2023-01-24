#!/bin/bash

#SBATCH --job-name=find_barcodes
#SBATCH --output=/home/users/alyulina/pdac/bc_counts/novogene_01.09.2023_X202SC22123847-Z01-F001/slurm/slurm_cluster_barcodes_%A-%a.out
#SBATCH --partition hns,dpetrov,normal
#SBATCH --time=1-00:00:00
#SBATCH --mem=8G
#SBATCH --array=1-44
#SBATCH --no-requeue
#SBATCH --mail-user=alyulina@stanford.edu
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=1

module load python/3.6.1
module load py-numpy/1.18.1_py36

echo $SLURM_ARRAY_TASK_ID
srun --cpu_bind=verbose python3 cluster_barcodes.py -s $(head -$SLURM_ARRAY_TASK_ID samples.txt | tail -1)
