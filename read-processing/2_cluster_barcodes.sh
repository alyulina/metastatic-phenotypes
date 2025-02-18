#!/bin/bash

#SBATCH --job-name=cluster_barcodes
#SBATCH --output=./slurm/slurm_find_barcodes_%A-%a.out
#SBATCH --partition hns,dpetrov,normal
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --array=1-191
#SBATCH --no-requeue
#SBATCH --mail-user=alyulina@stanford.edu
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=1

module load python/3.6.1
module load py-numpy/1.18.1_py36

echo $SLURM_ARRAY_TASK_ID

dir=$PWD
echo $dir

data="/scratch/users/alyulina/raw-reads"
echo "${data}"

# read sample id
sample_id=$(head -$SLURM_ARRAY_TASK_ID samples.txt | tail -1)

srun --cpu_bind=verbose python3 cluster_barcodes.py -s "$sample_id" -p "${data}"/
