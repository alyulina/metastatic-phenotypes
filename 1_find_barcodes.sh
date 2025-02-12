#!/bin/bash

#SBATCH --job-name=find_barcodes
#SBATCH --output=./slurm/slurm_find_barcodes_%A-%a.out
#SBATCH --partition hns,dpetrov,normal
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
#SBATCH --array=1-42
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

data="/scratch/users/alyulina/novogene_01.09.2023_X202SC22123847-Z01-F001/01.RawData"
echo "${data}"

# read sample id
sample=$(head -$SLURM_ARRAY_TASK_ID samples.txt | tail -1)

# go to where the data is
cd "${data}"/"$sample"/

# unzip .fq files
for i in *.fq.gz; do gzip -cd "$i" > "${i%.*}"; done

# rename unzipped .fq files
for i in *.fq; do cp $i $(echo $i | awk '{split($1,a,/_/); print a[1]"_"a[2]"_R"a[6]}'); done

# go back
cd "$dir"

srun --cpu_bind=verbose find_barcodes.py -s "$sample" -p "${data}"/

cp "${data}"/"$sample"/*_clID_bc_extracted.txt ./out
cp "${data}"/"$sample"/*_failed_clIDs.txt ./out
cp "${data}"/"$sample"/*_find_barcodes_stats.txt ./out
