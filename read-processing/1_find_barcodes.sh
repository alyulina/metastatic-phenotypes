#!/bin/bash

#SBATCH --job-name=find_barcodes
#SBATCH --output=./slurm/slurm_find_barcodes_%A-%a.out
#SBATCH --partition hns,dpetrov,normal
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
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

# go to where the data is
cd "${data}"/"$sample_id"/

# note that this only works for some names but not all
ls *_1.fq.gz | cut -d _ -f 1-3 | uniq \
    | while read id; do \
        cat $id*_*_1.fq.gz > ${id%_*}_R1.fq.gz; \
        cat $id*_*_2.fq.gz > ${id%_*}_R2.fq.gz; \
    done
    
# unzip all .fq.gz files
for i in *.fq.gz; do
    gzip -cd "$i" > "${i%.gz}"
done

# clean up: remove original unmerged .fq files
# rm *_L*_*.fq

# go back
cd "$dir"

srun --cpu_bind=verbose python3 find_barcodes.py -s "$sample_id" -p "${data}"/

# cp "${data}"/"$sample_id"/*_clID_bc_extracted.txt ./out
# cp "${data}"/"$sample_id"/*_failed_clIDs.txt ./out
# cp "${data}"/"$sample_id"/*_find_barcodes_stats.txt ./out
