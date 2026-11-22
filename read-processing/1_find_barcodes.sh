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

# merge lanes: works for any file naming, as long as mates end in _1.fq.gz / _2.fq.gz
shopt -s nullglob
r1_files=( *_1.fq.gz )
if [ ${#r1_files[@]} -eq 0 ]; then
    echo "No *_1.fq.gz files in $PWD" >&2; exit 1
fi
r2_files=()
for f in "${r1_files[@]}"; do
    r2="${f%_1.fq.gz}_2.fq.gz"
    [ -f "$r2" ] || { echo "Missing mate file: $r2" >&2; exit 1; }
    r2_files+=( "$r2" )
done
# decompress and concatenate lanes in matching order
gzip -cd "${r1_files[@]}" > "${sample_id}_R1.fq"
gzip -cd "${r2_files[@]}" > "${sample_id}_R2.fq"

# clean up: remove original unmerged .fq files
# rm *_L*_*.fq

# go back
cd "$dir"

srun --cpu_bind=verbose python3 find_barcodes.py -s "$sample_id" -p "${data}"/

# cp "${data}"/"$sample_id"/*_clID_bc_extracted.txt ./out
# cp "${data}"/"$sample_id"/*_failed_clIDs.txt ./out
# cp "${data}"/"$sample_id"/*_find_barcodes_stats.txt ./out
