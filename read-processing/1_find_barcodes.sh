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
sample=$(head -$SLURM_ARRAY_TASK_ID samples.txt | tail -1)

# go to where the data is
cd "${data}"/"$sample"/

#!/bin/bash

# unzip all .fq.gz files
for i in *.fq.gz; do
    gzip -cd "$i" > "${i%.gz}"
done

# get unique sample prefixes
for file in *_L*_*.fq; do
    sample=$(echo "$file" | awk -F'_' '{OFS="_"; print $1, $2}' | sed 's/_$//')
    
    # merge reads from multiple lanes
    for readnum in 1 2; do
        cat ${sample}_*_L*_${readnum}.fq > "${sample}_R${readnum}.fq"
    done
done

# clean up: remove original unmerged .fq files
# rm *_L*_*.fq

# go back
cd "$dir"

srun --cpu_bind=verbose find_barcodes.py -s "$sample" -p "${data}"/

cp "${data}"/"$sample"/*_clID_rBC_extracted.txt ./out
cp "${data}"/"$sample"/*_failed_clIDs.txt ./out
cp "${data}"/"$sample"/*_find_barcodes_stats.txt ./out
