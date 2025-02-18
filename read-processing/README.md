## Extracting clonal and cell line barcodes

**Requirements:** `python 3.6.1` or above, `numpy 1.18.1` or above, `re 2.2.1` or above, `bartender` (Zhao *et al.*, 2017).

#### 1. Finding and extracting the barcode regions from raw reads  
The code in `find_barcodes.py` finds and extracts the barcode region from raw reads, if the following conditions are satisfied:
1. The forward and reverse reads match the regular expressions `TT([ATGC]{4})' + '([ATGC]{3}G[ATGC]{3}G[ATGC]{3}G[ATGC]{3}G)(GAAAC){e<2}` and `(GTTTC){e<2}(C[ATGC]{3}C[ATGC]{3}C[ATGC]{3}C[ATGC]{3})([ATGC]{4})AA`, respectively.
2. The barcode region was sequenced with the average error rate of at most 1 error per 1000 nucleotides.
3. The four-nucleotide cell-line identifier (clID) bolongs to those known *a priori*.  
4. There is a perfect match between the two reads.

Running this code requires a list of demultiplexed samples (see `samples.txt` for an example) + paired-end reads for each sample (see `raw-reads`). If running on a cluster with the `slurm` job manager, simply execute 
```
sbatch 1_find_barcodes.sh
```
The script will first unzip the `.fq.gz` files, merge reads that came from from different lanes, rename the two resulting files to `sample_R1.fq` and `sample_R2.fq`, and finally, run `find_barcodes.py`. This should produce the following output files for each sample:  
`./out/sample_find_barcodes_stats.txt` with stats on the number of reads that were filtered out;  
`./out/sample_failed_clIDs.txt` containing reads that had a mismatch in clID but satisfied all other requirements;  
`./out/sample_clIDs_rBC_extracted.txt`, which contains clIDs and clonal barcodes from all reads;  
`./out/sample_clID_rBC_extracted.txt`, a file with clonal barcodes for each clID. 
    
#### 2. Clustering and counting clonal barcodes within each cell line
The code in `cluster_barcodes.py` performs clusters extracted clonal barcodes within each cell line by running `bartender` with the following parameters: `-c 1 -s 1 -l 5 -z -1 -d 2`. Again, this can be done by executing
```
sbatch 2_cluster_barcodes.sh
```
This should produce the followint output: `sample_clID_rBC_cluster_counts.txt`, a file that contains how many times each barcode cluster was seen. You can find these counts in `../data/bc-counts` for all samples listed in `../data/metadata.csv`.
