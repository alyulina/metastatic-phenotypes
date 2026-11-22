## Extracting clonal and cell line barcodes

**Requirements:** `python 3.6.1` or above, `numpy 1.18.1` or above, `re 2.2.1` or above, `bartender` (Zhao *et al.*, 2017).

#### 1. Finding and extracting the barcode regions from raw reads  
The code in `find_barcodes.py` finds and extracts the barcode region from raw reads, if the following conditions are satisfied:
1. The forward and reverse reads match the regular expressions `TT([ATGC]{4})' + '([ATGC]{3}G[ATGC]{3}G[ATGC]{3}G[ATGC]{3}G)(GAAAC){e<2}` and `(GTTTC){e<2}(C[ATGC]{3}C[ATGC]{3}C[ATGC]{3}C[ATGC]{3})([ATGC]{4})AA`, respectively.
2. The mean Phred quality score over the barcode region is at least 30 in both the forward and the reverse read.
3. The four-nucleotide cell-line identifier (clID) bolongs to those known *a priori*.  
4. There is a perfect match between the two reads.

**Note on clID mismatches.** Expected mismatches: three clIDs are deliberately not in the list. Two belong to spike-ins that are not used (only one spike-in, `GATC`, is kept), and the third belongs to a cell line that was accidentally given the same clID as one of those excluded spike-ins (mT4-2D), so its reads cannot be told apart and it was excluded as well. In some samples these reads make up a large share of the reads that satisfy all other conditions.

Running this code requires a list of demultiplexed samples (see `samples.txt` for an example) + paired-end reads for each sample (see `raw-reads` for an example). If running on a cluster with the `slurm` job manager, simply execute 
```
sbatch 1_find_barcodes.sh
```
The script will first unzip the `.fq.gz` files, merge reads that came from from different lanes, rename the two resulting files to `sample_R1.fq` and `sample_R2.fq`, and finally, run `find_barcodes.py`. This should produce the following output files for each sample:  
`./out/sample_find_barcodes_stats.txt` with stats on the number of reads that were filtered out;  
`./out/sample_failed_clIDs.txt` containing reads that had a mismatch in clID but satisfied all other requirements;  
`./out/sample_clIDs_rBC_extracted.txt`, which contains clIDs and clonal barcodes from all reads;  
`./out/sample_clID_rBC_extracted.txt`, a file with clonal barcodes for each clID. 
    
#### 2. Clustering and counting clonal barcodes within each cell line
Sequencing and PCR errors can introduce small changes in barcode sequences, causing reads originating from the same clone to appear as distinct barcodes. We therefore group highly similar barcode sequences that differ substantially in abundance, since a low-abundance barcode that closely resembles a much more abundant barcode is likely to represent a spurious clone.

The code in `cluster_barcodes.py` performs clusters extracted clonal barcodes within each cell line by running `bartender` with the following parameters: `-c 1 -s 1 -l 5 -z -1 -d 2`. Again, this can be done by executing
```
sbatch 2_cluster_barcodes.sh
```
This should produce the following output: `sample_clIDs_rBC_cluster_counts.txt`, a file that contains how many times each barcode cluster was seen. You can find these counts in `../data/bc-counts-clustered` for all samples listed in `../data/metadata.csv`.

#### 3. Removing contaminating barcodes
Some barcodes in the _in vivo_ samples cannot come from tumors in that mouse. `remove_contaminating_barcodes.ipynb` flags these sample–barcode pairs and removes them from the clustered counts. It uses two rules:

1. **Too large too soon:** at day two, a barcode can have at most (cells injected) × (its frequency in the pre-injection pool) × 4 cells, i.e. two doublings. A barcode is flagged if its observed size is significantly larger than this maximum (one-sided Poisson test on the number of DNA molecules amplified, estimated from the spike-in; $p < 10^{-4}$).
2. **Co-travelling barcodes:** a single barcode shared between mice is often expected by chance, but a set of barcodes that appears in several samples at the same relative frequencies points to cross-contamination. If sample $m$ leaked into sample $n$, then $f_{i,n} \geq \alpha_{n,m} f_{i,m}$ for every barcode $i$ of sample $m$, where $f_{i,n}$ is the frequency of barcode $i$ in sample $n$; equality holds unless $n$ also has its own tumor with barcode $i$. Every pair of samples from different mice with ≥ 20 shared barcodes is tested:
   - $\alpha_{n,m}$ is the median of $f_{i,n} / f_{i,m}$ over shared barcodes (recomputed without barcodes more than 10× above the first estimate, which are likely tumors of sample $n$ itself);
   - $m$ is a source for $n$ if most barcodes of sample $m$ are found in $n$ (≥ 80% of those expected at this $\alpha_{n,m}$), and few shared barcodes (< 10%) are significantly less frequent in $n$ than $\alpha_{n,m} f_{i,m}$ (binomial test, $p < 0.01$);
   - barcodes shared with a source are removed from $n$ unless they are significantly more frequent than $\alpha_{n,m} f_{i,m}$ ($p < 10^{-4}$) or more than $10 \times \alpha_{n,m} f_{i,m}$, i.e. likely its own tumors;
   - **spike-in contamination:** barcodes present in ≥80% of samples with spike-in, ≥100× less frequent in the pre-injection pools, and correlated with spike-in reads across samples (log scale, r ≥ 0.8) are removed from every sample of that experiment, before the tests above.

Pairs in which the source is from experiment 2 and the recipient from experiment 1 are tested as a negative control (experiment 2 was processed later); none passes. The thresholds were chosen so that the results agree with x–y plots of sample pairs inspected by eye, and each was varied in a sensitivity analysis (see the notebook).

The spike-in's own barcodes are never removed; _in vitro_ and pre-injection samples are unchanged. You can find the cleaned counts, `sample_clIDs_rBC_cluster_counts_cleaned.txt`, in `../data/bc-counts-cleaned` for all samples listed in `../data/metadata.csv`; they are used for all figures. The flagged pairs, with the rule that flagged each, are in `../data/flagged_barcodes.tsv`.
