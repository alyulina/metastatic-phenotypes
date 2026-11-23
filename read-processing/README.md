## Extracting clonal and cell line barcodes

Processing has three steps: extracting barcodes from raw reads, clustering them within each cell line, and removing contaminating barcodes. Example inputs and outputs are provided here (`samples.txt`, `raw-reads`, `out`); the full raw reads are too large for this repository and can be accessed as described in the manuscript. The outputs of step 2 (`../data/bc-counts-clustered`) and step 3 (`../data/bc-counts-cleaned`, the removed barcodes, and the per-cell-line tables) are in `../data` for all samples in `../data/metadata.csv`.

**Requirements:** `python 3.6.1` or above (tested with `python 3.12.16`), `numpy 1.18.1` or above (tested with `numpy 2.5.5`), `re 2.2.1` or above, `pandas 3.1.1` or above, `scipy 1.18.1` or above, `matplotlib 3.11.2` or above, `bartender` (Zhao *et al.*, 2017).

#### 1. Finding and extracting the barcode regions from raw reads  
The code in `find_barcodes.py` finds and extracts the barcode region from raw reads, if the following conditions are satisfied:
1. The forward and reverse reads match the regular expressions `TT([ATGC]{4})' + '([ATGC]{3}G[ATGC]{3}G[ATGC]{3}G[ATGC]{3}G)(GAAAC){e<2}` and `(GTTTC){e<2}(C[ATGC]{3}C[ATGC]{3}C[ATGC]{3}C[ATGC]{3})([ATGC]{4})AA`, respectively.
2. The mean Phred quality score over the barcode region is at least 30 in both the forward and the reverse read.
3. The four-nucleotide cell-line identifier (clID) belongs to those known *a priori*.  
4. There is a perfect match between the two reads.

**Note on clID mismatches.** Expected mismatches: three clIDs are deliberately not in the list. Two belong to spike-ins that are not used (only one spike-in, `GATC`, is kept), and the third belongs to a cell line that was accidentally given the same clID as one of those excluded spike-ins (mT4-2D), so its reads cannot be told apart and it was excluded as well. In some samples these reads make up a large share of the reads that satisfy all other conditions.

Running this code requires a list of demultiplexed samples (see `samples.txt` for an example) + paired-end reads for each sample (see `raw-reads` for an example). If running on a cluster with the `slurm` job manager, simply execute 
```
sbatch 1_find_barcodes.sh
```
The script will first unzip the `.fq.gz` files, merge reads that came from different lanes, rename the two resulting files to `sample_R1.fq` and `sample_R2.fq`, and finally, run `find_barcodes.py`. For each sample, it writes to `out` the extracted clIDs and clonal barcodes (one file with all reads, and one file per clID, used in step 2), the reads with an unknown clID, and the number of reads filtered out at each step. 
    
#### 2. Clustering and counting clonal barcodes within each cell line
Sequencing and PCR errors can introduce small changes in barcode sequences, causing reads originating from the same clone to appear as distinct barcodes. We therefore group highly similar barcode sequences that differ substantially in abundance, since a low-abundance barcode that closely resembles a much more abundant barcode is likely to represent a spurious clone.

The code in `cluster_barcodes.py` clusters extracted clonal barcodes within each cell line by running `bartender` with the following parameters: `-c 1 -s 1 -l 5 -z -1 -d 2`. Again, this can be done by executing
```
sbatch 2_cluster_barcodes.sh
```
This should produce the following output: `sample_clIDs_rBC_cluster_counts.txt`, a file that contains how many times each barcode cluster was seen. You can find these counts in `../data/bc-counts-clustered` for all samples listed in `../data/metadata.csv`.

#### 3. Removing contaminating barcodes
Clustered counts still include barcodes that are not tumors of the sampled mouse: material from other samples, from the spike-in stock, or added during library preparation. `remove_contaminating_barcodes.ipynb` finds these sample–barcode pairs and removes them. It uses the following rules:

1. **Too large or too widespread at day two:** at day two, a barcode can have at most (cells injected) × (its frequency in the pre-injection pool) × 4 cells, i.e. two doublings. A barcode is flagged if its observed size is significantly larger than this maximum (one-sided Poisson test on the number of DNA molecules amplified, estimated from the spike-in; $p < 10^{-4}$).

   A barcode can also be found in too many mice: it is expected in a sample with probability $1 - e^{-\lambda}$, where $\lambda$ = (molecules of its cell line in that sample) × (its share of that cell line in the pre-injection pool; a barcode absent from the pool counts as one molecule); one found in more mice than expected (Poisson-binomial test, $p < 10^{-6}$) is removed from these samples.
2. **Co-travelling barcodes:** a single barcode shared between mice is often expected by chance, but a set of barcodes that appears in several samples at the same relative frequencies points to cross-contamination. If sample $m$ leaked into sample $n$, then $f_{i,n} \geq \alpha_{n,m} f_{i,m}$ for every barcode $i$ of sample $m$, where $f_{i,n}$ is the frequency of barcode $i$ in sample $n$; equality holds unless $n$ also has its own tumor with barcode $i$. Every pair of libraries from different mice with ≥ 20 shared barcodes is tested (each library on its own, since contamination can enter one library and not another):
   - $\alpha_{n,m}$ is the median of $f_{i,n} / f_{i,m}$ over shared barcodes (recomputed without barcodes more than 10× above the first estimate, which are likely tumors of sample $n$ itself);
   - $m$ is a source for $n$ if most barcodes of sample $m$ are found in $n$ (≥ 80% of those expected at this $\alpha_{n,m}$), and few shared barcodes (< 10%) are significantly less frequent in $n$ than $\alpha_{n,m} f_{i,m}$ (binomial test, $p < 0.01$);
   - the source has to stand out: $n$ should contain more barcodes from $m$ than from other samples of the same experiment (≥ 3× the median; a sample with most barcodes of every other sample, e.g. a day-two lung, has them by chance, since many more cells were injected);
   - if $m$ itself received material from another sample $k$, the barcodes of $k$ are left out of $m$ for this test, so that material from $k$ in both samples does not make $m$ look like a source (it is removed through $k$);
   - barcodes shared with a source are removed from $n$ unless they are significantly more frequent than $\alpha_{n,m} f_{i,m}$ ($p < 10^{-4}$) or more than $10 \times \alpha_{n,m} f_{i,m}$, i.e. likely its own tumors;
   - **spike-in contamination:** barcodes present in ≥80% of samples with spike-in, ≥100× less frequent in the pre-injection pools, and correlated with spike-in reads across samples (log scale, r ≥ 0.8) are removed from every sample of that experiment, before the tests above.

The latter shows that the experiment 1 spike-in stock probably contained other cell lines. As a result, only 29,600 of the 50,000 cells added were `GATC`, so cell numbers are scaled accordingly.

The tests above also identify cell line 7160c2 as a likely contaminant in experiment 1 (added during library preparation), so it is removed from all _in vivo_ samples of that experiment.

Pairs in which the source is from experiment 2 and the recipient from experiment 1 are tested as a negative control (experiment 2 was processed later); none passes. The thresholds were chosen so that the results agree with pairwise frequency plots inspected by eye, and each was varied in a sensitivity analysis.

The spike-in's own barcodes are never removed; _in vitro_ and pre-injection samples are unchanged. The cleaned counts (`sample_clIDs_rBC_cluster_counts_cleaned.txt`) can be found in `../data/bc-counts-cleaned` for all samples listed in `../data/metadata.csv`; they are used for all figures. The flagged pairs, with the rule(s) that flagged each and, for co-travelling barcodes, the source sample(s), are in `../data/flagged_barcodes.tsv`.

The notebook also builds `../data/cell-line_read_counts.csv` and `../data/cell-line_cell_counts.csv` from the cleaned counts. Figures are saved in `figs`.
