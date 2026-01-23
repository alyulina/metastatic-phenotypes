#### This folder contains the following data:
- A table with **information for each sample in the study** (`metadata.csv`) with the following columns:  
    - `sample id`, a unique identifier given to each sample;  
   	- `experiment`, indicates whether this sample was from _in vitro_ passaging, _in vivo_ series 1, or _in vivo_ series 2 (the two _in vivo_ series differ primarily in how long cells were passaged before injected, beyond injection methods – see manuscript for details);
   	- `time point, d`, growth time in days (beginning at 0 days for pre-injection samples in both series);  
   	- `injection method`, can be `intravenous`, `intraperitoneal` (series 1), `intrasplenic` (series 2), `intrapancreatic` (series 2), and left blank (_in vitro_);  
   	- `tissue`, tissue collected (left blank for _in vitro_ samples);  
   	- `genotype`, mouse genotype (left blank for _in vitro_ samples);  
   	- `eartag (in vivo) / passage (in vitro)`, mouse eartag or passage replicate;
   	- `initial number of cells`, number of cells injected (_in vivo_) or transferred (_in vitro_); note that for day 0 samples, the initial cell number is _not_ the number of cells subsequently injected!
   	- `lung weight, g (in vivo) / cells (in vitro)`, measured lung weight in grams (_in vivo_) or estimated number of cell at the end of a passage (_in vitro_);  
   	- `spike-in added`, indicates whether a spike-in cell line was added during library preparation (see manuscript for details);
   	- `notes`, mostly used to indicate the replicate structure between samples (e.g. same library but different sequencing runs or same tissue but different libraries) and will include the colon-separated merged sample name that should be used for downstream processing.  

- **Clonal barcode counts** for all samples (`bc-counts/sample_clIDs_rBC_cluster_counts.txt`; see `metastatic-phenotypes/read-processing` for details on raw read processing).  

- A table containing **the number of reads associated with each cell line** in each sample, obtained by summing over clonal barcodes (`cell-line_read_counts.csv`).

- A table with **estimated cell counts for each cell line** (`cell-line_cell_counts.csv`). Barcode read counts were converted to cell numbers by normalizing to a barcoded spike-in cell line added at a known number during library preparation (50,000 cells per sample). Information from `metadata.csv` is additionally required to correctly merge reads for replicate samples. For details, see `utils.convert_barcode_reads_to_cell_counts()` in `../figures/utils.py`.


#### Raw data will be made available upon manuscript submission. Accession numbers and repository details will be provided in the manuscript.
