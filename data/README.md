## This folder contains the data needed to generate Figures 2–4 in the main text
- **Sample metadata** (`metadata.csv`) contains information for each sample, including:  
    - `sample id`, a unique identifier;  
    - `experiment`, whether the sample is _in vitro_ passaging, _in vivo_ series 1, or _in vivo_ series 2 (the two series differ mainly in pre-injection culturing, see manuscript);  
   	- `time point, d`, growth time in days (pre-injection samples are day 0);  
   	- `injection method`, `tissue`, mouse `genotype`, `eartag (in vivo) / passage (in vitro)`;
   	- `initial number of cells`, number of cells injected (_in vivo_) or transferred (_in vitro_); note that for day 0 samples, this is _not_ the number of cells subsequently injected!
   	- `lung weight, g (in vivo) / cells (in vitro)`, measured lung weight in grams (_in vivo_) or estimated number of cells at the end of a passage (_in vitro_);  
   	- `spike-in added`, indicates whether a spike-in cell line was added during library preparation (see manuscript for details);
   	- `notes`, mostly used to indicate the replicate structure between samples (e.g. same library but different sequencing runs or same tissue but different libraries) and will include the colon-separated merged sample name that should be used for downstream processing.  

- **Clonal barcode counts** for all samples (`bc-counts/sample_clIDs_rBC_cluster_counts.txt`). All upstream data processing – from raw sequencing reads to clustered clonal barcode counts – is contained in `../read-processing`.  

These files are sufficient to generate most figures in the manuscript, except for i. two supplementary panels associated with Figure 2 that describe initial library cloning and transduction (these panels were generated separately by Saswati) and ii. two supplementary panels associated with Figure 3 that compare clustered and unclustered reads, which additionally require the raw data. All downstream analyses and figure generation are contained in the `../figures` folder. Commonly used downstream are  

- **The number of reads associated with each cell line** in each sample, obtained by summing over clonal barcodes (`cell-line_read_counts.csv`).

- **Estimated cell counts for each cell line** (`cell-line_cell_counts.csv`). Read counts were converted to cell numbers by normalizing to a barcoded spike-in cell line added at a known abundance during library preparation (50,000 cells per sample). Information from `metadata.csv` is additionally required to correctly merge reads for replicate samples. For details, see `utils.convert_barcode_reads_to_cell_counts()` in `../figures/utils.py`.

_Raw data will be made available upon manuscript submission. Accession numbers and repository details will be provided in the manuscript._
