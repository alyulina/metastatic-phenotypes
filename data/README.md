#### This folder contains the following data:
- A table with information for each sample in the study (`metadata.csv`) with the following columns:  
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

- Clonal barcode counts for all samples (`bc-counts/sample_clIDs_rBC_cluster_counts.txt`, see `metastatic-phenotypes/read-processing` for details);  

- A table that contains the number of reads associated with each cell line in a sample (`cell-line_read_counts.csv`), as well as a table with estimated cell counts (`cell-line_cell_counts.csv`);  

- Histological (`histology_tumor_area.csv`) and lung weight (`tissue_weight.csv`) data for each cell line from single-cell-line transplant experiments (see Figure 1 in the manuscript for details).

#### Raw data will be made available upon manuscript submission. Accession numbers and repository details will be provided in the manuscript.
