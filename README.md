**Extracting barcodes and cell line IDs:**  
run `1_find_barcodes.sh` to run `find_barcodes.py` to find and extract barcodes and cell line IDs. *Only pairs reads if both the barcode and the clID in both reads match.*  

When running on the cluster, change the data directory in the `1_find_barcodes.sh` script. Running this requires a list of demultiplexed samples (`samples.txt`), a folder for slurm output called `slurm`, and a folder for output called `out`. The `1_find_barcodes.sh` will first unzip the `.fq.gz` and then rename them to follow names from the list of samples and an `_R1` or `_R2` suffix – make sure that these steps work correctly! The renaming step only works for certain kinds of names so check the slurm output files to make sure that that step worked. If running for the first time, install regex first w/ 
`pip3 install --user regex`.  

Running this should produce four kinds of output files:  
    (i) `sample_find_barcodes_stats.txt` w/ stats on the number of reads that were processed and filtered out;  
    (ii) `sample_failed_clIDs.txt` w/ output for reads that passed qc but had a mismatch in clID;  
    (iii) `sample_clID_bc_extracted.txt` –- a file with clIDs and barcodes from all reads (`.fq entry number \t clID \t barcode \n`);  
    (iv) `sample_clID_bc_extracted.txt` -- same as above but split into separate files for each clID.  
    
The first three files should be copied `scratch` to `home` by the `1_find_barcodes.sh` script but this does not always work so make sure to double-check. All of the output is in the data folder.  
  
  
**Clustering and counting barcodes:**  
work in progress!
