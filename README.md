**Not pairing paired-end reads:**  
reads will be paired only in both the barcode and the clID in both reads match; this is now done in `find_barcodes.py`. When running on the cluster, change the data directory in the `1_find_barcodes.sh` script. Requires a list of samples called `samples.txt, a folder for slurm output called `slurm`, and a folder for output called `out`. Only works for certain kinds of names (for renaming in `1_find_barcodes.sh` to work).

Running this should produce four kinds of output files:  
    (i) `sample_find_barcodes_stats.txt` w/ stats on the number of reads that were processed and filtered out;  
    (ii) `sample_failed_clIDs.txt` w/ output for reads that passed qc but had a mismatch in clID;  
    (iii) `sample_clID_bc_extracted.txt` –- this is what you want -- it has clID and barcode extracted from reads w/ all clIDs in one file;  
    (iv) `sample_clID_bc_extracted.txt` -- same as above but split into separate files for each clID.  
    
The first three files should be copied `scratch` to `home` by the `1_find_barcodes.sh` script but this does not always work so make sure to double-check. All of the output is in the data folder.  
  

**Extracting barcodes and cell line IDs:**  
run `1_find_barcodes.sh` to run `find_barcodes.py` to find and extract barcodes and cell line IDs;
will generate four sets of output files: (i) a file with clIDs and barcodes from all reads (`.fq entry number \t clID \t barcode \n`); 
(ii) a file per cell line with just the barcodes within that cell line, which will be consequently used by bartender (unless there were no barcodes found);
(iii) a file with inferred clIDs that did not match those in the list of known clIDs; and (iv) a file that has information about the number of reads that were filtered out for different reasons. Running this requires a list of demultiplexed samples (`samples.txt`). The `1_find_barcodes.sh` will first unzip the `.fq.gz` and then rename them to follow names from the list of samples and an `_R1` or `_R2` suffix – make sure that these steps work correctly! Check the slurm output files to make sure that there was enough memory for the code to run. If running for the first time, install regex first w/ 
`pip3 install --user regex`.

**Clustering and counting barcodes:**  
work in progress!
