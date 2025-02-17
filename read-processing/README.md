**Extracting barcodes and cell line IDs:**  
(change this to a readable readme)
To find and extract clonal barcodes and cell line IDs, execute `1_find_barcodes.sh`, which in turn will run `find_barcodes.py`.  

Running this code requires a list of demultiplexed samples (`samples.txt`), a folder for slurm output called `slurm`, and a folder for output called `out`. The `1_find_barcodes.sh` will first unzip the `.fq.gz` and then rename them to follow names from the list of samples and an `_R1` or `_R2` suffix – make sure that these steps work correctly! The renaming step only works for certain kinds of names so check the slurm output files to make sure that that step worked. Also, make sure to change the paths to data accordingly. If running for the first time, install regex first w/  
```
pip3 install --user regex
```

Doing the above should produce four kinds of output files:  
    (i) `<sample>_find_barcodes_stats.txt` w/ stats on the number of reads that were processed and filtered out;  
    (ii) `<sample>_failed_clIDs.txt` w/ output for reads that passed qc but had a mismatch in clID;  
    (iii) `<sample>_clID_bc_extracted.txt` –- a file with clIDs and barcodes from all reads (`.fq entry number \t clID \t barcode \n`);  
    (iv) `<sample>_clID_bc_extracted.txt` -- same as above but split into separate files for each clID.  
    
The first three files should be copied from `scratch` to `home` by the `1_find_barcodes.sh` script but this does not always work so make sure to double-check. All of the output is in the data folder.  

*Only pairs reads if both the barcode and the clID in both reads match.*
  
**Clustering and counting barcodes:**  
(make a readme)
