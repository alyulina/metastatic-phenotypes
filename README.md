**Not pairing paired-end reads:**  
reads will be paired only in both the barcode and the clID in both reads match; this is now done in `find_barcodes.py`.

**Extracting barcodes and cell line IDs:**  
run `1_find_barcodes.sh` to run `find_barcodes.py` to find and extract barcodes and cell line IDs;
will generate four sets of output files: (i) a file with clIDs and barcodes from all reads (`.fq entry number \t clID \t barcode \n`); 
(ii) a file per cell line with just the barcodes within that cell line, which will be consequently used by bartender (unless there were no barcodes found);
(iii) a file with inferred clIDs that did not match those in the list of known clIDs; and (iv) a file that has information about the number of reads that were filtered out for different reasons. Running this requires a list of demultiplexed samples (`samples.txt`). The `1_find_barcodes.sh` will first unzip the `.fq.gz` and then rename them to follow names from the list of samples and an `_R1` or `_R2` suffix – make sure that these steps work correctly! Check the slurm output files to make sure that there was enough memory for the code to run. If running for the first time, install regex first w/ 
`pip3 install --user regex`.

**Clustering and counting barcodes:**  
work in progress!
