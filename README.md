**Pairing paired-end reads:**  
run `0_pair_reads.sh` to run pear (https://github.com/tseemann/PEAR); 
requires a list of demultiplexed samples (`samples.txt`).

**Extracting barcodes and cell line IDs:**  
run `1_find_barcodes.sh` to run `find_barcodes.py` to find and extract barcodes and cell line IDs;
will generate two sets of output files: one file with all the extracted information 
and one file per cell line with just the barcodes within that cell line, which will be consequently used by bartender (unless there were no barcodes found); requires a list of demultiplexed samples (`samples.txt`); if running for the first time, install regex first w/ 
`pip3 install --user regex`. 

