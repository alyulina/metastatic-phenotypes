import os
import numpy as np
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-s",
                  "--sample",
                  type="string",
                  help="sample ID as in samples.txt",
                  dest="sample")

(options, args) = parser.parse_args()
sample = options.sample

# path to the folder w/ raw data for this sample
sample_path = '/scratch/users/alyulina/novogene_01.09.2023_X202SC22123847-Z01-F001/01.RawData/' + sample + '/'
bartender_path = '/home/groups/dpetrov/SOFTWARE/bartender-1.1-master/'

# cell line IDs, same as before, i.e.
# minumum Hamming distance between any two of them is 2
# excluded: all but one spike-ins & mt4-2D, which had the same clIDs
clIDs = ['GATC', 'AAGG', 'ACAC', 'ACCT', 'ACGA', 'ACTG', 'AGAG', 'AGCA', 'AGGT', 'AGTC', 'ATCG', 'ATGC', 'CAAC', 'CACT', 'CAGA', 'CATG', 'CCAA', 'CCTT', 'CGAT', 'CGTA', 'CTGT', 'CTTC', 'GAAG', 'GCAT', 'GCTA', 'GGAA', 'GGTT', 'GTAC', 'GTGA', 'GTTG', 'TCCA', 'TGAC', 'TTCC', 'TTGG']


clID_bc_count = [] # a list w/ clID_bc \t count \n

# want to have: a merged file for all cell lines for a sample w/ cell line ID, cluster center, and count for this cluster within each cell line

for i in clIDs:
    barcodes_path = sample_path + sample + '_' + i + '_rBC_extracted.txt'
    out_prefix = sample_path + sample + '_' + i + '_bartender'
    
    if not os.path.exists(barcodes_path):
        continue

    # else run bartender with the command line, no cutoff to not mess up with the line number -c 1, disable statistical test  -z -1, set up Hamming distance at 2 -d 2, seed length 8 -l, 8 threads
    os.system(f'{bartender_path}bartender_single_com -f {barcodes_path} -o {out_prefix} -c 1 -s 1 -l 5 -z -1 -d 2')

    cluster_id_cluster_seq_i = {} # a dictionary w/ cluster IDs as keys as cluster centers as values
    cluster_seq_count_i = {}  # a dictionary w/ cluster IDs as keys and counts as values

    # reading bartender output files
    for line in open(f'{out_prefix}_cluster.csv', 'r'):
        cluster_id, consensus_seq, score, t = line.strip().split(',')

        if cluster_id.isdigit():  # skipping the header line
            cluster_id_cluster_seq_i[cluster_id] = consensus_seq

    for line in open(f'{out_prefix}_barcode.csv', 'r'):
        bc, count, cluster_id = line.strip().split(',')

        if cluster_id.isdigit(): # skipping the header line
            consensus_seq = cluster_id_cluster_seq_i.get(cluster_id)
            if consensus_seq in cluster_seq_count_i:
                cluster_seq_count_i[consensus_seq] += int(count)
            else:
                cluster_seq_count_i[consensus_seq] = int(count)

    clID_bc_count.extend([i + '_' + x[0] + '\t' + str(x[1]) + '\n' for x in cluster_seq_count_i.items()])

with open(sample_path + sample + '_clIDs_rBC_cluster_counts.txt', 'w+') as o:
    o.writelines(clID_bc_count)
