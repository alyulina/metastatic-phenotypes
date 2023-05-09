import regex as re
import numpy as np
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-s",
                  "--sample",
                  type="string",
                  help="sample ID as in samples.txt",
                  dest="sample")

parser.add_option("-p",
                  "--path",
                  type="string",
                  help="path to raw reads (without the last backslash)",
                  dest="path")

(options, args) = parser.parse_args()
sample = options.sample

# path to the folder w/ raw data for this sample
sample_path = path + '/' + sample + '/'

def avg_qscore(x):
    return np.mean([ord(i) - 33 for i in x])

# regular expression pattern for the barcode region: 
# regex_fw = re.compile('TT' + '(clID)' + '(barcode)' + '(GAAAC){e<2}'), 
# regex_rv is the reverse complement
# allowing at most one error in the subsequent region and no errors in the preceeding region + no Ns in the barcode or clID
regex_fw = re.compile('TT' + '([ATGC]{4})' + '([ATGC]{3}G[ATGC]{3}G[ATGC]{3}G[ATGC]{3}G)' + '(GAAAC){e<2}')
regex_rv = re.compile('(GTTTC){e<2}' + '(C[ATGC]{3}C[ATGC]{3}C[ATGC]{3}C[ATGC]{3})' + '([ATGC]{4})' + 'AA')

# cell line IDs; 
# minumum Hamming distance between any two of them is 2
# keeping TGCT but the cell line is not there
clIDs = ['AACC', 'AAGG', 'ACAC', 'ACCT', 'ACGA', 'ACTG', 'AGAG', 'AGCA', 'AGGT', 'AGTC', 'ATCG', 'ATGC', 'CAAC', 'CACT',
         'CAGA', 'CATG', 'CCAA', 'CCTT', 'CGAT', 'CGTA', 'CTGT', 'CTTC', 'GAAG', 'GATC', 'GCAT', 'GCTA', 'GGAA', 'GGTT',
         'GTAC', 'GTGA', 'GTTG', 'TCCA', 'TGAC', 'TGCT', 'TGTG', 'TTCC', 'TTGG']

# input files w/ raw reads
reads_r1 = open(sample_path + sample + '_R1.fq', 'r')
reads_r2 = open(sample_path + sample + '_R2.fq', 'r')

# output files
clID_bc_out = open(sample_path + sample + '_clID_bc_extracted.txt', 'w+')
failed_clIDs_out = open(sample_path + sample + '_failed_clIDs.txt', 'w+')

# dict. w/ clIDs as keys and bartender input as values
clID_bc_dict = dict(zip(clIDs, [[] for i in range(len(clIDs))]))

i = 0 # counting reads
m = 0 # counting reads that had a regex match and hig enough q-score but did not pair well
o = 0 # counting reads that passed qc but did not have a clID match
n = 0 # counting reads that passed qc and had a clID match

while True:
    
    # reading one .fq entry at a time
    # r1 reads:
    line_0_r1 = reads_r1.readline() 
    if not line_0_r1: # stopping if no more reads
        break
        
    line_1_r1 = reads_r1.readline() # has the sequence
    line_2_r1 = reads_r1.readline()
    line_3_r1 = reads_r1.readline() # has the q-score
    
    # r2 reads:
    line_0_r2 = reads_r2.readline() 
    if not line_0_r2: 
        break
        
    line_1_r2 = reads_r2.readline() 
    line_2_r2 = reads_r2.readline()
    line_3_r2 = reads_r2.readline() 

    # assuming that r1 read is forward
    if regex_fw.search(line_1_r1) and regex_rv.search(line_1_r2): # if the pattern can be matched in both reads
        match_fw = regex_fw.search(line_1_r1) # regular expression search match
        match_rv = regex_rv.search(line_1_r2)
        
        start_r1, end_r1 = match_fw.span(2) 
        avg_q_r1 = avg_qscore(line_3_r1[start_r1 : end_r1]) # avg. q-score in forward read
        
        start_r2, end_r2 = match_rv.span(2) 
        avg_q_r2 = avg_qscore(line_3_r2[start_r2 : end_r2]) # avg. q-score in reverse read
        
        # if avg. q-score is at least 30 
        if avg_q_r1 >= 30 and avg_q_r2 >= 30:
        
            # if clID and barcode match between reads
            if match_fw.group(1) + match_fw.group(2) == rv_cmp(match_rv.group(2) + match_rv.group(3)):
        
                clID = match_fw.group(1) # cell line ID        
                bc = match_fw.group(2) # barcode
                            
                # and cell line ID matches one of those in the list
                if clID in clIDs:
                    n += 1
                
                    # .fastq entry number, cell line ID, barcode
                    clID_bc_out.write('\t'.join([str(i + 1), clID, bc]) + '\n') 
                
                    # clID: 'bc,n\n'
                    clID_bc_dict[clID].append('{},{}\n'.format(bc, str(i + 1)))

                else:
                    o += 1
                    
                    # .fastq entry number, cell line ID, barcode
                    failed_clIDs_out.write('\t'.join([str(i // 4 + 1), clID, bc]) + '\n')
            
            else:
                m += 1
                
                
    # assuming that r2 read is forward        
    elif regex_fw.search(line_1_r2) and regex_rv.search(line_1_r1):
        match_fw = regex_fw.search(line_1_r2)
        match_rv = regex_rv.search(line_1_r1)
        
        start_r2, end_r2 = match_fw.span(2) 
        avg_q_r2 = avg_qscore(line_3_r2[start_r2 : end_r2]) 
        
        start_r1, end_r1 = match_rv.span(2) 
        avg_q_r1 = avg_qscore(line_3_r1[start_r1 : end_r1])
        
        if avg_q_r2 >= 30 and avg_q_r1 >= 30:
        
            if match_fw.group(1) + match_fw.group(2) == rv_cmp(match_rv.group(2) + match_rv.group(3)):
        
                clID = match_fw.group(1) 
                bc = match_fw.group(2)
                            
                if clID in clIDs:
                    n += 1
                    clID_bc_out.write('\t'.join([str(i + 1), clID, bc]) + '\n') 
                    clID_bc_dict[clID].append('{},{}\n'.format(bc, str(i + 1)))

                else:
                    o += 1
                    failed_clIDs_out.write('\t'.join([str(i // 4 + 1), clID, bc]) + '\n')
                                
            else:
                m += 1
 
    i += 1

reads_r1.close()
reads_r2.close()
clID_bc_out.close()
failed_clIDs_out.close()

# writing stats out
stats_out = open(sample_path + sample + '_find_barcodes_stats.txt', 'w+')
stats_out.write('Total reads: ' + str(i) + '\n')
stats_out.write('Had a regex match and a high q-score but did not pair well: ' + str(o) + '\n')
stats_out.write('Passed qc but no clID match: ' + str(m) + '\n')
stats_out.write('Passed qc and had a known clID: ' + str(n) + '\n')
stats_out.close()


for i in clIDs:
    if len(clID_bc_dict[i]) == 0:
        continue # writing output only if there are barcodes for this cell line
    with open(sample_path + sample + '_' + i + '_bc_extracted.txt', 'w+') as o:
        o.writelines(clID_bc_dict[i])