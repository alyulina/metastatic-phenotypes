import regex as re
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

def avg_qscore(x):
    return np.mean([ord(i) - 33 for i in x])

# regular expression pattern for the barcode region:
# regex = re.compile('TT' + '(clID)' + '(barcode)' + '(GAAAC){e<2}'), 
# allowing at most one error in the subsequent region and no errors in the preceeding region 
regex = re.compile('TT' + '(....)' + '(...G...G...G...G)' + '(GAAAC){e<2}')

# cell line IDs; 
# minumum Hamming distance between any two of them is 2
clIDs = ['AACC', 'AAGG', 'ACAC', 'ACCT', 'ACGA', 'ACTG', 'AGAG', 'AGCA', 'AGGT', 'AGTC', 'ATCG', 'ATGC', 'CAAC', 'CACT',
         'CAGA', 'CATG', 'CCAA', 'CCTT', 'CGAT', 'CGTA', 'CTGT', 'CTTC', 'GAAG', 'GATC', 'GCAT', 'GCTA', 'GGAA', 'GGTT',
         'GTAC', 'GTGA', 'GTTG', 'TCCA', 'TGAC', 'TGCT', 'TGTG', 'TTCC', 'TTGG']

with open(sample_path + sample + '_merged.assembled.fastq', 'r') as f:
    lines = f.readlines()

clID_bc_out = open(sample_path + sample + '_clID_bc_extracted.txt', 'w+')

# dict. w/ clIDs as keys and bartender input as values
clID_bc_dict = dict(zip(clIDs, [[] for i in range(len(clIDs))]))

for i in range(len(lines)):
    if i % 4 == 0:  # for each merged read

        if regex.search(lines[i + 1]):  # if the pattern can be matched
            match = regex.search(lines[i + 1])  # regular expression search match
            clID = match.group(1)  # cell line ID
            bc = match.group(2)  # barcode
            start, end = match.span(2)
            avg_q = avg_qscore(lines[i + 3][start: end])  # avg. q-score

            # if avg. q-score is at least 30 and cell line ID matches one of those in the list
            if set(bc).issubset({'A', 'T', 'G', 'C', 'N'}) == True and avg_q >= 30 and clID in clIDs:
                # .fastq entry number, cell line ID, barcode, avg. q-score, and first .fastq entry line
                clID_bc_out.write('\t'.join([str(i // 4 + 1), clID, bc, str(avg_q), lines[i].strip('\n')]) + '\n')

                # clID: 'bc,n\n'
                clID_bc_dict[clID].append('{},{}\n'.format(bc, str(i // 4 + 1)))

print('Fraction of mapped reads: ' + str(len(found) / (len(lines) // 4)))  # how many reads had a barcode match

for i in clIDs:
    if len(clID_bc_dict[i]) == 0:
        continue # writing output only if there are barcodes for this cell line
    with open(sample_path + sample + '_' + i + '_bc_extracted.txt', 'w+') as o:
        o.writelines(clID_bc_dict[i])