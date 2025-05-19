# work in progress!

import os
import re
import math
import random
import pandas as pd
import numpy as np
import scipy as sp
from scipy import stats
import itertools
from itertools import combinations
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import cm
import matplotlib.font_manager

def versions():
    #print('./utils.py')
    for name, module in list(globals().items()):
        if hasattr(module, '__version__'):
            print(f'{name}={module.__version__}')


# setting up matplotlib:
font = matplotlib.font_manager.FontProperties(fname='/System/Library/Fonts/Helvetica.ttc') # font (although this does not matted because of what's next)
matplotlib.rcParams['pdf.fonttype'] = 42; matplotlib.rcParams['ps.fonttype'] = 42 # these two lines make text editable later if needed
matplotlib.rcParams['axes.facecolor'] = 'none'; matplotlib.rcParams['figure.facecolor'] = 'none' # removing figure background color

# where to find sample information & sequencing data
# excluded samples: SK1_8 (very low coverage), SK1_11 (duplicate of SK1_10), LS01_17 (no spike-in detected), LS01_32 (no reads)
metadata = pd.read_csv('/Users/alyulina/Research/metastatic-phenotypes/data/metadata.csv', index_col='sample id') # all sequencing samples
read_counts_by_cell_line = pd.read_csv('/Users/alyulina/Research/metastatic-phenotypes/data/cell-line_read_counts.csv', index_col='clID')
cell_counts_by_cell_line = pd.read_csv('/Users/alyulina/Research/metastatic-phenotypes/data/cell-line_cell_counts.csv', index_col='clID')
path_to_reads = '/Users/alyulina/Research/metastatic-phenotypes/data/barcode_counts/'
path_to_bootstraps = '/Users/alyulina/Research/metastatic-phenotypes/data/bootstraps/'

# where to save plots
path_to_save_figs = '/Users/alyulina/Research/metastatic-phenotypes/figures/'

# cell line IDs & corresponding cell line labels, spike-in must be first
# excluded: other spike-ins & mt4-2D, which had the same clIDs
clIDs = ['GATC', 'AAGG', 'ACAC', 'ACCT', 'ACGA', 'ACTG', 'AGAG', 'AGCA', 'AGGT', 'AGTC', 'ATCG', 'ATGC', 'CAAC', 'CACT', 'CAGA', 'CATG', 'CCAA', 'CCTT', 'CGAT', 'CGTA', 'CTGT', 'CTTC', 'GAAG', 'GCAT', 'GCTA', 'GGAA', 'GGTT', 'GTAC', 'GTGA', 'GTTG', 'TCCA', 'TGAC', 'TTCC', 'TTGG']
cell_line_labels = ['Spike-in 1', '7160c2', 'FC1199', '6694c2', '7160c5', 'BF857', 'BF1987', 'Panc02', '0688M', 'FC1245', 'BF4326', 'BF2117', '6419c5', 'BF4326Rep', 'KPC960', '0755P', 'BF1836', 'BF2014', 'BF5960', '6499c4', '6422c5', 'FC1245Rep', 'KPC-JH', 'KPC960Rep', 'BF2153', 'KC6141', 'mT3-2D', 'mT5-2D', '0764P', '6421c2', '0755A', '2838c3', 'FC1242', 'KPC961']
clID__label = dict(zip(clIDs, cell_line_labels))
clIDs_no_spikein = clIDs[1:]
clIDs_not_rejected = [x for x in clIDs_no_spikein if clID__label[x] not in ['KPC960', 'KPC960Rep', 'KPC961']] # and no spike-in
clIDs_rep = [['AAGG', 'ACGA'], ['AGTC', 'CTTC'], ['ATCG', 'CACT'], ['CAGA', 'GCAT']]

# colors
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('custom_cmap', ['#c13e2e', '#e86b5a', '#f9a28d', '#e9e0c8', '#a3c3c3', '#5e9faf', '#00768d'][::-1])
colors = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=1), cmap=cmap)

def map_colors(clID__value):
    return {clID: colors.to_rgba(np.log10(value)) for clID, value in clID__value.items()}

# colors based on 3w F1 lung expansion from experiment one
clID__color_in_vivo_exp1 = {'AAGG': (0.6137485582468281, 0.7514186851211072, 0.7573241061130334, 1.0),
                            'ACAC': (0.6392156862745098, 0.7647058823529411, 0.7647058823529411, 1.0),
                            'ACCT': (0.7928489042675894, 0.28465974625144175, 0.22099192618223756, 1.0),
                            'ACGA': (0.6779700115340254, 0.7807612456747405, 0.7674740484429066, 1.0),
                            'ACTG': (0.7813148788927335, 0.8235755478662052, 0.7748558246828142, 1.0),
                            'AGAG': (0.9325490196078431, 0.493194925028835, 0.42117647058823526, 1.0),
                            'AGCA': (0.0, 0.4627450980392157, 0.5529411764705883, 1.0),
                            'AGGT': (0.9617070357554787, 0.6925028835063437, 0.6073817762399077, 1.0),
                            'AGTC': (0.9592156862745098, 0.5794694348327566, 0.5011764705882353, 1.0),
                            'ATCG': (0.2862283737024221, 0.5875893886966551, 0.6564705882352941, 1.0),
                            'ATGC': (0.5755478662053056, 0.7314878892733564, 0.7462514417531718, 1.0),
                            'CAAC': (0.956078431372549, 0.5693194925028835, 0.49176470588235294, 1.0),
                            'CACT': (0.40364475201845446, 0.6417993079584775, 0.6964244521337947, 1.0),
                            'CAGA': (0.02602076124567474, 0.47409457900807384, 0.5623529411764706, 1.0),
                            'CATG': (0.7490196078431373, 0.8101960784313725, 0.7725490196078431, 1.0),
                            'CCAA': (0.3845444059976932, 0.6318339100346021, 0.6908881199538639, 1.0),
                            'CCTT': (0.9469434832756631, 0.7497116493656286, 0.6618223760092272, 1.0),
                            'CGAT': (0.7038062283737023, 0.7914648212226066, 0.7693194925028835, 1.0),
                            'CGTA': (0.7296424452133794, 0.8021683967704729, 0.7711649365628604, 1.0),
                            'CTGT': (0.9701960784313726, 0.6149942329873125, 0.5341176470588236, 1.0),
                            'CTTC': (0.940392156862745, 0.5185697808535178, 0.4447058823529411, 1.0),
                            'GAAG': (0.8329873125720877, 0.8449826989619377, 0.7785467128027681, 1.0),
                            'GCAT': (0.11275663206459056, 0.5119261822376009, 0.5937254901960785, 1.0),
                            'GCTA': (0.21683967704728951, 0.5573241061130334, 0.6313725490196078, 1.0),
                            'GGAA': (0.7296424452133794, 0.8021683967704729, 0.7711649365628604, 1.0),
                            'GGTT': (0.7102652825836216, 0.7941407151095732, 0.7697808535178777, 1.0),
                            'GTAC': (0.6844290657439446, 0.783437139561707, 0.7679354094579007, 1.0),
                            'GTGA': (0.8975778546712803, 0.8717416378316032, 0.7831603229527104, 1.0),
                            'GTTG': (0.7102652825836216, 0.7941407151095732, 0.7697808535178777, 1.0),
                            'TCCA': (0.7231833910034601, 0.7994925028835063, 0.7707035755478662, 1.0),
                            'TGAC': (0.891118800461361, 0.8690657439446366, 0.7826989619377163, 1.0),
                            'TTCC': (0.6973471741637831, 0.7887889273356401, 0.7688581314878893, 1.0),
                            'TTGG': (0.13877739331026528, 0.523275663206459, 0.6031372549019608, 1.0)}

# colors based on 3w F1 lung expansion from experiment two
clID__color_in_vivo_exp2 = {'AAGG': (0.645674740484429, 0.7673817762399077, 0.7651672433679354, 1.0),
                            'ACAC': (0.6201153402537485, 0.7547404844290657, 0.7591695501730104, 1.0),
                            'ACCT': (0.9137254901960784, 0.4322952710495962, 0.3647058823529411, 1.0),
                            'ACGA': (0.7877739331026528, 0.8262514417531718, 0.7753171856978085, 1.0),
                            'ACTG': (0.7877739331026528, 0.8262514417531718, 0.7753171856978085, 1.0),
                            'AGAG': (0.956078431372549, 0.5693194925028835, 0.49176470588235294, 1.0),
                            'AGCA': (0.0, 0.4627450980392157, 0.5529411764705883, 1.0),
                            'AGGT': (0.9607843137254902, 0.5845444059976933, 0.5058823529411766, 1.0),
                            'AGTC': (0.9262745098039216, 0.4728950403690887, 0.40235294117647047, 1.0),
                            'ATCG': (0.12143021914648212, 0.5157093425605537, 0.5968627450980393, 1.0),
                            'ATGC': (0.473679354094579, 0.6783391003460207, 0.716724336793541, 1.0),
                            'CAAC': (0.9764705882352941, 0.6352941176470588, 0.5529411764705883, 1.0),
                            'CACT': (0.0, 0.4627450980392157, 0.5529411764705883, 1.0),
                            'CAGA': (0.0, 0.4627450980392157, 0.5529411764705883, 1.0),
                            'CATG': (0.7619377162629757, 0.8155478662053056, 0.7734717416378316, 1.0),
                            'CCAA': (0.30357554786620533, 0.5951557093425606, 0.6627450980392157, 1.0),
                            'CCTT': (0.9529411764705882, 0.5591695501730103, 0.48235294117647054, 1.0),
                            'CGAT': (0.581914648212226, 0.7348096885813148, 0.7480968858131487, 1.0),
                            'CGTA': (0.8459054209919261, 0.8503344867358708, 0.7794694348327567, 1.0),
                            'CTGT': (0.960230680507497, 0.6982237600922722, 0.6128258362168397, 1.0),
                            'CTTC': (0.9262745098039216, 0.4728950403690887, 0.40235294117647047, 1.0),
                            'GAAG': (0.9277508650519031, 0.8240830449826989, 0.7325951557093425, 1.0),
                            'GCAT': (0.0, 0.4627450980392157, 0.5529411764705883, 1.0),
                            'GCTA': (0.44184544405997694, 0.6617301038062283, 0.7074971164936563, 1.0),
                            'GGAA': (0.8782006920415224, 0.8637139561707036, 0.7817762399077278, 1.0),
                            'GGTT': (0.9218454440599769, 0.8469665513264129, 0.7543713956170703, 1.0),
                            'GTAC': (0.7167243367935409, 0.7968166089965397, 0.7702422145328719, 1.0),
                            'GTGA': (0.7877739331026528, 0.8262514417531718, 0.7753171856978085, 1.0),
                            'GTTG': (0.9307035755478662, 0.812641291810842, 0.7217070357554787, 1.0),
                            'TCCA': (0.6392156862745098, 0.7647058823529411, 0.7647058823529411, 1.0),
                            'TGAC': (0.7038062283737023, 0.7914648212226066, 0.7693194925028835, 1.0),
                            'TTCC': (0.7361014994232987, 0.8048442906574395, 0.7716262975778546, 1.0),
                            'TTGG': (0.02602076124567474, 0.47409457900807384, 0.5623529411764706, 1.0)}

# highlighting replicate or related cell lines in plots
clID_replicates__linestyle = {('AAGG', 'ACGA'): ':', ('AGTC', 'CTTC'): '-', ('ATCG', 'CACT'): '-', ('CAGA', 'GCAT'): '-'}

# now, some helpful functions:
def convert_barcode_reads_to_cell_counts(merged_sample_ids,
                                         cell_line_ids=clIDs,
                                         path_to_data=path_to_reads,
                                         read_counts=read_counts_by_cell_line,
                                         metadata=metadata):
    """
    Converts barcode reads to estimated cell counts
    based on spiked in cell line read count.

    Args:
        merged_sample_ids (list): List of samples; samples that should be merged
                           must be supplied comma-separated as listed in the metadata 'notes' column.
        cell_line_ids (list): Four-nucleotide cell line IDs, optional.
        path_to_data (str): Path to the folder containing files with read counts, optional (path_to_reads defined above by default).
        read_counts (pandas.DataFrame): A precomputed table with read counts, optional (read_counts_by_cell_line.csv by default).
        metadata (pandas.DataFrame): A table containing experimental details, including initial cell numbers, optional (metadata.csv by default).

    Returns:
        tumor_sizes (pandas.DataFrame): A table with 'mouse', 'clID', 'barcode', 'size' (in cells) columns.
    """
    sample__counts = {}

    for sample in merged_sample_ids:
        # when we need to merge reads first before converting to counts
        clID_bc__n_reads = {}; R = 0; one = False

        sub_samples = sample.split(':')

        for sub_sample in sub_samples:
            file_path = os.path.join(path_to_data, f'{sub_sample}_clIDs_rBC_cluster_counts.txt')
            if not os.path.exists(file_path):
                print(f'File {file_path} not found.')
                continue

            if metadata.loc[sub_sample, 'spike-in added'] == False:
                print(f'No spike-in added in {sub_sample} ({metadata.loc[sub_sample, 'experiment']}, {metadata.loc[sub_sample, 'injection method']} {metadata.loc[sub_sample, 'tissue']}).')
                continue

            # checking if sample might need to be merged
            if len(sub_samples) == 1: # nothing to merge, will convert reads to counts using spike-in, add to out
                sample__counts[sub_sample] = dict(zip(clIDs, [{} for i in clIDs]))

                with open(file_path, 'r') as f:
                    for line in f:
                        clID_barcode, read_count = line.strip().split('\t')
                        clID, barcode = clID_barcode.split('_'); read_count = int(read_count)
                        if clID not in clIDs: # this should not happen if barcodes were processed correctly
                            continue
                        cell_count = 50000 / (1 + read_counts.loc['GATC', sub_sample]) * read_count
                        sample__counts[sub_sample][clID][barcode] = cell_count # adding cell counts!

            else:
                spikein = metadata.loc[sub_sample, 'notes'].split()[-2]
                if spikein == 'separate': # go in loop and add cell counts

                    if sample not in sample__counts:
                        sample__counts[sample] = dict(zip(clIDs, [{} for i in clIDs]))

                    with open(file_path, 'r') as f:
                        for line in f:
                            clID_barcode, read_count = line.strip().split('\t')
                            clID, barcode = clID_barcode.split('_'); read_count = int(read_count)
                            if clID not in clIDs:
                                continue
                            cell_count = 50000 / (1 + read_counts.loc['GATC', sub_sample]) * read_count

                            if barcode in sample__counts[sample][clID]:
                                sample__counts[sample][clID][barcode] += cell_count # adding cell counts!
                            else:
                                sample__counts[sample][clID][barcode] = cell_count

                elif spikein == 'one': # go in loop and merge reads first, then convert to counts

                    one = True
                    R += read_counts.loc['GATC', sub_sample] # retrieving spike-in reads

                    with open(file_path, 'r') as f:
                        for line in f:
                            clID_barcode, read_count = line.strip().split('\t')

                            if clID_barcode in clID_bc__n_reads:
                                clID_bc__n_reads[clID_barcode] += int(read_count)
                            else:
                                clID_bc__n_reads[clID_barcode] = int(read_count)

        if one:
            sample__counts[sample] = dict(zip(clIDs, [{} for i in clIDs]))
            for clID_bc, n_reads in clID_bc__n_reads.items():
                clID, barcode = clID_bc.split('_')
                if clID not in clIDs: # this should not happen if barcodes were processed correctly
                    continue
                sample__counts[sample][clID][barcode] = 50000 / (1 + R) * n_reads # 50,000 cells were spiked-in

    records = []
    for mouse, clID__counts in sample__counts.items():
        for clID, barcode__counts in clID__counts.items():
            for barcode, size in barcode__counts.items():
                records.append((mouse, clID, barcode, size))

    df = pd.DataFrame(records, columns=['mouse', 'clID', 'barcode', 'size'])

    return df.loc[df['clID'].isin(cell_line_ids)]


def get_burden_and_n(merged_sample_ids,
                     cell_line_ids=clIDs,
                     min_tumor_size=100,
                     path_to_data=path_to_reads,
                     read_counts=read_counts_by_cell_line,
                     cell_counts=cell_counts_by_cell_line,
                     metadata=metadata):
    """
    Calculates the two statistics on a set of samples: relative burden and number of tumors.

    Args:
        merged_sample_ids (list): Samples to include in the analysis;
                                  merged samples need to be colon-separated.
        cell_line_ids (list): Four-nucleotide cell line IDs, optional. Will use all cell lines except for spike-in.
        min_tumor_size (int): Minimum tumor size, cells, optional (default=100).
        path_to_data (str): Path to the folder containing files with read counts, optional (path_to_reads defined above by default).
        read_counts (pandas.DataFrame): A precomputed table with read counts, optional (read_counts_by_cell_line.csv by default).
        cell_counts (pandas.DataFrame): A precomputed table with cell counts, optional (cell_counts_by_cell_line.csv by default).
        metadata (pandas.DataFrame): A table containing experimental details, including initial cell numbers, optional (metadata.csv by default).

    Returns:
        clID__relative_burden (dict): Relative burden per cell line.
        clID__n_tumors (dict): Number of tumors per cell line.
    """

    clID__relative_burden = {}
    clID__n_tumors = {}

    # just making sure that the spike-in is not there
    cell_line_ids_no_spikein = [clID for clID in cell_line_ids if clID != 'GATC']
    n_mice = len(merged_sample_ids)

    # precompute pre-injection cell counts for each clID first,
    # by finding their fractions and multiplying by the total number of cells injected
    experiment = metadata.loc[[s.split(':')[0] if ':' in s else s for s in merged_sample_ids], 'experiment'].iloc[0]
    pre_injection_samples = metadata[(metadata['time point, d'] == 0) & (metadata['experiment'] == experiment)]
    clID__n_0 = {clID: np.mean([read_counts[sample_id][clID] / sum(read_counts[sample_id][1:]) * np.mean([metadata.loc[x]['initial number of cells'] for x in sample_id.split(':')]) for sample_id in pre_injection_samples.index.tolist()]) for clID in cell_line_ids_no_spikein} # excluding spike-in

    # gather tumors for each sample + precompute mouse weights
    # just aggregarind all tumor sizes, without the size cutoff
    tumor_sizes = convert_barcode_reads_to_cell_counts(merged_sample_ids, cell_line_ids=cell_line_ids_no_spikein, # keeping all clIDs but spike-in
                                                       path_to_data=path_to_data, read_counts=read_counts_by_cell_line, metadata=metadata)

    clID__burden = {}
    for clID in cell_line_ids_no_spikein:

        clID_tumors = tumor_sizes[tumor_sizes['clID'] == clID] # filter tumors of this cell line
        clID__n_tumors[clID] = clID_tumors[clID_tumors['size'] >= min_tumor_size].shape[0]
        clID__burden[clID] = clID_tumors['size'].sum() / (clID__n_0.get(clID, 1) * n_mice)

    average_burden = np.mean(list(clID__burden.values()))
    clID__relative_burden = {clID: clID__burden[clID] / average_burden for clID in clID__burden.keys()}

    return clID__relative_burden, clID__n_tumors


# let's define one bootstrap function for everything! will bootstrap mice, tumors, weight samples by number of tumors above a certain size
# return: bootstrap distributions for per-cell-line normalized burden relative to average across all cell lines and number of tumors
def bootstrap_burden_n_sizes(merged_sample_ids,
                             cell_line_ids=clIDs,
                             B_mice=100,
                             B_tumors=100,
                             min_tumor_size=100,
                             weighted=True,
                             path_to_data=path_to_reads,
                             read_counts=read_counts_by_cell_line,
                             cell_counts=cell_counts_by_cell_line,
                             metadata=metadata):
    """
    Bootstraps mice and tumors to calculate various statistics.

    Args:
        merged_sample_ids (list): Samples to include in the analysis;
                                  merged samples need to be colon-separated.
        cell_line_ids (list): Four-nucleotide cell line IDs, optional. Will use all cell lines except for spike-in.
        B_mice (int): Number of bootstrap resampling iterations for mice, optional (default=1000).
        B_tumors (int): Number of bootstrap resampling iterations for tumors, optional (default=1000).
        min_tumor_size (int): Minimum tumor size, cells, optional (default=100).
        weighted (bool): Whether to weight bootstrap inerations by the total number tumors larger than 100 cells
                         when generating the bootstrap distribution, optional (default=True).
        path_to_data (str): Path to the folder containing files with read counts, optional (path_to_reads defined above by default).
        read_counts (pandas.DataFrame): A precomputed table with read counts, optional (read_counts_by_cell_line.csv by default).
        cell_counts (pandas.DataFrame): A precomputed table with cell counts, optional (cell_counts_by_cell_line.csv by default).
        metadata (pandas.DataFrame): A table containing experimental details, including initial cell numbers, optional (metadata.csv by default).

    Returns:
        clID__relative_burden_bootstraps (dict): A dictionary containing B_mice x B_tumors bootstrap replicates of the relative tumor burden.
        clID__n_tumors_bootstraps (dict): A dictionary containing B_mice x B_tumors bootstrap replicates of the number of tumors above the minimum size.
        clID__sizes_bootstraps (dict): A dictionary containing B_mice x B_tumors bootstrap replicates of the tumor sizes above the minimum size.
    """

    # used to have / might be nice to have (but right now not necessary):
    # clID__tumors_bootstraps (dict): A dictionary containing B_mice x B_tumors lists with tumor sizes above the minimum size.
    
    # just making sure that the spike-in is not there
    cell_line_ids_no_spikein = [clID for clID in cell_line_ids if clID != 'GATC']
    clID__id = {clID: i for i, clID in enumerate(cell_line_ids_no_spikein)} # let's map clIDs to indices
    id__clID = {i: clID for i, clID in enumerate(cell_line_ids_no_spikein)} # and back
    n_clIDs = len(cell_line_ids_no_spikein)

    # precompute pre-injection cell counts for each clID first,
    # by finding their fractions and multiplying by the total number of cells injected
    experiment = metadata.loc[[s.split(':')[0] if ':' in s else s for s in merged_sample_ids], 'experiment'].iloc[0]
    pre_injection_samples = metadata[(metadata['time point, d'] == 0) & (metadata['experiment'] == experiment)]
    # working w/ clID indices from now on
    clID__n_0 = {clID__id[clID]: np.mean([read_counts[sample_id][clID] / sum(read_counts[sample_id][1:]) * np.mean([metadata.loc[x]['initial number of cells'] for x in sample_id.split(':')]) for sample_id in pre_injection_samples.index.tolist()]) for clID in cell_line_ids_no_spikein} # excluding spike-in
    n0_per_clID = np.array([clID__n_0[i] for i in range(n_clIDs)])  # shape (n_clIDs,)

    # gather tumors for each sample + precompute mouse weights
    # just aggregarind all tumor sizes, without the size cutoff
    tumor_sizes = convert_barcode_reads_to_cell_counts(merged_sample_ids, cell_line_ids=cell_line_ids_no_spikein, # keeping all clIDs but spike-in
                                                       path_to_data=path_to_data, read_counts=read_counts_by_cell_line, metadata=metadata)

    # converting clIDs to indices
    tumor_sizes['clID'] = tumor_sizes['clID'].map(clID__id)

    # group by mouse, precompute per-mouse arrays of [clID, size]
    mouse_tumors = tumor_sizes.groupby('mouse')[['clID', 'size']].apply(lambda x: x.to_numpy())

    # calculating weights for mouse resampling: total number of tumors above the size cutoff per mouse
    # using all cell lines, not just those requested, to get weights
    mouse__n_tumors = tumor_sizes[(tumor_sizes['size'] >= min_tumor_size)].groupby('mouse').size().to_dict()
    weights_not_normalized_yet = [mouse__n_tumors.get(mouse, 0) + 1 if weighted else 1 for mouse in merged_sample_ids]
    weights = np.array(weights_not_normalized_yet) / np.sum(weights_not_normalized_yet)
    n_mice = len(merged_sample_ids) # number of mice

    clID__relative_burden_bootstraps = {clID: [] for clID in cell_line_ids_no_spikein}
    clID__n_tumors_bootstraps = {clID: [] for clID in cell_line_ids_no_spikein}
    clID__sizes_bootstraps = {clID: [] for clID in cell_line_ids_no_spikein}
    clID__sizes_stats_bootstraps = {clID: [] for clID in cell_line_ids_no_spikein}
    for _ in range(B_mice):
        sampled_mice = np.random.choice(merged_sample_ids, size=n_mice, replace=True, p=weights)

        # now, want to get all tumors from the sampled mice & keep the clID info
        # initialize an empty list to collect the tumor sizes along with their clID information for the sampled mice
        sampled_sizes_with_clID = np.concatenate([mouse_tumors[mouse] for mouse in sampled_mice])
        sampled_sizes_with_clID[:, 1] = sampled_sizes_with_clID[:, 1].astype(float)
        n_tumors = len(sampled_sizes_with_clID) # all tumors, across all cell lines

        # resample indices once: shape (B_tumors, n_tumors)
        resampled_indices = np.random.randint(0, n_tumors, size=(B_tumors, n_tumors))

        # create resampled arrays: shape (B_tumors, n_tumors, 2) & split columns into two arrays
        resampled_sizes_with_clID = sampled_sizes_with_clID[resampled_indices]
        clIDs_resampled = resampled_sizes_with_clID[:, :, 0].astype(int)  # clID indices: shape (B_tumors, n_tumors)
        sizes_resampled = resampled_sizes_with_clID[:, :, 1].astype(float)  # tumor sizes: shape (B_tumors, n_tumors)

        
        # collect tumor sizes above threshold per clID per bootstrap
        for b in range(B_tumors):
            clIDs_b = clIDs_resampled[b]
            sizes_b = sizes_resampled[b]
            for i, clID in id__clID.items():
                mask = (clIDs_b == i) & (sizes_b >= min_tumor_size)
                clID__sizes_bootstraps[clID] += list(sizes_b[mask])
        
        # we want to compute two things for each of the B_tumors samples:
        # sum of tumor sizes per clID and count of tumors ≥ minimum size per clID
        # we'll be clever and count bins!
        sum_per_clID = np.array([np.bincount(clIDs_resampled[b], weights=sizes_resampled[b], minlength=n_clIDs) for b in range(B_tumors)])
        count_per_clID = np.array([np.bincount(clIDs_resampled[b][sizes_resampled[b] >= min_tumor_size], minlength=n_clIDs) for b in range(B_tumors)])

        burdens = sum_per_clID / (n0_per_clID[None, :] * n_mice) # adding an extra dimension to an array for broadcasting
        relative_burdens = burdens / burdens.mean(axis=1, keepdims=True)

        for i, clID in id__clID.items():
            clID__relative_burden_bootstraps[clID] += list(relative_burdens[:, i])
            clID__n_tumors_bootstraps[clID] += list(count_per_clID[:, i])

    for clID in cell_line_ids_no_spikein:
        sizes = clID__sizes_bootstraps[clID]
        if sizes:
            clID__sizes_stats_bootstraps[clID] = {
                2.5: np.quantile(sizes, 0.025),
                50: np.median(sizes),
                97.5: np.quantile(sizes, 0.975)
            }
        else:
            clID__sizes_stats_bootstraps[clID] = {2.5: np.nan, 50: np.nan, 97.5: np.nan}

    return clID__relative_burden_bootstraps, clID__n_tumors_bootstraps, clID__sizes_stats_bootstraps


# getting the bootstrap distribution of y/x statistics
def bootstrap_index(clID__relative_burden_bootstraps_x, clID__relative_burden_bootstraps_y, epsilon=1e-6):
    out = {}
    for clID in clID__relative_burden_bootstraps_x.keys():
        numerator = np.array(clID__relative_burden_bootstraps_y[clID])
        denominator = np.array(clID__relative_burden_bootstraps_x[clID])

        # check for zeros before division and print the clID
        if np.any(denominator == 0):
            print(f"Warning: clID {clID} has a zero denominator.")

        out[clID] = (numerator + epsilon * (denominator == 0)) / (denominator + epsilon * (denominator == 0))  # add epsilon only for zeros
    return out


# functions to make the many plots:
def correlation_plot(cell_line_ids,
                     clID__x, clID__xerr,
                     clID__y, clID__yerr,
                     clID__color=clID__color_in_vivo_exp1, # do not forget to switch to the right experiment
                     clID__label=clID__label, # labels
                     label_immune_rejected=True, # whether to label immune-rejected cell lines
                     xtitle='metric', ytitle='metric', # metrics to be displayed as axis labels
                     xlims=[1e-3, 1e3], xticks=[1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3], xlabels=['$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$1$', '$10^{1}$', '$10^{2}$', '$10^{3}$'],
                     ylims=[1e-3, 1e3], yticks=[1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3], ylabels=['$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$1$', '$10^{1}$', '$10^{2}$', '$10^{3}$']):
    """
    Makes a scatter plot comparing two experiments or conditions,
    and calculates Spearman and Pearson correlation coefficients with p-values.

    Makes a bar plot with some kind of metric with on a log-scale y-axis, and cell lines on the x-axis.

    Parameters:
        cell_line_ids (list): Cell line IDs, in any order; all of these cell lines will be used to compute correlation coefficients.
        clID__x, clID__y (dict): Dictionaries with cell line IDs as keys and mean metric values as values (e.g. 'AAGT': 0.02).
        clID__xerr, clID__yerr (dict): Dictionaries with cell line IDs as keys
                                       and -/+ errors from the mean as values (e.g. 'AAGT': [0.01, 0.01]).
        clID__color (dict): Dictionary mapping cell line IDs to colors (by default, the one defined above based on the first series of experiments).
        clID__label (dict): Dictionary mapping cell line IDs to label names (by default, the one defined above).
        label_immune_rejected (bool): Whether to label immune-rejected cell lines (on by default).
        xtitle, ytitle (string): Labels for the axes.
        xlims, ylims (list): Limits for the axes (default lims=[1e-3, 1e3]).
        xticks, yticks (list): Locations for axes ticks (default ticks=[1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]).
        xlabels, ylabels (list): Locations for axes tick labels (default labels=['1e-3', '1e-2', '1e-1', '1e0', '1e1', '1e2', '1e3']).

    Returns:
        ax (matplotlib.axes): A plot!
    """

    cell_line_ids_no_spikein = [clID for clID in cell_line_ids if clID != 'GATC']

    fig, ax = plt.subplots(figsize=(4, 4))

    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['bottom', 'left']].set_linewidth(0.5)
    ax.tick_params(width=0.5)

    xs = [clID__x[x] for x in cell_line_ids_no_spikein]
    ys = [clID__y[y] for y in cell_line_ids_no_spikein]

    ax.scatter(xs, ys, color=[clID__color[z] for z in cell_line_ids_no_spikein], s=60, edgecolor='black', linewidth=0.5, alpha=1, zorder=1)
    ax.errorbar(xs, ys, xerr=[[clID__xerr[x][0] for x in cell_line_ids_no_spikein], [clID__xerr[x][1] for x in cell_line_ids_no_spikein]], yerr=[[clID__yerr[y][0] for y in cell_line_ids_no_spikein], [clID__yerr[y][1] for y in cell_line_ids_no_spikein]], color='black', linewidth=0.5, ls='none', zorder=-2)

    ax.plot([min(xlims[0], ylims[0]), max(xlims[1], ylims[1])], [min(xlims[0], ylims[0]), max(xlims[1], ylims[1])], color='#c1c1c1', linewidth=0.5, zorder=-3)  # diagonal line
    ax.plot(xlims, [1e0, 1e0], color='#c1c1c1', linewidth=0.5, alpha=0.3, zorder=-3)  # horizontal line
    ax.plot([1e0, 1e0], ylims, color='#c1c1c1', linewidth=0.5, alpha=0.3, zorder=-3)  # vertical line

    spearmans_rho, spearmans_pval = sp.stats.spearmanr(xs, ys)  # rank correlation coefficient and p-value
    pearsons_r, pearsons_pval = sp.stats.pearsonr(xs, ys)  # Pearson correlation coefficient and p-value
    ax.text(2 * xlims[0], 0.3 * ylims[1], f'ρ = {spearmans_rho:.2f}, p = {spearmans_pval:.0e}', size=10)
    ax.text(2 * xlims[0], 0.15 * ylims[1], f'r = {pearsons_r:.2f}, p = {pearsons_pval:.0e}', size=10)

    if label_immune_rejected == True:
        for clID in cell_line_ids_no_spikein:
            if clID in clIDs_not_rejected:
                continue
            plt.text(clID__x[clID], clID__y[clID], clID__label[clID]) # labeling immune rejected clIDs 

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.minorticks_off()
    ax.set_xlim(xlims[0], xlims[1])
    ax.set_ylim(ylims[0], ylims[1])
    ax.set_xticks(xticks, xlabels, size=12)
    ax.set_yticks(yticks, ylabels, size=12)

    ax.set_xlabel(xtitle, fontsize=14, labelpad=6)
    ax.set_ylabel(ytitle, fontsize=14, labelpad=6)

    return ax
