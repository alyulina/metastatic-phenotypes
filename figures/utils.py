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
path_to_reads = '/Users/alyulina/Research/metastatic-phenotypes/data/bc_counts/'

# where to save plots
path_to_save_figs = '/Users/alyulina/Research/metastatic-phenotypes/figures/'

# cell line IDs, spike-in must be first
# excluded: other spike-ins & mt4-2D, which had the same clIDs
clIDs = ['GATC', 'AAGG', 'ACAC', 'ACCT', 'ACGA', 'ACTG', 'AGAG', 'AGCA', 'AGGT', 'AGTC', 'ATCG', 'ATGC', 'CAAC', 'CACT', 'CAGA', 'CATG', 'CCAA', 'CCTT', 'CGAT', 'CGTA', 'CTGT', 'CTTC', 'GAAG', 'GCAT', 'GCTA', 'GGAA', 'GGTT', 'GTAC', 'GTGA', 'GTTG', 'TCCA', 'TGAC', 'TTCC', 'TTGG']

# corresponding cell line labels
cell_line_labels = ['Spike-in 1', '7160c2', 'FC1199', '6694c2', '7160c5', 'BF857', 'BF1987', 'Panc02', '0688M', 'FC1245', 'BF4326', 'BF2117', '6419c5', 'BF4326Rep', 'KPC960', '0755P', 'BF1836', 'BF2014', 'BF5960', '6499c4', '6422c5', 'FC1245Rep', 'KPC-JH', 'KPC960Rep', 'BF2153', 'KC6141', 'mT3-2D', 'mT5-2D', '0764P', '6421c2', '0755A', '2838c3', 'FC1242', 'KPC961']

clID__label = dict(zip(clIDs, cell_line_labels))

# colors based on 3w F1 lung expansion from experiment one
clID__color_in_vivo_exp1 = {'ACCT': (0.7568627450980392, 0.24313725490196078, 0.1803921568627451, 1.0),
                            'AGAG': (0.9341176470588235, 0.4982698961937715, 0.4258823529411764, 1.0),
                            'AGTC': (0.9388235294117646, 0.5134948096885812, 0.43999999999999995, 1.0),
                            'CTGT': (0.9654901960784313, 0.5997693194925028, 0.52, 1.0),
                            'CTTC': (0.9670588235294117, 0.6048442906574394, 0.5247058823529412, 1.0),
                            'CAAC': (0.9670588235294117, 0.6048442906574394, 0.5247058823529412, 1.0),
                            'TGAC': (0.9454671280276816, 0.7554325259515571, 0.6672664359861592, 1.0),
                            'CCTT': (0.9336562860438292, 0.801199538638985, 0.7108189158016147, 1.0),
                            'AGGT': (0.8975778546712803, 0.8717416378316032, 0.7831603229527104, 1.0),
                            'GTGA': (0.865282583621684, 0.8583621683967705, 0.7808535178777393, 1.0),
                            'GAAG': (0.8071510957324106, 0.8342791234140715, 0.7767012687427912, 1.0),
                            'ACTG': (0.7554786620530565, 0.8128719723183391, 0.7730103806228373, 1.0),
                            'CGTA': (0.7554786620530565, 0.8128719723183391, 0.7730103806228373, 1.0),
                            'GTTG': (0.7296424452133794, 0.8021683967704729, 0.7711649365628604, 1.0),
                            'CATG': (0.6908881199538639, 0.7861130334486736, 0.768396770472895, 1.0),
                            'TTCC': (0.6521337946943483, 0.7700576701268742, 0.7656286043829296, 1.0),
                            'GGTT': (0.645674740484429, 0.7673817762399077, 0.7651672433679354, 1.0),
                            'ACGA': (0.6328489042675893, 0.7613840830449826, 0.7628604382929642, 1.0),
                            'GGAA': (0.6328489042675893, 0.7613840830449826, 0.7628604382929642, 1.0),
                            'TCCA': (0.6264821222606689, 0.7580622837370241, 0.7610149942329872, 1.0),
                            'AAGG': (0.6201153402537485, 0.7547404844290657, 0.7591695501730104, 1.0),
                            'CGAT': (0.5946482122260669, 0.7414532871972318, 0.7517877739331026, 1.0),
                            'GTAC': (0.5882814302191465, 0.7381314878892733, 0.7499423298731257, 1.0),
                            'ATGC': (0.5691810841983851, 0.7281660899653979, 0.7444059976931948, 1.0),
                            'ACAC': (0.5373471741637832, 0.7115570934256055, 0.7351787773933103, 1.0),
                            'CCAA': (0.29490196078431374, 0.5913725490196079, 0.659607843137255, 1.0),
                            'TTGG': (0.24286043829296425, 0.5686735870818915, 0.6407843137254902, 1.0),
                            'GCTA': (0.21683967704728951, 0.5573241061130334, 0.6313725490196078, 1.0),
                            'GCAT': (0.19949250288350634, 0.5497577854671281, 0.6250980392156863, 1.0),
                            'CACT': (0.19949250288350634, 0.5497577854671281, 0.6250980392156863, 1.0),
                            'CAGA': (0.12143021914648212, 0.5157093425605537, 0.5968627450980393, 1.0),
                            'ATCG': (0.06071510957324106, 0.4892272202998847, 0.5749019607843138, 1.0),
                            'AGCA': (0.0, 0.4627450980392157, 0.5529411764705883, 1.0)}

clID__color_in_vivo_exp2 = {'ACCT': (0.7568627450980392, 0.24313725490196078, 0.1803921568627451, 1.0),
                            'AGTC': (0.9341176470588235, 0.4982698961937715, 0.4258823529411764, 1.0),
                            'CTTC': (0.9388235294117646, 0.5134948096885812, 0.43999999999999995, 1.0),
                            'AGAG': (0.9654901960784313, 0.5997693194925028, 0.52, 1.0),
                            'CCTT': (0.9670588235294117, 0.6048442906574394, 0.5247058823529412, 1.0),
                            'AGGT': (0.9670588235294117, 0.6048442906574394, 0.5247058823529412, 1.0),
                            'CTGT': (0.9454671280276816, 0.7554325259515571, 0.6672664359861592, 1.0),
                            'CAAC': (0.9336562860438292, 0.801199538638985, 0.7108189158016147, 1.0),
                            'GGTT': (0.8975778546712803, 0.8717416378316032, 0.7831603229527104, 1.0),
                            'GAAG': (0.865282583621684, 0.8583621683967705, 0.7808535178777393, 1.0),
                            'GTTG': (0.8071510957324106, 0.8342791234140715, 0.7767012687427912, 1.0),
                            'CGTA': (0.7554786620530565, 0.8128719723183391, 0.7730103806228373, 1.0),
                            'GGAA': (0.7554786620530565, 0.8128719723183391, 0.7730103806228373, 1.0),
                            'GTGA': (0.7296424452133794, 0.8021683967704729, 0.7711649365628604, 1.0),
                            'ACTG': (0.6908881199538639, 0.7861130334486736, 0.768396770472895, 1.0),
                            'ACGA': (0.6521337946943483, 0.7700576701268742, 0.7656286043829296, 1.0),
                            'CATG': (0.645674740484429, 0.7673817762399077, 0.7651672433679354, 1.0),
                            'TGAC': (0.6328489042675893, 0.7613840830449826, 0.7628604382929642, 1.0),
                            'TTCC': (0.6328489042675893, 0.7613840830449826, 0.7628604382929642, 1.0),
                            'GTAC': (0.6264821222606689, 0.7580622837370241, 0.7610149942329872, 1.0),
                            'AAGG': (0.6201153402537485, 0.7547404844290657, 0.7591695501730104, 1.0),
                            'TCCA': (0.5946482122260669, 0.7414532871972318, 0.7517877739331026, 1.0),
                            'ACAC': (0.5882814302191465, 0.7381314878892733, 0.7499423298731257, 1.0),
                            'CGAT': (0.5691810841983851, 0.7281660899653979, 0.7444059976931948, 1.0),
                            'GCTA': (0.5373471741637832, 0.7115570934256055, 0.7351787773933103, 1.0),
                            'ATGC': (0.29490196078431374, 0.5913725490196079, 0.659607843137255, 1.0),
                            'CCAA': (0.24286043829296425, 0.5686735870818915, 0.6407843137254902, 1.0),
                            'ATCG': (0.21683967704728951, 0.5573241061130334, 0.6313725490196078, 1.0),
                            'TTGG': (0.19949250288350634, 0.5497577854671281, 0.6250980392156863, 1.0),
                            'CACT': (0.19949250288350634, 0.5497577854671281, 0.6250980392156863, 1.0),
                            'CAGA': (0.12143021914648212, 0.5157093425605537, 0.5968627450980393, 1.0),
                            'GCAT': (0.06071510957324106, 0.4892272202998847, 0.5749019607843138, 1.0),
                            'AGCA': (0.0, 0.4627450980392157, 0.5529411764705883, 1.0)}

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
            file_path = os.path.join(path_to_data, f"{sub_sample}_clIDs_rBC_cluster_counts.txt")
            if not os.path.exists(file_path):
                print(f"File {file_path} not found.")
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

    return df



# let's define one bootstrap function for everything! will bootstrap mice, tumors, weight samples by number of tumors above a certain size
# return: bootstrap distributions for per-cell-line normalized burden relative to average across all cell lines, number of tumors, and their sizes
def bootstrap(merged_sample_ids,
              cell_line_ids=[clID for clID in clIDs if clID != 'GATC'],
              B_mice=1000,
              B_tumors=1000,
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
        cell_line_ids (list): Four-nucleotide cell line IDs, with spike-in excluded, optional.
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
    """

    # used to have / might be nice to have (but right now not necessary):
    # clID__tumors_bootstraps (dict): A dictionary containing B_mice x B_tumors lists with tumor sizes above the minimum size.

    # just making sure that the spike-in is not there
    cell_line_ids_no_spikein = [clID for clID in cell_line_ids if clID != 'GATC']

    # precompute pre-injection cell counts for each clID first,
    # by finding their fractions and multiplying by the total number of cells injected
    experiment = metadata.loc[[s.split(':')[0] if ':' in s else s for s in merged_sample_ids], 'experiment'].iloc[0]
    pre_injection = metadata[(metadata['time point, d'] == 0) & (metadata['experiment'] == experiment)]
    clID__n_0 = {clID: np.mean([read_counts[sample_id][clID] / sum(read_counts[sample_id][1:]) * pre_injection['initial number of cells'].mean()
                                for sample_id in pre_injection.index.tolist()]) for clID in cell_line_ids_no_spikein} # excluding spike-in


    # gather tumors for each sample + precompute mouse weights
    # just aggregarind all tumor sizes, without the size cutoff
    tumor_sizes = convert_barcode_reads_to_cell_counts(merged_sample_ids, cell_line_ids=[clID for clID in clIDs if clID !='GATC'], # keeping all clIDs but spike-in
                                                       path_to_data=path_to_data, read_counts=read_counts_by_cell_line, metadata=metadata)

    # calculating weights for mouse resampling: total number of tumors above the size cutoff per mouse
    # using all cell lines, not just those requested, to get weights
    tumor_sizes_grouped_by_mouse = tumor_sizes.groupby('mouse')
    mouse__n_tumors = tumor_sizes[(tumor_sizes['size'] >= min_tumor_size)].groupby('mouse').size().to_dict()
    weights_not_normalized_yet = [mouse__n_tumors.get(mouse, 0) + 1 if weighted else 1 for mouse in merged_sample_ids]
    weights = np.array(weights_not_normalized_yet) / np.sum(weights_not_normalized_yet)
    n_mice = len(merged_sample_ids) # number of mice

    clID__relative_burden_bootstraps = {clID: [] for clID in cell_line_ids_no_spikein}
    clID__n_tumors_bootstraps = {clID: [] for clID in cell_line_ids_no_spikein}
    for b_m in range(B_mice):
        sampled_mice = np.random.choice(merged_sample_ids, size=n_mice, replace=True, p=weights)

        # now, want to get all tumors from the sampled mice & keep the clID info
        # initialize an empty list to collect the tumor sizes along with their clID information for the sampled mice
        sampled_sizes_with_clID = []
        for mouse in sampled_mice:
            mouse_tumors = tumor_sizes_grouped_by_mouse.get_group(mouse)
            sampled_sizes_with_clID.extend(zip(mouse_tumors['clID'], mouse_tumors['size'])) # adding (clID, size) pairs
        n_tumors = len(sampled_sizes_with_clID)

        # now, bootstrapping tumors
        for b_t in range(B_tumors):
            resampled_sizes = random.choices(sampled_sizes_with_clID, k=n_tumors)
            df = pd.DataFrame(resampled_sizes, columns=['clID', 'size']) # trying to speed up the whole thing by vectorizing

            # calculating relative tumor burden (sum of all sizes divided by total amount of cells injected, relative to the average across all cell lines)
            clID_sizes_summed = df.groupby('clID')['size'].sum()
            clID__burden = {clID: clID_sizes_summed.get(clID, 0) / (clID__n_0[clID] * n_mice) for clID in cell_line_ids_no_spikein} # initially all cell lines were present so no division error possible
            average_burden_across_cell_lines = np.mean(list(clID__burden.values()))

            for clID in cell_line_ids_no_spikein:
                clID__relative_burden_bootstraps[clID].append(clID__burden[clID] / average_burden_across_cell_lines) # this is the relative burden that goes to the output!

                # now, pooling tumors above minimim size
                clID__n_tumors_bootstraps[clID].append(df[(df['clID'] == clID) & (df['size'] >= min_tumor_size)].shape[0]) # only keeping the number of tumors for now

    return clID__relative_burden_bootstraps, clID__n_tumors_bootstraps


# getting the bootstrap distribution of y/x statistics
def bootstrap_index(clID__relative_burden_bootstraps_x, clID__relative_burden_bootstraps_y): # see above for definitions
    return {clID: np.array(clID__relative_burden_bootstraps_y[clID]) / np.array(clID__relative_burden_bootstraps_x[clID]) for clID in clID__relative_burden_bootstraps_x.keys()}

# functions to make the many plots:
def correlation_plot(cell_line_ids,
                     clID__x, clID__xerr,
                     clID__y, clID__yerr,
                     clID__color=clID__color_in_vivo_exp1, # do not forget to switch to the right experiment
                     clID__label=clID__label, # labels
                     xtitle='metric', ytitle='metric', # metrics to be displayed as axis labels
                     lims=[1e-3, 1e3], ticks=[1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3], labels=['1e-3', '1e-2', '1e-1', '1e0', '1e1', '1e2', '1e3']):
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
        xtitle, ytitle (string): Labels for the axes.
        lims (list): Limits for the axes (default lims=[1e-3, 1e3]).
        ticks (list): Locations for axes ticks (default ticks=[1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]).
        labels (list): Locations for axes tick labels (default labels=['1e-3', '1e-2', '1e-1', '1e0', '1e1', '1e2', '1e3']).

    Returns:
        ax (matplotlib.axes): A plot!
    """

    cell_line_ids_no_spikein = [clID for clID in clIDs if clID != 'GATC']

    fig, ax = plt.subplots(figsize=(4, 4))

    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['bottom', 'left']].set_linewidth(0.5)
    ax.tick_params(width=0.5)

    xs = [clID__x[x] for x in cell_line_ids_no_spikein]
    ys = [clID__y[y] for y in cell_line_ids_no_spikein]

    ax.scatter(xs, ys, color=[clID__color[z] for z in cell_line_ids_no_spikein], s=80, edgecolor='black', linewidth=0.5, alpha=1, zorder=1)
    ax.errorbar(xs, ys, xerr=[[clID__xerr[x][0] for x in cell_line_ids_no_spikein], [clID__xerr[x][1] for x in cell_line_ids_no_spikein]], yerr=[[clID__yerr[y][0] for y in cell_line_ids_no_spikein], [clID__yerr[y][1] for y in cell_line_ids_no_spikein]], color='black', linewidth=0.5, ls='none', zorder=-2)

    ax.plot(lims, lims, color='#c1c1c1', linewidth=0.5, zorder=-3)  # diagonal line
    ax.plot(lims, [1e0, 1e0], color='#c1c1c1', linewidth=0.5, alpha=0.3, zorder=-3)  # horizontal line
    ax.plot([1e0, 1e0], lims, color='#c1c1c1', linewidth=0.5, alpha=0.3, zorder=-3)  # vertical line


    spearmans_rho, spearmans_pval = sp.stats.spearmanr(xs, ys)  # rank correlation coefficient and p-value
    pearsons_r, pearsons_pval = sp.stats.pearsonr(xs, ys)  # Pearson correlation coefficient and p-value
    ax.text(2 * lims[0], 0.3 * lims[1], f'ρ = {spearmans_rho:.2f}, p = {spearmans_pval:.0e}', size=10)
    ax.text(2 * lims[0], 0.15 * lims[1], f'r = {pearsons_r:.2f}, p = {pearsons_pval:.0e}', size=10)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.minorticks_off()
    ax.set_xlim(lims[0], lims[1])
    ax.set_ylim(lims[0], lims[1])
    ax.set_xticks(ticks, labels, size=12)
    ax.set_yticks(ticks, labels, size=12)

    ax.set_xlabel(xtitle, fontsize=14, labelpad=6)
    ax.set_ylabel(ytitle, fontsize=14, labelpad=6)

    return ax

# will also define functions to separately calculate medians + confidence intervals based on the bootstrap output
# maybe also a function to calculate data means
