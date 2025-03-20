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

matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf');
font = matplotlib.font_manager.FontProperties(fname='/System/Library/Fonts/Helvetica.ttc')

# these two lines make text editable later if needed
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# cell line IDs, spike-in must be first
# excluded: other spike-ins & mt4-2D, which had the same clIDs
clIDs = ['GATC', 'AAGG', 'ACAC', 'ACCT', 'ACGA', 'ACTG', 'AGAG', 'AGCA', 'AGGT', 'AGTC', 'ATCG', 'ATGC', 'CAAC', 'CACT', 'CAGA', 'CATG', 'CCAA', 'CCTT', 'CGAT', 'CGTA', 'CTGT', 'CTTC', 'GAAG', 'GCAT', 'GCTA', 'GGAA', 'GGTT', 'GTAC', 'GTGA', 'GTTG', 'TCCA', 'TGAC', 'TTCC', 'TTGG']

# corresponding cell line labels
cell_line_labels = ['Spike-in 1', '7160c2', 'FC1199', '6694c2', '7160c5', 'BF857', 'BF1987', 'Panc02', '0688M', 'FC1245', 'BF4326', 'BF2117', '6419c5', 'BF4326Rep', 'KPC960', '0755P', 'BF1836', 'BF2014', 'BF5960', '6499c4', '6422c5', 'FC1245Rep', 'KPC-JH', 'KPC960Rep', 'BF2153', 'KC6141', 'mT3-2D', 'mT5-2D', '0764P', '6421c2', '0755A', '2838c3', 'FC1242', 'KPC961']

clID__label = dict(zip(clIDs, cell_line_labels))

# colors based on 3w F1 lung expansion from experiment one
clID__color_in_vivo = {'ACCT': (0.7568627450980392, 0.24313725490196078, 0.1803921568627451, 1.0),
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

# highlighting replicate or related cell lines in plots
clID_replicates__linestyle = {('AAGG', 'ACGA'): ':', ('AGTC', 'CTTC'): '-', ('ATCG', 'CACT'): '-', ('CAGA', 'GCAT'): '-'}

# now, some helpful functions:
def convert_barcode_reads_to_cell_counts(sample_ids, path_to_data,
                                         read_counts, metadata,
                                         clIDs=clIDs):
    """
    Converts barcode reads to estimated cell counts
    based on spiked in cell line read count.

    Args:
        sample_ids (list): List of samples; samples that should be merged
                           must be supplied together as in metadata 'notes' column.
        path_to_data (str): Path to the folder containing files with read counts.
        read_counts (pandas.DataFrame): A precomputed table with read counts (read_counts_by_cell_line.csv).
        metadata (pandas.DataFrame): A table containing experimental details, including initial cell numbers (metadata.csv).

    Returns:
        sample__counts (dict): Dictionary where samples are keys and values are dictionaries {clID: {barcode: cell count}}.
    """
    sample__counts = {}

    for sample in sample_ids:
        # when we need to merge reads first before converting to counts
        clID_bc__n_reads = {}; R = 0; one = False

        sub_samples = sample.split(':')

        for sub_sample in sub_samples:
            file_path = os.path.join(path_to_data, f"{sub_sample}_clIDs_rBC_cluster_counts.txt")
            if not os.path.exists(file_path):
                print(f"File {file_path} not found.")
                continue

            if metadata.loc[sub_sample, 'spike-in added'] == False:
                print(f"No spike-in added.")
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
                sample__counts[sample][clID][barcode] = 50000 / (1 + R) * n_reads

    return sample__counts


def bootstrap_relative_burden(samples, read_counts, cell_counts, metadata, clIDs=clIDs, B=10000):
    """
    Bootstrapps mice to estimate the relative metastatic tumor burden of different cell lines.

    Args:
        samples (list): Samples to include in the analysis;
                        merged samples need to be colon-separated.
        read_counts (pandas.DataFrame): A precomputed table with read counts (read_counts_by_cell_line.csv).
        cell_counts (pandas.DataFrame): A precomputed table with cell line counts (cell_counts_by_cell_line.csv).
        metadata (pandas.DataFrame): A table containing experimental details, including initial cell numbers (metadata.csv).
        clIDs (list): Four-nucleotide cell line IDs, optional.
        B (int): Number of bootstrap resampling iterations, optional (default=10000).

    Returns:
        clID__distr (dict): Bootstrapped distributions of normalized tumor burden (expansion)
                            across mice j for each cell line i (sum_{j} n_{ij} / n_{i}(0)).
        clID__distr_relative_to_mean (dict): Bootstrapped distributions of normalized tumor burden (expansion)
                                             across mice j for each cell line i (sum_{j} n_{ij} / n_{i}(0)),
                                             relative to the mean across all cell lines
                                             (normalized within each bootstrap realization).
        clID__exp_avg_relative_to_mean (dict): Average of the burden distribution of each cell line
                                               relative to the mean across all cell lines.
        clID__exp_err (dict): 0.95 confidence intervals for the relative burden estimates.
        clIDs_sorted (list): Cell line identifiers sorted by their relative burden.
    """

    clID__distr = {}  # bootstrap distributions of n_t / n_0
    clID__distr_relative_to_mean = {}
    clID__exp_avg_relative_to_mean = {} # weighted average across mice relative to average across all cell lines
    clID__exp_err = {} # 0.95 ci of the above

    n_samples = len(samples)
    experiment = metadata.loc[[s.split(':')[0] if ':' in s else s for s in samples], 'experiment'].iloc[0]
    pre_injection = metadata[(metadata['time point, d'] == 0) & (metadata['experiment'] == experiment)]

    # precompute pre-injection cell counts for each clID first,
    # by finding thir fractions and multiplying by the total number of cells injected
    n_0_dict = {clID: np.mean([read_counts[sample_id][clID] / sum(read_counts[sample_id][1:]) * pre_injection['initial number of cells'].mean()
                               for sample_id in pre_injection.index.tolist()]) for clID in clIDs if clID != 'GATC'}

    # bootstrapping mices to make one big mouse
    # this is equivalent to weighting each sample by the total number of cells and taking the avg
    bootstraps = []
    for clID, n_0 in n_0_dict.items():

        mouse_i_samples = np.random.choice(samples, size=(B, n_samples), replace=True)
        n_t = np.array([np.sum([cell_counts[mouse][clID] for mouse in mouse_i]) for mouse_i in mouse_i_samples])
        clID__distr[clID] = n_t / n_0
        bootstraps.append(clID__distr[clID])

    bootstraps = np.array(bootstraps)
    mean_across_cell_lines_per_bootstrap = np.mean(bootstraps, axis=0)

    original_mean_per_clID = {
        clID: np.mean([cell_counts[mouse][clID] / n_0_dict[clID] for mouse in samples])
        for clID in clIDs if clID != 'GATC'
    }
    original_mean = np.mean(list(original_mean_per_clID.values()))  # mean across cell lines for normalization


    for i, clID in enumerate(clID__distr):
        clID__distr_relative_to_mean[clID] = clID__distr[clID] / mean_across_cell_lines_per_bootstrap  # element-wise division

        if clID == 'GATC':
            continue

        ci_upper = np.percentile(clID__distr_relative_to_mean[clID], 97.5)
        ci_lower = np.percentile(clID__distr_relative_to_mean[clID], 2.5)

        clID__exp_err[clID] = [original_mean_per_clID[clID] / original_mean - ci_lower,
                               ci_upper - original_mean_per_clID[clID] / original_mean]
        clID__exp_avg_relative_to_mean[clID] = original_mean_per_clID[clID] / original_mean


    clIDs_sorted = [y[0] for y in sorted([[x, clID__exp_avg_relative_to_mean[x]] for x in clIDs if x != 'GATC'], key = lambda x: x[-1], reverse=True)]

    return clID__distr, clID__distr_relative_to_mean, clID__exp_avg_relative_to_mean, clID__exp_err, clIDs_sorted

# CHECK IF THIS IS RIGHT
def bootstrap_n_tumors(samples, path_to_data, read_counts, cell_counts, metadata, clIDs=clIDs, B=10000):
    """
    Bootstraps mice to estimate the relative seeding of different cell lines.

    Args:
        samples (list): Samples to include in the analysis;
                        merged samples need to be colon-separated.
        path_to_data (str): Path to the folder containing files with read counts.
        read_counts (pandas.DataFrame): A precomputed table with read counts (read_counts_by_cell_line.csv).
        cell_counts (pandas.DataFrame): A precomputed table with cell line counts (cell_counts_by_cell_line.csv).
        metadata (pandas.DataFrame): A table containing experimental details, including initial cell numbers (metadata.csv).
        clIDs (list): Four-nucleotide cell line IDs, optional.
        B (int): Number of bootstrap resampling iterations, optional (default=10000).

    Returns:
        clID__distr (dict): Bootstrapped distributions of normalized tumor number, relative to initial cell line proportion
                            across mice j for each cell line i (sum_{j} n_{ij} / n_{i}(0)).
        clID__distr_relative_to_mean (dict): Same as above but normalized within each bootstrap realization.
        clID__n_tumors_relative_to_mean (dict): Average of the tumor number distribution (bootstrapped) of each cell line
                                                relative to the mean across all cell lines.
        clID__err (dict): 0.95 confidence intervals for the relative burden estimates.
        clIDs_sorted (list): Cell line identifiers sorted by their relative burden.
    """

    # getting tumor sizes:
    sample__clID__barcode__count = convert_barcode_reads_to_cell_counts(samples, path_to_data=path_to_data, read_counts=read_counts, metadata=metadata)

    clID__n_tumors = {}
    for sample in samples:
        for i, clID in enumerate(clIDs[1:]):
            if clID not in clID__n_tumors:
                clID__n_tumors[clID] = []
            clID__n_tumors[clID].append(len([x for x in sample__clID__barcode__count[sample][clID].values()]))

    clID__distr = {}  # bootstrap distributions of the number of tumors
    clID__distr_relative_to_mean = {}
    clID__n_tumors_relative_to_mean = {} # weighted average across mice relative to average across all cell lines
    clID__err = {} # 0.95 ci of the above

    n_samples = len(samples)
    experiment = metadata.loc[[s.split(':')[0] if ':' in s else s for s in samples], 'experiment'].iloc[0]
    pre_injection = metadata[(metadata['time point, d'] == 0) & (metadata['experiment'] == experiment)]

    # precompute pre-injection cell counts for each clID first,
    # by finding their fractions and multiplying by the total number of cells injected
    n_0_dict = {clID: np.mean([read_counts[sample_id][clID] / sum(read_counts[sample_id][1:]) * pre_injection['initial number of cells'].mean()
                               for sample_id in pre_injection.index.tolist()]) for clID in clIDs if clID != 'GATC'}

    # bootstrapping mice to make one big mouse
    bootstraps = []
    n_mice = len(next(iter(clID__n_tumors.values())))
    for clID, n_0 in n_0_dict.items():

        mouse_i_samples = np.random.choice(range(n_mice), size=(B, n_mice), replace=True)
        n_t = np.array([sum(clID__n_tumors[clID][mouse] for mouse in mouse_i) for mouse_i in mouse_i_samples])
        clID__distr[clID] = n_t / n_0
        bootstraps.append(clID__distr[clID])

    bootstraps = np.array(bootstraps)

    # normalize each bootstrap sample (across cell lines) by the average of n_t / n_0
    normalized_bootstraps = bootstraps / np.mean(bootstraps, axis=1, keepdims=True)

    original_mean_per_clID = {
        clID: np.mean([len([x for x in sample__clID__barcode__count[sample][clID].values()]) for sample in samples])
        for clID in clIDs if clID != 'GATC'
    }
    original_mean = np.mean(list(original_mean_per_clID.values()))

    for i, clID in enumerate(clID__distr):
        clID__distr_relative_to_mean[clID] = normalized_bootstraps[i]

        if clID == 'GATC':
            continue

        ci_upper = np.percentile(clID__distr_relative_to_mean[clID], 97.5)
        ci_lower = np.percentile(clID__distr_relative_to_mean[clID], 2.5)

        clID__err[clID] = [np.mean(clID__distr_relative_to_mean[clID]) - ci_lower,
                           ci_upper - np.mean(clID__distr_relative_to_mean[clID])]

        clID__n_tumors_relative_to_mean[clID] = original_mean_per_clID[clID] / original_mean


    clIDs_sorted = [y[0] for y in sorted([[x, clID__n_tumors_relative_to_mean[x]] for x in clIDs if x != 'GATC'], key = lambda x: x[-1], reverse=True)]

    return clID__distr, clID__distr_relative_to_mean, clID__n_tumors_relative_to_mean, clID__err, clIDs_sorted


### CHECK THIS THIS MIGHT NOT BE RIGHT
def bootstrap_tumor_size_percentiles(samples, path_to_data, percentiles, read_counts, metadata, clIDs=clIDs, B=10000):
    """
    Bootstraps tumor sizes from multiple mice and creates an artificial big mouse with all tumor sizes pooled.
    Additionally, calculates percentiles of the distribution for each clID
    across pooled real mice and bootstraps to calculate 0.95 confidence intervals
    for those percentiles. This function also computes the cumulative sum of tumor sizes
    up to a given percentile and bootstraps the confidence intervals for these cumsums.

    Args:
        samples (list): Samples to include in the analysis;
                        merged samples need to be colon-separated.
        path_to_data (str): Path to the folder containing files with read counts.
        percentiles (list): List with percentiles to compute.
        read_counts (pandas.DataFrame): A precomputed table with read counts (read_counts_by_cell_line.csv).
        metadata (pandas.DataFrame): A table containing experimental details, including initial cell numbers (metadata.csv).
        clIDs (list): Four-nucleotide cell line IDs, optional.
        B (int): Number of bootstrap resampling iterations, optional (default=10000).

    Returns:
        clID__percentiles (dict): Percentiles of the tumor size distributions.
        clID__ci (dict): 0.95 confidence intervals for the percentiles.
        clID__cumsum_percentiles (dict): Cumulative sums up to the tumor size at each percentile.
        clID__cumsum_ci (dict): 0.95 confidence intervals for the cumulative sums up to each percentile.
    """
    clID__distr = {}
    clID__percentiles = {}
    clID__ci = {}
    clID__cumsum_percentiles = {}
    clID__cumsum_ci = {}

    # getting tumor sizes:
    sample__clID__barcode__count = convert_barcode_reads_to_cell_counts(samples, path_to_data=path_to_data, read_counts=read_counts, metadata=metadata)
    clID__tumor_sizes = {}

    for sample in samples:
        for i, clID in enumerate(clIDs[1:]):
            if clID not in clID__tumor_sizes:
                clID__tumor_sizes[clID] = {}
            clID__tumor_sizes[clID][sample] = [x for x in sample__clID__barcode__count[sample][clID].values()]

    for clID, sample__sizes in clID__tumor_sizes.items():
        clID__distr[clID] = [size for sample_sizes in sample__sizes.values() for size in sample_sizes]

        # compute percentiles
        clID__percentiles[clID] = np.percentile(clID__distr[clID], percentiles)

        # calculate the cumulative tumor burden up to each percentile
        sorted_sizes = np.sort(clID__distr[clID])

        cumulative_sum_up_to_percentile = []
        for p in clID__percentiles[clID]:
            sizes_up_to_percentile = sorted_sizes[sorted_sizes <= p]
            cumulative_sum_up_to_percentile.append(np.sum(sizes_up_to_percentile) / sum(sorted_sizes))

        clID__cumsum_percentiles[clID] = cumulative_sum_up_to_percentile

        # bootstrapping for percentiles:
        bootstrap_percentiles = []
        bootstrap_cumsums = []

        for _ in range(B):
            resample = np.random.choice(clID__distr[clID], size=len(clID__distr[clID]), replace=True)
            resample_sorted = np.sort(resample)
            bootstrap_percentiles.append(np.percentile(resample, percentiles))

            bootstrap_cumsum_up_to_percentile = []
            for p in bootstrap_percentiles[-1]:
                # cumulative sum of sizes up to this percentile size in the resampled data
                resample_sizes_up_to_percentile = resample_sorted[resample_sorted <= p]

                bootstrap_cumsum_up_to_percentile.append(sum(resample_sizes_up_to_percentile) / sum(resample_sorted))

            bootstrap_cumsums.append(bootstrap_cumsum_up_to_percentile)

        ci_lower = np.percentile(bootstrap_percentiles, 2.5, axis=0)
        ci_upper = np.percentile(bootstrap_percentiles, 97.5, axis=0)
        clID__ci[clID] = (ci_lower, ci_upper)

        cumsum_lower = np.percentile(bootstrap_cumsums, 2.5, axis=0)
        cumsum_upper = np.percentile(bootstrap_cumsums, 97.5, axis=0)
        clID__cumsum_ci[clID] = (cumsum_lower, cumsum_upper)

    return clID__percentiles, clID__ci, clID__cumsum_percentiles, clID__cumsum_ci



# functions to make the many plots:
def metric_bar_plot(clIDs_ordered, clID__y, clID__err, ytitle='metric', # metric (to be displayed on y-axis)
                    clID__color=clID__color_in_vivo, clID__label=clID__label, # colors / labels
                    clID_replicates__linestyle=clID_replicates__linestyle, # how to highlight replicates
                    ylims=[1e-2, 1e1], yticks=[1e-2, 1e-1, 1e0, 1e1], ylabels=['0.01', '0.1', '1', '10']):
    """
    Makes a bar plot with some kind of metric on a log-scale y-axis, and cell lines on the x-axis.

    Args:
        clIDs_ordered (list): Ordered list of cell line IDs, must exclude spike-in.
        clID__y (dict): Dictionary with cell line IDs as keys and mean metric values as values (e.g. 'AAGT': 0.02).
        clID__err (dict): Dictionary with cell line IDs as keys
                          and -/+ errors from the mean as values (e.g. 'AAGT': [0.01, 0.01]).
        ytitle (string): Label for the y-axis.
        clID__color (dict): Dictionary mapping cell line IDs to colors (by default, the one defined above).
        clID__label (dict): Dictionary mapping cell line IDs to label names.
        clID_replicates__linestyle (dict): Dictionary mapping cell line replicates to connecting line style.
        ylims (list): Limits for the y-axis (default ylims=[1e-2, 1e1]).
        yticks (list): Locations for y-axis ticks (default yticks=[1e-2, 1e-1, 1e0, 1e1]).
        ylabels (list): Locations for y-axis tick labels (default ylabels=['0.01', '0.1', '1', '10']).

    Returns:
        ax (matplotlib.axes): A plot!
    """

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines['bottom'].set_position(('axes', -0.025))
    ax.spines['left'].set_position(('axes', -0.025))

    ax.spines[['bottom', 'left']].set_linewidth(0.5)
    ax.yaxis.set_tick_params(width=0.5)

    heights, bottoms, ys, errs, cs = [], [], [], [[], []], []

    for i in clIDs_ordered:
        cs.append(clID__color[i])
        ys.append(clID__y[i])
        errs[0].append(clID__err[i][0])
        errs[1].append(clID__err[i][1])

        if clID__y[i] > 1:
            bottoms.append(1)
            heights.append(clID__y[i] - 1)
        else:
            bottoms.append(clID__y[i])
            heights.append(1 - clID__y[i])

    ax.bar(range(len(clIDs_ordered)), heights, width=0.8, bottom=bottoms, color=cs,
           edgecolor='black', linewidth=0.5, clip_on=False)

    ax.errorbar(range(len(clIDs_ordered)), ys, yerr=errs, color='black', linewidth=0.5,
                alpha=1, zorder=1, ls='none', clip_on=False)

    ax.set_yscale('log')
    ax.minorticks_off()
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, size=12)
    ax.set_ylim(ylims[0], ylims[1])

    ax.set_xlim(0, len(clIDs_ordered) - 1)
    ax.set_xticks(range(len(clIDs_ordered)))
    ax.set_xticklabels([clID__label[x] for x in clIDs_ordered], size=12, rotation=90)

    ax.set_ylabel(ytitle, fontsize=14, labelpad=6)

    # connecting replicate cell lines by lines
    for key, value in clID_replicates__linestyle.items():
        y = ylims[0] * 0.08
        x = sorted([clIDs_ordered.index(key[0]), clIDs_ordered.index(key[1])])
        ax.annotate('', xy=(x[0], y), xycoords='data', xytext=(x[1], y), textcoords='data',
                    arrowprops=dict(arrowstyle='-', linestyle=value, linewidth=0.5, color='0',
                    patchA=None, patchB=None, connectionstyle='bar,angle=180,fraction=-0.5'),
                    annotation_clip=False);
    return ax

def comparison_scatter_plot(clIDs, # always exclude spike-in, sometimes Panc02 as well
                            clID__x, clID__xerr,
                            clID__y, clID__yerr,
                            xtitle='metric', ytitle='metric', # metric (to be displayed on y-axis)
                            clID__color=clID__color_in_vivo, clID__label=clID__label, # colors / labels
                            lims=[1e-3, 1e3], ticks=[1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3], labels=['0.001', '0.01', '0.1', '1', '10', '100', '1000']):
    """
    Makes a scatter plot comparing two experiments or conditions,
    and calculates Spearman and Pearson correlation coefficients with p-values.

    Makes a bar plot with some kind of metric with on a log-scale y-axis, and cell lines on the x-axis.

    Parameters:
        clIDs (list): Cell line IDs, in any order; all of these cell lines will be used to compute correlations.
        clID__x, clID__y (dict): Dictionaries with cell line IDs as keys and mean metric values as values (e.g. 'AAGT': 0.02).
        clID__xerr, clID__yerr (dict): Dictionaries with cell line IDs as keys
                                       and -/+ errors from the mean as values (e.g. 'AAGT': [0.01, 0.01]).
        clID__color (dict): Dictionary mapping cell line IDs to colors (by default, the one defined above).
        clID__label (dict): Dictionary mapping cell line IDs to label names.
        xtitle, ytitle (string): Labels for the axes.
        lims (list): Limits for the axes (default lims=[1e-3, 1e3]).
        ticks (list): Locations for axes ticks (default ticks=[1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]).
        ylabels (list): Locations for axes tick labels (default ylabels=labels=['0.001', '0.01', '0.1', '1', '10', '100', '1000']).

    Returns:
        ax (matplotlib.axes): A plot!
    """

    fig, ax = plt.subplots(figsize=(4, 4))

    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['bottom', 'left']].set_linewidth(0.5)
    ax.tick_params(width=0.5)

    xs = [clID__x[x] for x in clIDs]
    ys = [clID__y[y] for y in clIDs]

    ax.scatter(xs, ys, color=[clID__color[z] for z in clIDs], s=80, edgecolor='black', linewidth=0.75, alpha=1, zorder=1)
    ax.errorbar(xs, ys, xerr=[[clID__xerr[x][0] for x in clIDs], [clID__xerr[x][1] for x in clIDs]], yerr=[[clID__yerr[y][0] for y in clIDs], [clID__yerr[y][1] for y in clIDs]], color='black', linewidth=0.75, ls='none', zorder=-2)

    ax.plot(lims, lims, color='#c1c1c1', zorder=-3)  # diagonal Line
    ax.plot(lims, [1e0, 1e0], color='#c1c1c1', alpha=0.3, zorder=-3)  # horizontal Line
    ax.plot([1e0, 1e0], lims, color='#c1c1c1', alpha=0.3, zorder=-3)  # vertical Line


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
