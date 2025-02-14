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
# excluded: other spike-ins & mt4-2D, which hade the same clIDs
clIDs = ['GATC', 'AAGG', 'ACAC', 'ACCT', 'ACGA', 'ACTG', 'AGAG', 'AGCA', 'AGGT', 'AGTC', 'ATCG', 'ATGC', 'CAAC', 'CACT', 'CAGA', 'CATG', 'CCAA', 'CCTT', 'CGAT', 'CGTA', 'CTGT', 'CTTC', 'GAAG', 'GCAT', 'GCTA', 'GGAA', 'GGTT', 'GTAC', 'GTGA', 'GTTG', 'TCCA', 'TGAC', 'TTCC', 'TTGG']

# corresponding cell line labels
cell_line_labels = ['Spike-in 1', '7160c2', 'FC1199', '6694c2', '7160c5', 'BF857', 'BF1987', 'Panc02', '0688M', 'FC1245', 'BF4326', 'BF2117', '6419c5', 'BF4326', 'KPC960', '0755P', 'BF1836', 'BF2014', 'BF5960', '6499c4', '6422c5', 'FC1245', 'KPC-JH', 'KPC960', 'BF2153', 'KC6141', 'mT3-2D', 'mT5-2D', '0764P', '6421c2', '0755A', '2838c3', 'FC1242', 'KPC961']

clID__label = dict(zip(clIDs, cell_line_labels))

# colors based on in vitro expansion
clID__color = {'AGGT': (0.7568627450980392, 0.24313725490196078, 0.1803921568627451, 1.0),
               'AGCA': (0.7712572087658592, 0.2597462514417531, 0.19663206459054208, 1.0),
               'CTTC': (0.9545098039215686, 0.5642445213379469, 0.48705882352941177, 1.0),
               'AGTC': (0.9592156862745098, 0.5794694348327566, 0.5011764705882353, 1.0),
               'ATCG': (0.973517877739331, 0.6467358708189158, 0.5638292964244522, 1.0),
               'GAAG': (0.9646597462514417, 0.6810611303344867, 0.5964936562860439, 1.0),
               'CACT': (0.9498961937716263, 0.7382698961937716, 0.6509342560553634, 1.0),
               'ACCT': (0.9498961937716263, 0.7382698961937716, 0.6509342560553634, 1.0),
               'CCTT': (0.9439907727797, 0.7611534025374858, 0.6727104959630913, 1.0),
               'CATG': (0.9351326412918108, 0.7954786620530565, 0.7053748558246828, 1.0),
               'GGAA': (0.9174163783160323, 0.8641291810841983, 0.7707035755478662, 1.0),
               'ACAC': (0.9144636678200692, 0.8755709342560554, 0.7815916955017301, 1.0),
               'GTTG': (0.8975778546712803, 0.8717416378316032, 0.7831603229527104, 1.0),
               'AGAG': (0.8975778546712803, 0.8717416378316032, 0.7831603229527104, 1.0),
               'TTGG': (0.8265282583621683, 0.8423068050749711, 0.7780853517877739, 1.0),
               'TCCA': (0.7748558246828141, 0.8208996539792387, 0.77439446366782, 1.0),
               'CTGT': (0.768396770472895, 0.8182237600922722, 0.7739331026528258, 1.0),
               'TTCC': (0.7361014994232987, 0.8048442906574395, 0.7716262975778546, 1.0),
               'CAGA': (0.7296424452133794, 0.8021683967704729, 0.7711649365628604, 1.0),
               'GCAT': (0.7167243367935409, 0.7968166089965397, 0.7702422145328719, 1.0),
               'CGTA': (0.6521337946943483, 0.7700576701268742, 0.7656286043829296, 1.0),
               'GGTT': (0.6328489042675893, 0.7613840830449826, 0.7628604382929642, 1.0),
               'ACGA': (0.6264821222606689, 0.7580622837370241, 0.7610149942329872, 1.0),
               'AAGG': (0.5564475201845444, 0.7215224913494809, 0.740715109573241, 1.0),
               'GTAC': (0.4800461361014994, 0.6816608996539792, 0.7185697808535179, 1.0),
               'CGAT': (0.4609457900807381, 0.6716955017301037, 0.7130334486735871, 1.0),
               'GCTA': (0.44184544405997694, 0.6617301038062283, 0.7074971164936563, 1.0),
               'CAAC': (0.44184544405997694, 0.6617301038062283, 0.7074971164936563, 1.0),
               'TGAC': (0.40364475201845446, 0.6417993079584775, 0.6964244521337947, 1.0),
               'ATGC': (0.39727797001153403, 0.638477508650519, 0.6945790080738178, 1.0),
               'ACTG': (0.3845444059976932, 0.6318339100346021, 0.6908881199538639, 1.0),
               'GTGA': (0.1734717416378316, 0.5384083044982699, 0.615686274509804, 1.0),
               'CCAA': (0.10408304498269896, 0.5081430219146482, 0.5905882352941176, 1.0)}

# now, some helpful functions:
def convert_barcode_reads_to_cell_counts(sample_ids, path_to_data,
                                         read_counts, metadata,
                                         clIDs=clIDs):
    """
    Converts barcode reads to estimated cell counts based on spiked in cell line read count.
    
    Args:
        sample_ids (list): List of samples; 
                  samples that should be merged must be supplied together as in metadata 'notes' column.
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
            file_path = os.path.join(path_to_data, f"{sub_sample}_merged_clIDs_bc_clusters_counts.txt")
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
    Bootstrapps mice to estimate the relative metastatic or tumor burden of different cell lines.

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
        clID__exp_avg_relative_to_mean (dict): Average of the burden distribution of each cell line
                                               relative to the mean across all cell lines.
        clID__exp_err (dict): 0.95 confidence intervals for the relative burden estimates.
        clIDs_sorted (list): Cell line identifiers sorted by their relative burden.
    """
    
    clID__distr = {}  # bootstrap distributions of n_t / n_0
    clID__exp_avg_relative_to_mean = {} # weighted average across mice relative to average across all cell lines
    clID__exp_err = {} # 0.95 ci of the above
    
    n_samples = len(samples)
    experiment = metadata.loc[[s for s in samples if ':' not in s], 'experiment'].iloc[0]    
    pre_injection = metadata[(metadata['time point, d'] == 0) & (metadata['experiment'] == experiment)]

    # precompute pre-injection cell counts for each clID first,
    # by finding thir fractions and multiplying by the total number of cells injected
    n_0_dict = {clID: np.mean([read_counts[sample_id][clID] / sum(read_counts[sample_id][1:]) * pre_injection['initial number of cells'].mean() 
                               for sample_id in pre_injection.index.tolist()]) for clID in clIDs if clID != 'GATC'}

    # bootstrapping mices to make one big mouse
    # this is equivalent to weighting each sample by the total number of cells and taking the avg
    avgs = []
    for clID, n_0 in n_0_dict.items():
        
        mouse_i_samples = np.random.choice(samples, size=(B, n_samples), replace=True)
        n_t = np.array([sum(cell_counts[mouse][clID] for mouse in mouse_i) for mouse_i in mouse_i_samples])
        clID__distr[clID] = n_t / (n_0 * B)
        avgs.append(clID__distr[clID])
        
    mean_across_cell_lines = np.mean(avgs)

    for clID in clIDs:
        if clID == 'GATC':
            continue
        y = np.mean(clID__distr[clID])
    
        ci_upper = np.percentile(clID__distr[clID], 97.5)
        ci_lower = np.percentile(clID__distr[clID], 2.5)
        clID__exp_err[clID] = [(y - ci_lower) / mean_across_cell_lines, (ci_upper - y) / mean_across_cell_lines]
        clID__exp_avg_relative_to_mean[clID] = y / mean_across_cell_lines
    
    clIDs_sorted = [y[0] for y in sorted([[x, clID__exp_avg_relative_to_mean[x]] for x in clIDs if x != 'GATC'], key = lambda x: x[-1], reverse=True)]

    return clID__distr, clID__exp_avg_relative_to_mean, clID__exp_err, clIDs_sorted


# functions to make the many plots:
def metric_bar_plot(clIDs_ordered, clID__y, clID__err, ytitle='metric', # metric (to be displayed on y-axis)
                    clID__color=clID__color, clID__label=clID__label, # colors / labels
                    ylims=[1e-2, 1e1], yticks=[1e-2, 1e-1, 1e0, 1e1], ylabels=['0.01', '0.1', '1', '10']):
    """
    Makes a bar plot with some kind of metric with on a log-scale y-axis, and cell lines on the x-axis.
    
    Args:
        clIDs_ordered (list): Ordered list of cl IDs, must exclude spike-in.
        clID__y (dict): Dictionary with cell line IDs as keys and mean metric values as values (e.g. 'AAGT': 0.02).
        clID__err (dict): Dictionary with cell line IDs as keys 
                          and -/+ errors from the mean as values (e.g. 'AAGT': [0.01, 0.01]).
        ytitle (string): Label for the y-axis.
        clID__color (dict): Dictionary mapping cell line IDs to colors (by default, the one defined above).
        clID__label (dict): Dictionary mapping cell line IDs to label names.
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
           edgecolor='black', linewidth=0.75, clip_on=False)
    
    ax.errorbar(range(len(clIDs_ordered)), ys, yerr=errs, color='black', linewidth=0.75, 
                alpha=1, zorder=1, ls='none', clip_on=False)

    ax.set_yscale('log')
    ax.minorticks_off()
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, size=12)
    ax.set_ylim(ylims[0], ylims[1])

    ax.set_xlim(0, len(clIDs_ordered) - 1)
    ax.set_xticks(range(len(clIDs_ordered)))
    ax.set_xticklabels([clID__label[x].replace('_', ' ') for x in clIDs_ordered], size=12, rotation=90)

    ax.set_ylabel(ytitle, fontsize=14, labelpad=6)

    return ax

def comparison_scatter_plot(clIDs, # always exclude spike-in, sometimes Panc02 as well
                            clID__x, clID__xerr, 
                            clID__y, clID__yerr, 
                            xtitle='metric', ytitle='metric', # metric (to be displayed on y-axis)
                            clID__color=clID__color, clID__label=clID__label, # colors / labels
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
    
    xs = [clID__x[x] for x in clIDs]
    ys = [clID__y[y] for y in clIDs]
    
    ax.scatter(xs, ys, color=[clID__color[z] for z in clIDs], s=80, edgecolor='black', linewidth=0.75, alpha=1, zorder=1)
    ax.errorbar(xs, ys, xerr=[[clID__xerr[x][0] for x in clIDs], [clID__xerr[x][1] for x in clIDs]], yerr=[[clID__yerr[y][0] for y in clIDs], [clID__yerr[y][1] for y in clIDs]], color='black', linewidth=0.75, ls='none', zorder=-2)

    ax.plot(lims, lims, color='#c1c1c1', zorder=-3)  # diagonal Line
    ax.plot(lims, [1e0, 1e0], color='#c1c1c1', alpha=0.3, zorder=-3)  # horizontal Line
    ax.plot([1e0, 1e0], lims, color='#c1c1c1', alpha=0.3, zorder=-3)  # vertical Line


    spearmans_rho, spearmans_pval = sp.stats.spearmanr(xs, ys)  # rank correlation coefficient and p-value
    pearsons_r, pearsons_pval = sp.stats.pearsonr(xs, ys)  # Pearson correlation coefficient and p-value
    ax.text(2 * lims[0], 0.3 * lims[1], f'ρ = {spearmans_rho:.2f}, p = {spearmans_pval:.2f}', size=10)
    ax.text(2 * lims[0], 0.15 * lims[1], f'r = {pearsons_r:.2f}, p = {pearsons_pval:.2f}', size=10)

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
