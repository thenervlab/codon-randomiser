import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap
import matplotlib.patches as mpatches
import math
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import pingouin as pg

from loguru import logger
logger.info('Import OK')

syn_codon_usage = 'Results_cleaned/codon-usage-syn.csv'
original_codon_usage = 'Results_cleaned/codon-usage-original-sequence.csv'

syn_protein_codon= pd.read_csv(syn_codon_usage)
original_protein_codon = pd.read_csv(original_codon_usage)

#------------------ro protein codon usage------------------------
# comparison of all optimal frequencies
#drop stop codons
syn_protein_codon = syn_protein_codon[~syn_protein_codon['aminoacid'].str.contains('stp', case=False, na=False)]

#get absolute difference of optimal_freq-ro
syn_protein_codon['syn_codon_usage'] = syn_protein_codon['syn_nucleotide_identity']/syn_protein_codon['optimal_frequency']

ro_codon_comparison = syn_protein_codon[['syn_codon', 'aminoacid', 'syn_codon_usage', 'protein']].copy()

ro_dict = {'ENV':'Final ENV', 'GAG':'Final GAG', 'POL':'Final POL', 'PRO':'Final PRO'}

ro_codon_comparison['sequence'] = ro_codon_comparison['protein'].map(ro_dict) 


#------------------original protein codon usage------------------------
# comparison of all optimal frequencies
#drop stop codons
original_protein_codon = original_protein_codon[~original_protein_codon['aminoacid'].str.contains('stp', case=False, na=False)]

#get absolute difference of optimal_freq-ro
original_protein_codon['original_codon_usage'] = original_protein_codon['nucleotide_identity']/original_protein_codon['optimal_frequency']
original_codon_comparison = original_protein_codon[['original_codon', 'aminoacid', 'original_codon_usage', 'protein']].copy()
original_dict = {'ENV':'ENVcon', 'GAG':'GAGcon', 'POL':'POLcon', 'PRO':'PROcon'}

original_codon_comparison['sequence'] = original_codon_comparison['protein'].map(original_dict) 

#read in rEGFP codon usage
rENV_comparison = pd.read_csv('Results_cleaned/Misc-sequence-codon-usages/codon-usage-randomised_ENV-sequence.csv')
oENV_TUB_comparison = pd.read_csv('Results_cleaned/Misc-sequence-codon-usages/codon-usage-TUBA1A-oENV.csv')
EGFP_comparison = pd.read_csv('Results_cleaned/Misc-sequence-codon-usages/codon-usage-EGFP-sequences.csv')

concat_comparison = pd.concat([ro_codon_comparison, original_codon_comparison, rENV_comparison, oENV_TUB_comparison, EGFP_comparison], ignore_index=True)

concat_comparison['all_comparisons'] = concat_comparison['syn_codon_usage'].fillna(0) + concat_comparison['original_codon_usage'].fillna(0)+ concat_comparison['rand_RSCU'].fillna(0) + concat_comparison['RSCU'].fillna(0) + concat_comparison['new_RSCU'].fillna(0)

concat_comparison.to_csv('codon-optimacy-comparison.csv', index=False)

# -------------Highest optimal codon usage value for amino acid in sequence-------------------
idx = concat_comparison.groupby(['sequence', 'aminoacid'])['all_comparisons'].idxmax()

# Convert idx to a DataFrame containing the relevant rows
aminoacid_max_values = concat_comparison.loc[idx].reset_index(drop=True)

#---------------------------------Statistical analysis---------------------------------------
pairs = [
    ('Final ENV', 'ENVcon'),
    ('Final GAG', 'GAGcon'),
    ('Final POL', 'POLcon'),
    ('Final PRO', 'PROcon'),
    ('rENV', 'Final ENV'),
    ('oENV', 'rENV'),
    ('TUBA1A', 'Final ENV'),
    ('rEGFP', 'EGFP')
]

grouped = aminoacid_max_values.groupby('sequence')['all_comparisons'].apply(list)
aov = pg.anova(data=aminoacid_max_values, dv='all_comparisons', between='sequence', detailed=True)

pairwise = pg.pairwise_tukey(
    data=aminoacid_max_values, 
    dv='all_comparisons', 
    between='sequence', 
    effsize='hedges'
)

pairwise = pairwise[(pairwise['p-tukey']<0.05)]
#------------------------------------------Visualisation----------------------------------------
# Palette
NRV_palette_custom = {
    'Final POL': '#520B57',  
    'Final GAG': '#80A4DA', 
    'Final PRO': '#0BBDA3', 
    'Final ENV': '#FF960A',
    'POLcon': '#000000',
    'ENVcon': '#000000',
    'PROcon': '#000000',
    'GAGcon': '#000000',
    'rENV': '#EA4A07',
    'oENV':'#A00935',
    'TUBA1A': '#915A95',
    'EGFP': "#17DE24",
    'rEGFP': "#2F1969"
}

# Define Barplot function
def bardotplot(data, xcol, ycol, order, hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, palette=False, xlabel='', ylabel=False, pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.2, cap_width=2):
    if ax is None:
        fig, ax = plt.subplots()
    if hue == None:
        dodge = False
    else:
        dodge = True
    sns.barplot(
        data=data,
        x=xcol,
        y=ycol,
        hue=hue,
        palette=palette,
        capsize=cap_size,
        errwidth=cap_width,
        ax=ax,
        dodge=dodge,
        order=order,
        hue_order=hue_order,
        edgecolor='white'
    )
    sns.stripplot(
        data=data,
        x=xcol,
        y=ycol,
        hue=scat_hue,
        palette=palette,
        ax=ax,
        edgecolor='#fff',
        linewidth=1,
        s=dot_size,
        order=order,
        hue_order=scat_hue_order,
        dodge=dodge,
    )

    ax.set(ylabel=ylabel)

    ax.set_xlabel(xlabel)
    if xticks:
        ax.set_xticks(xticks)
        ax.set_xticklabels(hue_order*len(order))
    if groups:
        for group_label, (x0, x1, x2) in groups.items():
            ax.annotate(group_label, xy=(x0, group_label_y),
                        xycoords='data', ha='center', annotation_clip=False)
            trans = ax.get_xaxis_transform()
            ax.plot([x1, x2], [group_line_y, group_line_y],
                    color="black", transform=trans, clip_on=False)

    if legend == '':
        ax.legend('', frameon=False)
    else:    
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
    
    return ax

plt.rcParams['svg.fonttype']= 'none'

# Create plot
fig, ax = plt.subplots(figsize=(7.7, 2.6))
bardotplot(
    data=aminoacid_max_values, 
    xcol='sequence', 
    ycol='all_comparisons', 
    order=['GAGcon','Final GAG','PROcon','Final PRO','POLcon', 'Final POL','ENVcon','Final ENV','rENV', 'oENV', 'EGFP', 'rEGFP', 'TUBA1A'], 
    hue=None, 
    hue_order=None, 
    scat_hue=None, 
    scat_hue_order=None, 
    palette=NRV_palette_custom, 
    xlabel='Sequence', 
    ylabel='Delta optimal codon usage', 
    xticks=None, 
    groups=None, 
    group_label_y=-0.18, 
    group_line_y=-0.05, 
    legend='', 
    dot_size=5, 
    cap_size=0.2, 
    cap_width=1,
    ax=ax)

plt.ylim(   )
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
ax.set_xlabel('Sequence', fontsize=10)
ax.set_ylabel('Maximum Aminoacid suitability', fontsize=10)

plt.savefig('Relative-adaptiveness-codon-usage.svg')

