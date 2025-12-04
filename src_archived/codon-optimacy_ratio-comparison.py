import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap
from loguru import logger
import matplotlib.patches as mpatches
from statannotations.Annotator import Annotator



logger.info('Import OK')

ro_codon_usage = 'results/codon-usage.csv'
original_codon_usage = 'results/codon-usage-original-sequence.csv'



ro_protein_codon= pd.read_csv(ro_codon_usage)
original_protein_codon = pd.read_csv(original_codon_usage)






#colour palette

NRV_palette = {
    '': '#000000', 
    '': '#520B57', 
    '': '#915A95', 
    '': '#80A4DA', 
    '': '#0BBDA3', 
    '': '#FF960A',
    '':'#EA4A07',
    '':'#A00935'
}

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
}


#------------------ro protein codon usage------------------------
# comparison of all optimal frequencies
#drop stop codons
ro_protein_codon = ro_protein_codon[~ro_protein_codon['aminoacid'].str.contains('stp', case=False, na=False)]

#get absolute difference of optimal_freq-ro
ro_protein_codon['ro_RSCU'] = abs(ro_protein_codon['optimal_freq-ro'])

ro_codon_comparison = ro_protein_codon[['ro_codon', 'aminoacid', 'ro_RSCU', 'protein']].copy()


ro_dict = {'ENV':'Final ENV', 'GAG':'Final GAG', 'POL':'Final POL', 'PRO':'Final PRO'}

ro_codon_comparison['sequence'] = ro_codon_comparison['protein'].map(ro_dict) 


#------------------original protein codon usage------------------------
# comparison of all optimal frequencies
#drop stop codons

import math
original_protein_codon = original_protein_codon[~original_protein_codon['aminoacid'].str.contains('stp', case=False, na=False)]


#get absolute difference of optimal_freq-ro
original_protein_codon['original_RSCU'] = abs(original_protein_codon['optimal_freq-original'])

original_codon_comparison = original_protein_codon[['original_codon', 'aminoacid', 'original_RSCU', 'protein']].copy()

original_dict = {'ENV':'ENVcon', 'GAG':'GAGcon', 'POL':'POLcon', 'PRO':'PROcon'}


original_codon_comparison['sequence'] = original_codon_comparison['protein'].map(original_dict) 

#read in rEGFP codon usage
rENV_comparison = pd.read_csv('codon-usage-rENV-sequence.csv')

concat_comparison = pd.concat([ro_codon_comparison, original_codon_comparison, rENV_comparison], ignore_index=True)


concat_comparison['all_comparisons'] = concat_comparison['ro_RSCU'].fillna(0) + concat_comparison['original_RSCU'].fillna(0)+ concat_comparison['optimal_freq-new'].fillna(0)



# Get the index of the highest all_comparisons value for each sequence and aminoacid
idx = concat_comparison.groupby(['sequence', 'aminoacid'])['all_comparisons'].idxmax()

# Create new DataFrame with those rows
highest_comparison_df = concat_comparison.loc[idx].reset_index(drop=True)




#----------------------REMOVE CODONS THAT MIGHT BE USED IN TGTG --------------


high_diff_comparison = concat_comparison[concat_comparison['all_comparisons'] > 0.3].copy()


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

    if pairs:
        annotator = Annotator(
            ax=ax, pairs=pairs, data=data, x=xcol, y=ycol, order=order, hue=hue, hue_order=hue_order)
        annotator.configure(test='t-test_ind', text_format='star',
                            loc='inside', comparisons_correction=correction, line_width=0.5)
        annotator.apply_and_annotate()

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


import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype']= 'none'

# Create plot
fig, ax = plt.subplots(figsize=(7.7, 2.6))


pairs = [
    ('Final ENV', 'ENVcon'),
    ('Final GAG', 'GAGcon'),
    ('Final POL', 'POLcon'),
    ('Final PRO', 'PROcon'),
    ('rENV', 'Final ENV'),
]

bardotplot(
    data=highest_comparison_df, 
    xcol='sequence', 
    ycol='all_comparisons', 
    order=['GAGcon','Final GAG','PROcon','Final PRO','POLcon', 'Final POL','ENVcon','Final ENV','rENV'], 
    hue=None, 
    hue_order=None, 
    scat_hue=None, 
    scat_hue_order=None, 
    palette=NRV_palette_custom, 
    xlabel='Sequence', 
    ylabel='Delta optimal codon usage', 
    pairs=pairs, 
    correction='holm-bonferroni', 
    xticks=None, 
    groups=None, 
    group_label_y=-0.18, 
    group_line_y=-0.05, 
    legend='', 
    dot_size=5, 
    cap_size=0.2, 
    cap_width=1,
    ax=ax)

plt.ylim()
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
ax.set_xlabel('Sequence', fontsize=8)
ax.set_ylabel('Codon relative adaptiveness', fontsize=8)

plt.savefig('Relative-adaptiveness-codon-usage.svg')

