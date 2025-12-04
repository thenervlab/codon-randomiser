import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap
from loguru import logger
import matplotlib.patches as mpatches

logger.info('Import OK')

input_folder = 'results/codon-usage.csv'

protein_codon_all = pd.read_csv(input_folder)

# Palette
NRV_palette_custom = {
    'POL': '#520B57',  
    'GAG': '#80A4DA', 
    'PRO': '#0BBDA3', 
    'ENV': '#FF960A',
}

# comparison of all optimal frequencies
# Remove stop codon
protein_codon_all = protein_codon_all[~protein_codon_all['aminoacid'].str.contains('stp', case=False, na=False)]

#calculate log2(ro/o)
protein_codon_all['comparison'] = protein_codon_all['ro_nucleotide_identity'] / protein_codon_all['optimal_frequency']
protein_codon_all['log2_comparison'] = np.log2(protein_codon_all['comparison'])
aminoacid_codon_usages = [group.copy() for _, group in protein_codon_all.groupby('aminoacid')]

#Data sets for subplots
subplot_dictionary = {
    'Arg': 'subplot_1', 'Ala': 'subplot_1', 'Asn': 'subplot_1', 'Asp': 'subplot_1', 'Cys': 'subplot_1',
    'Gln': 'subplot_2', 'Glu': 'subplot_2', 'Gly': 'subplot_2', 'His': 'subplot_2', 'Leu': 'subplot_2',
    'Ile': 'subplot_3', 'Lys': 'subplot_3', 'Met': 'subplot_3', 'Pro': 'subplot_3', 'Ser': 'subplot_3',
    'Phe': 'subplot_4', 'Thr': 'subplot_4', 'Trp': 'subplot_4', 'Tyr': 'subplot_4', 'Val': 'subplot_4'
}


for df in aminoacid_codon_usages:
    aa = df['aminoacid'].iloc[0]
    df['subplot'] = subplot_dictionary.get(aa, None)


#groupby subplots
all_codon_usages = pd.concat(aminoacid_codon_usages, ignore_index=True)
subplot_groups = dict(tuple(all_codon_usages.groupby('subplot')))

all_codon_usages.sort_values(['subplot', 'aminoacid', 'ro_codon'], inplace=True)

# codon_order = all_codon_usages['ro_codon'].unique()
# all_codon_usages['ro_codon'] = pd.Categorical(
#     all_codon_usages['ro_codon'],
#     categories=codon_order,
#     ordered=True
# )

#PLOTTING


fig, axes = plt.subplots(4,1, figsize=(7, 2.4*4), gridspec_kw={'hspace': 0.4})
subplot_order = ['subplot_1', 'subplot_2', 'subplot_3', 'subplot_4']

for ax, (subplot_name, df) in zip(axes, all_codon_usages.groupby('subplot')):
    
# proteins = df['protein'].unique()
    # color = NRV_palette_custom.get(protein, '#333333')
    sns.stripplot(
        data=df,
        x='ro_codon',
        y='log2_comparison',
        palette=NRV_palette_custom,
        size=8,
        jitter=True,
        ax=ax,
        hue='protein'
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('$\log_{2}$(Final/Optimal)')
    ax.set_xlabel('')
    ax.axhline(y=0, color='black', linestyle='dotted', linewidth=1)
    if subplot_name == 'subplot_4':
        codon_labels = list(df['ro_codon'].unique())
        # Pad with empty strings to reach 16 labels
        while len(codon_labels) < 16:
            codon_labels.append('')
        ax.set_xticks(range(16))
        ax.set_xticklabels(codon_labels, rotation=0, fontsize=10)
    else:
        ax.legend(False)
# After plotting all subplots, set the same y-limits for all axes
# 1. Find the global min and max y-values across all subplot data
ymin = min(df['log2_comparison'].min() for df in subplot_groups.values())
ymax = max(df['log2_comparison'].max() for df in subplot_groups.values())

# Optionally, round or pad the limits for aesthetics
pad = 0.2 * (ymax - ymin)
ymin -= pad
ymax += pad

# 2. Set the same y-limits for all axes
for ax in axes:
    ax.set_ylim(ymin, ymax)

plt.savefig('codon-optimacy.svg', bbox_inches='tight', pad_inches=0)

plt.tight_layout()
plt.show()


