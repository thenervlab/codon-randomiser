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
    'POL': '#520B57',  
    'GAG': '#80A4DA', 
    'PRO': '#0BBDA3', 
    'ENV': '#FF960A',
}

# comparison of all optimal frequencies
#drop stop codons
protein_codon_all = protein_codon_all[~protein_codon_all['aminoacid'].str.contains('stp', case=False, na=False)]

#calculate log2(ro/o)
protein_codon_all['comparison'] = protein_codon_all['ro_nucleotide_identity'] / protein_codon_all['optimal_frequency']

protein_codon_all['log2_comparison'] = np.log2(protein_codon_all['comparison'])

aminoacid_codon_usages = [group.copy() for _, group in protein_codon_all.groupby('aminoacid')]

#Data sets for subplots
subplot_1_aa = ['Arg', 'Ala', 'Asn', 'Asp', 'Cys']
subplot_1 = [df for df in aminoacid_codon_usages if df['aminoacid'].iloc[0] in subplot_1_aa]
subplot_1_concat = pd.concat(subplot_1, ignore_index=True)

subplot_2_aa = ['Gln', 'Glu', 'Gly', 'His', 'Leu']
subplot_2 = [df for df in aminoacid_codon_usages if df['aminoacid'].iloc[0] in subplot_2_aa]
subplot_2_concat = pd.concat(subplot_2, ignore_index=True)

subplot_3_aa = ['Ile', 'Lys', 'Met', 'Pro', 'Ser']
subplot_3 = [df for df in aminoacid_codon_usages if df['aminoacid'].iloc[0] in subplot_3_aa]
subplot_3_concat = pd.concat(subplot_3, ignore_index=True)

subplot_4_aa = ['Phe', 'Thr', 'Trp', 'Tyr', 'Val']
subplot_4 = [df for df in aminoacid_codon_usages if df['aminoacid'].iloc[0] in subplot_4_aa]
subplot_4_concat = pd.concat(subplot_4, ignore_index=True)




subplot_concat_list = [
    ('Arg_Ala_Asn_Asp_Cys', subplot_1_concat),
    ('Gln_Glu_Gly_His_Leu', subplot_2_concat),
    ('Ile_Lys_Met_Pro_Ser', subplot_3_concat)
]

for fig_label, subplot_concat in subplot_concat_list:
    cm = 1/2.54  # centimeters in inches
    fig, ax = plt.subplots(figsize=(6.8, 2.3622))

    # Sort codons for consistent plotting
    protein_codon_all = subplot_concat.sort_values(['aminoacid', 'ro_codon'])

    # Create a categorical type for codons to ensure correct order
    codon_order = protein_codon_all['ro_codon'].unique()
    protein_codon_all['ro_codon'] = pd.Categorical(protein_codon_all['ro_codon'], categories=codon_order, ordered=True)

    # Plot all proteins, colored by protein
    proteins = protein_codon_all['protein'].unique()
    for protein in proteins:
        sub_df = protein_codon_all[protein_codon_all['protein'] == protein]
        color = NRV_palette_custom.get(protein, '#333333')
        sns.stripplot(
            data=sub_df,
            x='ro_codon',
            y='log2_comparison',
            color=color,
            size=8,
            jitter=True,
            ax=ax,
        )

    # Dotted horizontal line at y=0
    ax.axhline(y=0, color='black', linestyle='dotted', linewidth=1)

    # Set axis labels
    ax.set_xlabel('', fontsize=14)
    ax.set_ylabel('', fontsize=8)

    # Find the start and end indices for each amino acid group
    codon_to_aa = protein_codon_all.drop_duplicates('ro_codon').set_index('ro_codon')['aminoacid']
    amino_acids = [codon_to_aa[codon] for codon in codon_order]

    group_starts = []
    group_ends = []
    group_labels = []
    current_aa = amino_acids[0]
    start_idx = 0
    for idx, aa in enumerate(amino_acids + [None]):  # Add None to trigger last group
        if aa != current_aa:
            end_idx = idx - 1
            group_starts.append(start_idx)
            group_ends.append(end_idx)
            group_labels.append(current_aa)
            current_aa = aa
            start_idx = idx

    # Draw horizontal bars ("umbrellas") above each group, outside the graph area
    y_min, y_max = ax.get_ylim()
    y_bar = y_max + (y_max - y_min) * 0.05  # 8% above the top of the data
    for start, end in zip(group_starts, group_ends):
        ax.hlines(
            y=y_bar, xmin=start-0.4, xmax=end+0.4, color='black', linewidth=1, clip_on=False
        )

    # Adjust plot limits to make space for bars
    ax.set_ylim(y_min, y_bar + (y_max - y_min) * 0.12)
    ax.set_ylim(ax.get_ylim()[0], y_bar + 0.5)

    secax = ax.secondary_xaxis('top')
    secax.set_xticks([(s+e)/2 for s, e in zip(group_starts, group_ends)])
    secax.set_xticklabels(group_labels, rotation=0, fontsize=10)
    secax.set_xlabel('', fontsize=12)
    # Hide the top spine of the secondary axis robustly
    secax.spines['top'].set_color('none')
    secax.spines['top'].set_linewidth(0)

    plt.tight_layout()
    plt.savefig(f'codon-optimacy-graph-{fig_label}.svg', bbox_inches='tight')
    plt.close(fig)

















#legend plot
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(8, 2.15))
gs = gridspec.GridSpec(1, 2, width_ratios=[5, 0.5], wspace=0.01)  # even narrower right panel

# Main plot (squeezed columns)
ax = fig.add_subplot(gs[0])

protein_codon_all = subplot_4_concat.sort_values(['aminoacid', 'ro_codon'])
codon_order = protein_codon_all['ro_codon'].unique()
protein_codon_all['ro_codon'] = pd.Categorical(protein_codon_all['ro_codon'], categories=codon_order, ordered=True)

proteins = protein_codon_all['protein'].unique()
for protein in proteins:
    sub_df = protein_codon_all[protein_codon_all['protein'] == protein]
    color = NRV_palette_custom.get(protein, '#333333')
    sns.stripplot(
        data=sub_df,
        x='ro_codon',
        y='log2_comparison',
        color=color,
        size=8,           # smaller dots
        jitter=5,      # less jitter
        ax=ax,
        dodge=False       # no dodge, squeeze columns
    )

ax.axhline(y=0, color='black', linestyle='dotted', linewidth=1)
ax.set_xlabel('', fontsize=14)
ax.set_ylabel('', fontsize=8)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=9)

# Umbrella bars
codon_to_aa = protein_codon_all.drop_duplicates('ro_codon').set_index('ro_codon')['aminoacid']
amino_acids = [codon_to_aa[codon] for codon in codon_order]
group_starts, group_ends, group_labels = [], [], []
current_aa = amino_acids[0]
start_idx = 0
for idx, aa in enumerate(amino_acids + [None]):
    if aa != current_aa:
        end_idx = idx - 1
        group_starts.append(start_idx)
        group_ends.append(end_idx)
        group_labels.append(current_aa)
        current_aa = aa
        start_idx = idx
y_min, y_max = ax.get_ylim()
y_bar = y_max + (y_max - y_min) * 0.05
for start, end in zip(group_starts, group_ends):
    ax.hlines(y=y_bar, xmin=start-0.4, xmax=end+0.4, color='black', linewidth=1, clip_on=False)
ax.set_ylim(y_min, y_bar + (y_max - y_min) * 0.12)
ax.set_ylim(ax.get_ylim()[0], y_bar + 0.5)

secax = ax.secondary_xaxis('top')
secax.set_xticks([(s+e)/2 for s, e in zip(group_starts, group_ends)])
secax.set_xticklabels(group_labels, rotation=0, fontsize=9)
secax.set_xlabel('', fontsize=12)
secax.spines['top'].set_color('none')
secax.spines['top'].set_linewidth(0)

# Right panel for legend
ax_right = fig.add_subplot(gs[1])
ax_right.axis('off')
handles = [
    mlines.Line2D([], [], color=NRV_palette_custom.get(protein, '#333333'), marker='o', linestyle='None', markersize=8, label=protein)
    for protein in proteins
]
ax_right.legend(
    handles=handles,
    title='Protein',
    loc='center',
    fontsize=9,
    title_fontsize=10,
    frameon=False,
    bbox_to_anchor=(0.5, 0.5),
    borderaxespad=0.0,
    handletextpad=0.2,
    labelspacing=0.2,
    borderpad=0.1
)

plt.savefig('codon-optimacy-legend.svg', bbox_inches='tight', pad_inches=0)
plt.show()

