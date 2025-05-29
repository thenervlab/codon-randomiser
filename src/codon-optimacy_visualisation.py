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

#TEST GLYCINE

# glycine_df = protein_codon_all[protein_codon_all['aminoacid'] == 'Gly'].copy()

# glycine_df['comparison'] = glycine_df['ro_nucleotide_identity']/glycine_df['optimal_frequency']

# glycine_df['log2_comparison'] = np.log2(glycine_df['comparison'])

# # Get unique proteins and assign a color to each
# proteins = glycine_df['protein'].unique()

# plt.figure(figsize=(5,5))

# # dot plot for ro_nucleotide_identity per codon
# for protein in proteins:
#     sub_df = glycine_df[glycine_df['protein'] == protein]
#     sns.stripplot(
#         data=sub_df,
#         x='ro_codon',
#         y='log2_comparison',
#         color=NRV_palette_custom[protein],
#         size=8,
#         jitter=False,
#         label=protein
#     )

# #dotted horizontal line at y=1
# plt.axhline(y=0, color='black', linestyle='dotted', linewidth=1)

# #edit legend
# handles = [mpatches.Patch(color=NRV_palette_custom[protein], label=protein) for protein in proteins]

# plt.legend(handles=handles, title='Protein')
# plt.title('Glycine')
# plt.xlabel('Codon')
# plt.ylabel('RO_freq/Opt_freq')
# plt.tight_layout()
# plt.show()

# comparison of all optimal frequencies
#drop stop codons
protein_codon_all = protein_codon_all[~protein_codon_all['aminoacid'].str.contains('stp', case=False, na=False)]

#calculate log2(ro/o)
protein_codon_all['comparison'] = protein_codon_all['ro_nucleotide_identity'] / protein_codon_all['optimal_frequency']

protein_codon_all['log2_comparison'] = np.log2(protein_codon_all['comparison'])

aminoacid_codon_usages = [group.copy() for _, group in protein_codon_all.groupby('aminoacid')]

#plot figure
from matplotlib.gridspec import GridSpec
import matplotlib.transforms as mtransforms
cm = 1/2.54  # centimeters in inches
fig, axes = plt.subplots(7,3, figsize=(18*cm, 42 * cm), layout='constrained')
axes = axes.flatten()

for idx, df in enumerate(aminoacid_codon_usages):
    ax = axes[idx]
    # Plot all proteins for this amino acid, colored by protein
    proteins = df['protein'].unique()
    for protein in proteins:
        sub_df = df[df['protein'] == protein]
        color = color = NRV_palette_custom.get(protein, '#333333')
        sns.stripplot(
            data=sub_df,
            x='ro_codon',
            y='log2_comparison',
            color=color,
            size=8,
            jitter=False,
            ax=ax
        )
    # Dotted horizontal line at y=0
    ax.axhline(y=0, color='black', linestyle='dotted', linewidth=1)
    # Title and labels
    ax.set_title(df['aminoacid'].iloc[0])
    ax.set_xlabel('')
    ax.set_ylabel('')

# Set a single x and y label for the whole figure
# fig.text(0.5, 0.04, 'Codon', ha='center', fontsize=12)
# fig.text(0.04, 0.5, 'log2(RO_freq/Opt_freq)', va='center', rotation='vertical', fontsize=12)

# Hide any unused subplots
for ax in axes[len(aminoacid_codon_usages):]:
    ax.axis('off')

handles = [mpatches.Patch(color=NRV_palette_custom.get(protein, '#333333'), label=protein) for protein in NRV_palette_custom]

# Place the legend in the first unused (hidden) subplot, if there is one
if len(aminoacid_codon_usages) < len(axes):
    legend_ax = axes[len(aminoacid_codon_usages)]
    handles = [mpatches.Patch(color=NRV_palette_custom.get(protein, '#333333'), label=protein) for protein in NRV_palette_custom]
    legend_ax.legend(
        handles=handles,
        title='Protein',
        loc='center',
        fontsize=18,
        title_fontsize=20
    )
    legend_ax.axis('off')

# Save the figure
plt.savefig('codon-optimacy-graphs.svg')

plt.show()
