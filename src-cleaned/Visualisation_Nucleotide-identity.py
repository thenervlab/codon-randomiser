import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import seaborn as sns
from textwrap import wrap
from loguru import logger
from itertools import groupby
from operator import itemgetter

logger.info('Import OK')

input_path = 'results/'
input_codon_map = 'experimental_data/codon-table.csv'
output_folder = 'results/'

np.random.seed(105105)

#Read un data
sequences = pd.read_csv('final_sequences/final_synonymous_sequences/all_sequences.csv')
sequences['protein'] = sequences['data_source'].str.split('_').str[-1].str.split('.').str[0]
sequences['protein'] = sequences['protein'].str.lower()

#-------------------------------Define graph positions--------------------------------
graph_pos = {'env':10,'pol':18,'pro':26,'gag':34}

#Original codon position
sequences['original_codon_position'] = 0
sequences['original_codon_adjusted_positions'] = sequences['protein'].map(graph_pos) + sequences['original_codon_position']

#Syn identity position
sequences['syn_graph_negatives'] = - sequences['syn_nucleotide_identity']
sequences['adjusted_position_syn'] = sequences['protein'].map(graph_pos) + sequences['syn_graph_negatives']

# Optimised identity position
sequences['adjusted_position_o'] = sequences['protein'].map(graph_pos) + sequences['nucleotide_identity']
tick_labels = dict(sequences[['adjusted_position_syn', 'syn_graph_negatives']].values)
tick_labels.update(sequences[['adjusted_position_o', 'nucleotide_identity']].values)
abs_tick_labels = {k:abs(v) for k,v in tick_labels.items()}

# Create a color palette 
NRV_palette_custom = {
    'pol': '#520B57',  
    'gag': '#80A4DA', 
    'pro': '#0BBDA3', 
    'env': '#FF960A',
}
stem_circle_size = 1
circle_callout_size = 2

#-------------------------------Visualisation--------------------------------
#plot graph
fig, ax = plt.subplots(figsize=(7, 4.375))
for protein, df in sequences.groupby('protein'):
    df = df.reset_index(drop=True)
    # Use custom palette for this protein
    color = NRV_palette_custom.get(protein, '#897FB8')
    markerline, stemlines, baseline = ax.stem(
        range(len(df)),
        df['adjusted_position_syn'],
        bottom=graph_pos[protein],
        linefmt=color,
        markerfmt=color,
        basefmt='k'
    )
    markerline.set_markersize(stem_circle_size)
    plt.setp(stemlines, 'linewidth', 0.5)

    markerline, stemlines, baseline = ax.stem(range(len(df)), df['adjusted_position_o'], bottom=graph_pos[protein], linefmt='#858585')
    markerline.set_markersize(stem_circle_size)
    plt.setp(stemlines, 'linewidth', 0.5)  
            
    markerline, stemlines, baseline = ax.stem(range(len(df)), df['original_codon_adjusted_positions'], bottom=graph_pos[protein], linefmt='#915A95')
    markerline.set_markersize(stem_circle_size)
    plt.setp(stemlines, 'linewidth', 0.5)

    # Clustered points: both R_found_syn and R_found_orig are 'Y'
    mask_both = (df['R_found_syn'] == 'Y') & (df['R_found_orig'] == 'Y')
    clustered_points = df[mask_both]
    markerline, stemlines, baseline = ax.stem(
        clustered_points.index,
        clustered_points['adjusted_position_o'],
        bottom=graph_pos[protein],
        linefmt='#A00935',
        markerfmt='#A00935',
        basefmt='k'
    )
    markerline.set_markersize(circle_callout_size)

    # Clustered points: R_found_syn is 'Y' and R_found_orig is not 'Y'
    mask_new_only = (df['R_found_syn'] == 'Y') & (df['R_found_orig'] != 'Y')
    clustered_points = df[mask_new_only]
    markerline, stemlines, baseline = ax.stem(
        clustered_points.index,
        clustered_points['adjusted_position_o'],
        bottom=graph_pos[protein],
        linefmt='#A00935',
        markerfmt='#A00935',
        basefmt='k'
    )
    markerline.set_markersize(circle_callout_size)
    baseline.set_color('k')

    clustered_points = df[df['R_found_syn'] == 'Y']
    markerline, stemlines, baseline = ax.stem(
        clustered_points.index,
        clustered_points['adjusted_position_syn'],
        bottom=graph_pos[protein],
        linefmt='#A00935',
        markerfmt='#A00935',
        basefmt='k'
    )
    markerline.set_markersize(circle_callout_size)

    clustered_points = df[df['R_found_orig'] == 'Y']
    markerline, stemlines, baseline = ax.stem(
        clustered_points.index,
        clustered_points['original_codon_adjusted_positions'],
        bottom=graph_pos[protein],
        linefmt='#292929',
        markerfmt= '#292929',
        basefmt='k'
    )
    markerline.set_color('#292929')
    markerline.set_markersize(3) 
    
plt.xlabel('Codon',fontsize=8)
plt.ylabel('ΔResidue identity', labelpad= +10)
ax.set_yticks(list(abs_tick_labels.keys()),labels=list(abs_tick_labels.values()), fontsize=8)
ax.set_ylim(5,40)
ax.set_xticks(ax.get_xticks(), labels = [int(x) for x in ax.get_xticks()], fontsize=8)
ax.set_xlim(0,1000)
# Remove the top and right borders
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

#Place labels at the y positions for each protein
for protein, ypos in graph_pos.items():
    trans = mtransforms.ScaledTranslation(-35/72, -11/72, fig.dpi_scale_trans)
    ax.text(
        +38, ypos +0.6,
        protein,
        transform=ax.transData + trans,
        fontsize=8, va='bottom',
        rotation=90
    )

plt.savefig('figures-cleaned/nucleotide_identity.svg', bbox_inches='tight', pad_inches=0.05)



