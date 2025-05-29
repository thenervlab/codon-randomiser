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

input_path = 'ro-sequences_no_motifs/'
sequence_path = 'sequences/'
input_codon_map = 'experimental_data/codon-table.csv'
output_folder = 'results/'

np.random.seed(105105)
#read in original codons
sequences = pd.read_csv('ro-sequences_no_motifs/all_sequences.csv')
sequences.drop(
    columns=[
        'Unnamed: 0', 'R_found_orig', 'graph_negatives', 'T0_A_new', 'identified_repeats_new',
        'R_found_new', 'T0_A_final', 'identified_repeats_final', 'R_found_final', 'T0_A_orig',
        'identified_repeats_orig', 'R_found_orig', 'R_cluster', 'pre-ro_codon', 'T0_A_pre-ro', 'identified_repeats_pre-ro','R_found_pre-ro','T0_A_ro','identified_repeats_ro','R_found_ro','original_codon_position','original_codon_adjusted_positions','ro_graph_negatives','adjusted_position_ro','adjusted_position_o'
    ],
    inplace=True
)

codons = sequences

# read in codon map
input_codon_map = 'experimental_data/codon-table.csv'

codon_map = pd.read_csv(input_codon_map)

codon_map.sort_values(['AminoAcid'])


amino_dic = codon_map.set_index('Codon')['AminoAcid'].to_dict()


codons['aminoacid'] = codons['ro_codon'].map(amino_dic)

# add column with amino acid
    # group together codon and amino acid
    # count function via groupby



codons_by_protein = [group.copy() for _, group in codons.groupby('protein')]

protein_proportions = []
for codons in codons_by_protein:
    # Get the protein name from the group
    protein_name = codons['protein'].iloc[0]
    proportion = (
        codons[['aminoacid', 'ro_codon', 'ro_nucleotide_identity']]
        .groupby(['aminoacid', 'ro_codon'])
        .count() / codons[['aminoacid', 'ro_nucleotide_identity']]
        .groupby(['aminoacid'])
        .count()
    ).reset_index()

    codon_frac_dic = codon_map.set_index('Codon')['Fraction'].to_dict()
    proportion['optimal_frequency'] = proportion['ro_codon'].map(codon_frac_dic)
    proportion['optimal_freq-ro'] = proportion['optimal_frequency'] - proportion['ro_nucleotide_identity']
    proportion['protein'] = protein_name  # Add protein column
    protein_proportions.append(proportion)

#save protein proportions

protein_proportions = pd.concat(protein_proportions)

protein_proportions.to_csv(f'{output_folder}codon-usage.csv')