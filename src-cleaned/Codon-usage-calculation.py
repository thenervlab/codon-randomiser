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
input_codon_map = 'experimental_data/codon-table.csv'
output_folder = 'Results_cleaned/'

# Read in sequences
sequences = pd.read_csv('final_sequences/final_synonymous_sequences/all_sequences.csv')
sequences['protein'] = sequences['data_source'].str.split('_').str[-1].str.split('.').str[0]
sequences['protein'] = sequences['protein'].str.lower()

# read in codon map
input_codon_map = 'experimental_data/codon-table.csv'
codon_map = pd.read_csv(input_codon_map)
codon_map.sort_values(['AminoAcid'])
amino_dic = codon_map.set_index('Codon')['AminoAcid'].to_dict()
sequences['aminoacid'] = sequences['syn_codon'].map(amino_dic)


sequences_by_protein = [group.copy() for _, group in sequences.groupby('protein')]


#-------------------------------Calculate codon usage proportions--------------------------------
protein_proportions = []
for codon in sequences_by_protein:
    # Get the protein name from the group
    protein_name = codon['protein'].iloc[0]
    proportion = (
        codon[['aminoacid', 'syn_codon', 'syn_nucleotide_identity']]
        .groupby(['aminoacid', 'syn_codon'])
        .count() / codon[['aminoacid', 'syn_nucleotide_identity']]
        .groupby(['aminoacid'])
        .count()
    ).reset_index()

    codon_frac_dic = codon_map.set_index('Codon')['Fraction'].to_dict()
    proportion['optimal_frequency'] = proportion['syn_codon'].map(codon_frac_dic)
    proportion['optimal_freq-ro'] = proportion['optimal_frequency'] - proportion['syn_nucleotide_identity']
    proportion['protein'] = protein_name  # Add protein column
    protein_proportions.append(proportion)

#save protein proportions

protein_proportions = pd.concat(protein_proportions)

protein_proportions.to_csv(f'{output_folder}codon-usage-syn.csv')


#-------------------------------Repeat for original codons--------------------------------

protein_proportions_original = []
for codons in sequences_by_protein:
    # Get the protein name from the group
    protein_name = codons['protein'].iloc[0]
    proportion = (
        codons[['aminoacid', 'original_codon', 'nucleotide_identity']]
        .groupby(['aminoacid', 'original_codon'])
        .count() / codons[['aminoacid', 'nucleotide_identity']]
        .groupby(['aminoacid'])
        .count()
    ).reset_index()

    codon_frac_dic = codon_map.set_index('Codon')['Fraction'].to_dict()
    proportion['optimal_frequency'] = proportion['original_codon'].map(codon_frac_dic)
    proportion['optimal_freq-original'] = proportion['optimal_frequency'] - proportion['nucleotide_identity']
    proportion['protein'] = protein_name  # Add protein column
    protein_proportions_original.append(proportion)

#save protein proportions

protein_proportions_original = pd.concat(protein_proportions_original)

protein_proportions_original.to_csv(f'{output_folder}codon-usage-original-sequence.csv')
