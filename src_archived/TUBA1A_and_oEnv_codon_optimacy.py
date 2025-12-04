import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap

from loguru import logger

logger.info('Import OK')

input_path = 'experimental_data/TUBA1A.txt'
input_path_oENV = 'sequences/env_opt.txt'
input_codon_map = 'experimental_data/codon-table.csv'
output_folder = 'results/misc-codon-usage/'


# Read in sequence file -  from txt or fasta?
TUB_sequence = pd.read_table(input_path).columns.tolist()[0].upper()

TUB_split = [TUB_sequence[i:i+3] for i in range(0, len(TUB_sequence), 3)]


data = list(zip(TUB_split))
codons = pd.DataFrame(data, columns=['TUB_split'])
codons['protein'] = 'TUBA1A'

# read in codon map
input_codon_map = 'experimental_data/codon-table.csv'

codon_map = pd.read_csv(input_codon_map)

codon_map.sort_values(['AminoAcid'])

amino_dic = codon_map.set_index('Codon')['AminoAcid'].to_dict()

codons['aminoacid'] = codons['TUB_split'].map(amino_dic)

#Calculate proportion of codon usage for original codon

codons_by_protein = [group.copy() for _, group in codons.groupby('protein')]

protein_proportions_TUB = []
for codons in codons_by_protein:
    # Get the protein name from the group
    protein_name = codons['protein'].iloc[0]
    # Count codons per amino acid and codon
    codon_counts = codons.groupby(['aminoacid', 'TUB_split']).size().reset_index(name='codon_count')
    aa_counts = codons.groupby(['aminoacid']).size().reset_index(name='aa_count')
    # Merge to get proportion
    proportion = pd.merge(codon_counts, aa_counts, on='aminoacid')
    proportion['proportion'] = proportion['codon_count'] / proportion['aa_count']
    codon_frac_dic = codon_map.set_index('Codon')['Fraction'].to_dict()
    proportion['optimal_frequency'] = proportion['TUB_split'].map(codon_frac_dic)
    proportion['optimal_freq-original'] = proportion['optimal_frequency'] - proportion['proportion']
    proportion['protein'] = protein_name  # Add protein column
    protein_proportions_TUB.append(proportion)

protein_proportions_TUB = pd.concat(protein_proportions_TUB)
protein_proportions_TUB['sequence'] = 'TUBA1A'

#-----------------oEnv protein codon usage------------------------

# Read in sequence file
oENV_sequence = pd.read_table(input_path_oENV).columns.tolist()[0].upper()

oENV_split = [oENV_sequence[i:i+3] for i in range(0, len(oENV_sequence), 3)]


data = list(zip(oENV_split))
codons = pd.DataFrame(data, columns=['oENV_split'])
codons['protein'] = 'oENV'

# read in codon map
input_codon_map = 'experimental_data/codon-table.csv'

codon_map = pd.read_csv(input_codon_map)

codon_map.sort_values(['AminoAcid'])

amino_dic = codon_map.set_index('Codon')['AminoAcid'].to_dict()

codons['aminoacid'] = codons['oENV_split'].map(amino_dic)

#Calculate proportion of codon usage for original codon

codons_by_protein = [group.copy() for _, group in codons.groupby('protein')]

protein_proportions_oirginal_oENV = []
for codons in codons_by_protein:
    # Get the protein name from the group
    protein_name = codons['protein'].iloc[0]
    # Count codons per amino acid and codon
    codon_counts = codons.groupby(['aminoacid', 'oENV_split']).size().reset_index(name='codon_count')
    aa_counts = codons.groupby(['aminoacid']).size().reset_index(name='aa_count')
    # Merge to get proportion
    proportion = pd.merge(codon_counts, aa_counts, on='aminoacid')
    proportion['proportion'] = proportion['codon_count'] / proportion['aa_count']
    codon_frac_dic = codon_map.set_index('Codon')['Fraction'].to_dict()
    proportion['optimal_frequency'] = proportion['oENV_split'].map(codon_frac_dic)
    proportion['optimal_freq-original'] = proportion['optimal_frequency'] - proportion['proportion']
    proportion['protein'] = protein_name  # Add protein column
    protein_proportions_oirginal_oENV.append(proportion)

protein_proportions_oirginal_oENV = pd.concat(protein_proportions_oirginal_oENV)
protein_proportions_oirginal_oENV['sequence'] = 'oENV'

#------------------rEGFP protein codon usage------------------------
# comparison of all optimal frequencies


#concat both tables
concat_sequences = pd.concat([protein_proportions_oirginal_oENV, protein_proportions_TUB], ignore_index=True)

#drop stop codons
concat_sequences = concat_sequences[~concat_sequences['aminoacid'].str.contains('stp', case=False, na=False)]
#calculate RSCU
concat_sequences['RSCU'] = concat_sequences['proportion'] / concat_sequences['optimal_frequency']

concat_sequences.drop(columns=['codon_count', 'aa_count'], inplace=True)



concat_sequences.to_csv(f'{output_folder}codon-usage-TUBA1A-oENV.csv', index=False)
