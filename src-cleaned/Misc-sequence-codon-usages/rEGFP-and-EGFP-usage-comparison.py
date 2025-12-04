import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap

from loguru import logger

logger.info('Import OK')

input_path = 'experimental_data/EGFP.txt'
input_path_rEGFP = 'experimental_data/rEGFP.txt'
input_codon_map = 'experimental_data/codon-table.csv'
output_folder = 'Results_cleaned/Misc-sequence-codon-usages'

# Read in sequence file
EGFP_sequence = pd.read_table(input_path).columns.tolist()[0].upper()
rEGFP_sequence = pd.read_table(input_path_rEGFP).columns.tolist()[0].upper()

#Split into codons
EGFP_split = [EGFP_sequence[i:i+3] for i in range(0, len(EGFP_sequence), 3)]
rEGFP_split = [rEGFP_sequence[i:i+3] for i in range(0, len(rEGFP_sequence), 3)]

# Produce dataframe
data = list(zip(EGFP_split, rEGFP_split))
codons = pd.DataFrame(data, columns=['original_codon', 'new_codon'])
codons['protein'] = 'EGFP'

# Calculate % identity with original sequence
## At the codon level
## At the nucleotide level

codons['codon_identity'] = [1 if old == new else 0 for old, new in codons[['original_codon', 'new_codon']].values]
codon_identity = codons['codon_identity'].sum() / len(codons) * 100
logger.info(f'The proportion of matched codons is:{codon_identity} %')

nt_identities = []
for old, new in codons[['original_codon', 'new_codon']].values:
    if old == new:
        nt_identities.append(3)
    else:
        nt_identity = 0
        for nt_old, nt_new in zip(old, new):
            if nt_old == nt_new:
                nt_identity += 1
    
        nt_identities.append(nt_identity)
codons['nucleotide_identity'] = nt_identities
    
#--------------------------------Comparison of EGFP--------------------------------
# read in codon map
input_codon_map = 'experimental_data/codon-table.csv'
codon_map = pd.read_csv(input_codon_map)
codon_map.sort_values(['AminoAcid'])
amino_dic = codon_map.set_index('Codon')['AminoAcid'].to_dict()
codons['aminoacid'] = codons['original_codon'].map(amino_dic)

#Calculate proportion of codon usage for original codon
codons_by_protein = [group.copy() for _, group in codons.groupby('protein')]

protein_proportions_oirginal = []
for codons in codons_by_protein:
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
    proportion['protein'] = protein_name 
    protein_proportions_oirginal.append(proportion)

protein_proportions_oirginal = pd.concat(protein_proportions_oirginal)

# Clean comparison of all optimal frequencies
# Drop stop codons
protein_proportions_oirginal = protein_proportions_oirginal[~protein_proportions_oirginal['aminoacid'].str.contains('stp', case=False, na=False)]
protein_proportions_oirginal['original_RSCU'] = (protein_proportions_oirginal['nucleotide_identity'] / protein_proportions_oirginal['optimal_frequency'])
original_codon_comparison = protein_proportions_oirginal[['original_codon', 'aminoacid', 'original_RSCU', 'protein']].copy()
original_codon_comparison['sequence'] = 'EGFP'

#--------------------------------Repeat for comparison of rEGFP--------------------------------
# Map amino acids for new codons
codons['aminoacid'] = codons['new_codon'].map(amino_dic)

#Calculate proportion of codon usage for original codon
codons_by_protein = [group.copy() for _, group in codons.groupby('protein')]

protein_proportion_new = []
for codons in codons_by_protein:
    # Get the protein name from the group
    protein_name = codons['protein'].iloc[0]
    proportion = (
        codons[['aminoacid', 'new_codon', 'nucleotide_identity']]
        .groupby(['aminoacid', 'new_codon'])
        .count() / codons[['aminoacid', 'nucleotide_identity']]
        .groupby(['aminoacid'])
        .count()
    ).reset_index()

    codon_frac_dic = codon_map.set_index('Codon')['Fraction'].to_dict()
    proportion['optimal_frequency'] = proportion['new_codon'].map(codon_frac_dic)
    proportion['optimal_freq-new'] = proportion['optimal_frequency'] - proportion['nucleotide_identity']
    proportion['protein'] = protein_name 
    protein_proportion_new.append(proportion)

protein_proportion_new = pd.concat(protein_proportion_new)

# Clean comparison of all optimal frequencies
# Drop stop codons
protein_proportion_new = protein_proportion_new[~protein_proportion_new['aminoacid'].str.contains('stp', case=False, na=False)]
protein_proportion_new['new_RSCU'] = (protein_proportion_new['nucleotide_identity'] / protein_proportion_new['optimal_frequency'])
new_codon_comparison = protein_proportion_new[['new_codon', 'aminoacid', 'new_RSCU', 'protein']].copy()
new_codon_comparison['sequence'] = 'rEGFP'

# Save combined data
concat_sequences= pd.concat([new_codon_comparison, original_codon_comparison], ignore_index=True)
concat_sequences.to_csv(f'{output_folder}/codon-usage-EGFP-sequences.csv', index=False)

