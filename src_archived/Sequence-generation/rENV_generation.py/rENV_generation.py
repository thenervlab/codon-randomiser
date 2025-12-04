import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap

from loguru import logger

logger.info('Import OK')

input_path = 'sequences/env_original.txt'
input_codon_map = 'experimental_data/codon-table.csv'
output_folder = 'results/ENV-data/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Optional: Set seed
np.random.seed(10301)

# Optional: Set stringency
stringency = 'high'
# stringency = 'low'


# Read in codon table
codon_map = pd.read_csv(input_codon_map)
stop_codons = codon_map[codon_map['FullName'] == 'Stop'].copy()['Codon'].tolist()

# Read in sequence file -  from txt or fasta?
sequence = pd.read_table(input_path).columns.tolist()[0].upper()

# Verify sequence is in groups of 3 codons
if len(sequence) % 3 != 0:
    logger.info(f'Sequence provided contains an incomplete codon. {len(sequence) % 3} nucleotides will be discarded.')
    sequence = sequence[:-(len(sequence) % 3)]

# Locate start codon (ideally at the start, if not find)
start_pos = sequence.find('ATG')
if start_pos == -1:
    logger.info(f'No start codon detected.')
else:
    logger.info(f'Start codon detected at pos. {start_pos}.')

# Locate stop codon (ideally at end, if missing add warning)
stop_posns = {}
from_start = pd.DataFrame(wrap(sequence[start_pos:], 3)[:-1], columns=['original_codon'])

for stop in stop_codons:
    stop_pos = from_start.index[from_start['original_codon'] == stop].min()
    if not np.isnan(stop_pos):
        stop_posns[stop_pos] = stop
if len(stop_posns) == 0:
    logger.info(f'No stop codon detected.')
    stop_pos = len(sequence) - start_pos  # Use the rest of the sequence
else:
    stop, stop_pos = stop_posns[np.min(list(stop_posns.keys()))], np.min(list(stop_posns.keys()))
    logger.info(f'{stop} stop codon detected in position {start_pos+stop_pos}.')
    

# For each codon in the sequence:
## Identify AA coded
## Identify options for other codons
## Randomly select a codon
## Apppend to new sequence

codons = pd.DataFrame(wrap(sequence[start_pos:start_pos+stop_pos], 3)[:-1], columns=['original_codon'])

new_codons = []
for codon in codons['original_codon']:
    amino_acid = dict(codon_map[['Codon', 'AminoAcid']].values)[codon]
    
    available_codons = codon_map[codon_map['AminoAcid'] == amino_acid]['Codon'].tolist().copy()
        
    if stringency == 'high':
        if len(available_codons) > 1:
            available_codons = [val for val in available_codons if val != codon]
    
    new_codons.append(np.random.choice(available_codons))
    
codons['new_codon'] = new_codons  


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

#NEW CODON TG REPEATS
#Locate TG repeats in optimised codon
new_codon_seq = ''.join(codons['new_codon'].tolist())
T_1 = re.sub(r'T(?=G)', '1', new_codon_seq)
T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
T0 = re.sub(r'[A-Z]', '0', T_1_A)
codons['T0_A_new'] = wrap(T0, 3)

# Replace '1010' repeats with 'RRRR'
repeat_seq = ''.join(codons['T0_A_new'].tolist())
pattern = r"(1010)"
matches = re.sub(pattern, 'RRRR', repeat_seq)
codons['identified_repeats_new'] = wrap(matches, 3)

# Mark codons with 'R'
codons['R_found_new'] = ['Y' if 'R' in repeats else 'N' for repeats in codons['identified_repeats_new']]


#identify 1st codon of RRRR repeats

codons['R_cluster'] = (codons['R_found_new'] == 'Y') & (codons['R_found_new'].shift(-1) == 'Y')

#Randomise RRRR repeats

input_codon_map = 'experimental_data/codon-table.csv'

codon_map = pd.read_csv(input_codon_map)

# Remove rows where 'Codon' starts with 'TG' and ends with 'TG'
codon_map_no_TG = codon_map[
    ~(codon_map['Codon'].str.startswith('TG') | codon_map['Codon'].str.endswith('TG'))
]

stop_codons = codon_map[codon_map['FullName'] == 'Stop'].copy()['Codon'].tolist() 

stringency = 'high'
new_codons = []

def randomise_codons(i):
    final_codons = []
    for codon, cluster in zip(codons['new_codon'], codons['R_cluster']):
        if cluster==True:
            amino_acid = dict(codon_map[['Codon', 'AminoAcid']].values)[codon]

            available_codons = codon_map[codon_map['AminoAcid'] == amino_acid]['Codon'].tolist().copy()
            if stringency == 'high':
                if len(available_codons) > 1: 
                    available_codons = [val for val in available_codons if val != codon]
                if len(available_codons) > 1:
                    if len([val for val in available_codons if 'TG' not in val]) > 0:
                        available_codons = [val for val in available_codons if 'TG' not in val]
                final_codons.append(np.random.choice(available_codons))
        else:
            final_codons.append(codon)
    return final_codons

codons['pre-ro_codon'] = randomise_codons(codons)

# re-calculate codon identity

nt_identities = []
for old, new in codons[['original_codon', 'pre-ro_codon']].values:
    if old == new:
        nt_identities.append(3)
    else:
        nt_identity = 0
        for nt_old, nt_new in zip(old, new):
            if nt_old == nt_new:
                nt_identity += 1
    
        nt_identities.append(nt_identity)
codons['final_nucleotide_identity'] = nt_identities

#--------------------------------COMPARISON OF rENV--------------------------------


# read in codon map
input_codon_map = 'experimental_data/codon-table.csv'

codon_map = pd.read_csv(input_codon_map)

codon_map.sort_values(['AminoAcid'])

amino_dic = codon_map.set_index('Codon')['AminoAcid'].to_dict()

codons['aminoacid'] = codons['pre-ro_codon'].map(amino_dic)

codons['protein'] = 'rENV'

#Calculate proportion of codon usage for original codon


codons_by_protein = [group.copy() for _, group in codons.groupby('protein')]

protein_proportion_new = []
for codons in codons_by_protein:
    # Get the protein name from the group
    protein_name = codons['protein'].iloc[0]
    proportion = (
        codons[['aminoacid', 'pre-ro_codon', 'final_nucleotide_identity']]
        .groupby(['aminoacid', 'pre-ro_codon'])
        .count() / codons[['aminoacid', 'final_nucleotide_identity']]
        .groupby(['aminoacid'])
        .count()
    ).reset_index()

    codon_frac_dic = codon_map.set_index('Codon')['Fraction'].to_dict()
    proportion['optimal_frequency'] = proportion['pre-ro_codon'].map(codon_frac_dic)
    proportion['optimal_freq-new'] = proportion['optimal_frequency'] - proportion['final_nucleotide_identity']
    proportion['protein'] = protein_name  # Add protein column
    protein_proportion_new.append(proportion)


protein_proportion_new = pd.concat(protein_proportion_new)


#------------------rEGFP protein codon usage------------------------
# comparison of all optimal frequencies
#drop stop codons
protein_proportion_new = protein_proportion_new[~protein_proportion_new['aminoacid'].str.contains('stp', case=False, na=False)]

protein_proportion_new['rand_RSCU'] = protein_proportion_new['final_nucleotide_identity']/protein_proportion_new['optimal_frequency']

new_codon_comparison = protein_proportion_new[['pre-ro_codon', 'aminoacid', 'rand_RSCU', 'protein','optimal_freq-new']].copy()

new_codon_comparison['sequence'] = 'rENV'

new_codon_comparison.to_csv('codon-usage-rENV-sequence.csv', index=False)
