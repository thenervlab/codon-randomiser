import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap

from loguru import logger

logger.info('Import OK')

input_path = 'experimental_data/test-seq1.txt'
input_codon_map = 'experimental_data/codon-table.csv'
output_folder = 'results/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Optional: Set seed
np.random.seed(10301)

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

# Locate stop codon (ideally at end, if missing add warning)
for stop in stop_codons:
    stop_pos = sequence.find(stop)
    if stop_pos > -1:
        logger.info(f'{stop} stop codon detected in position {stop_pos}.')
        break
if stop_pos == -1:
    logger.info(f'No stop codon detected.')

# For each codon in the sequence:
## Identify AA coded
## Identify options for other codons
## Randomly select a codon
## Apppend to new sequence

codons = pd.DataFrame(wrap(sequence[start_pos:stop_pos], 3)[:-1], columns=['original_codon'])

new_codons = []
for codon in codons['original_codon']:
    amino_acid = dict(codon_map[['Codon', 'AminoAcid']].values)[codon]
    
    available_codons = codon_map[codon_map['AminoAcid'] == amino_acid]['Codon'].tolist().copy()
    
    new_codons.append(np.random.choice(available_codons))
    
codons['new_codon'] = new_codons  


# Calculate % identity with original sequence
## At the codon level
## At the nucleotide level

codons['codon_identity'] = [1 if old == new else 0 for old, new in codons[['original_codon', 'new_codon']].values]

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
    

# Optional: Visualise comparison 
# -> ball & stick plot where each codon is a position,
# and each nucleotide within that codon that is matched
# scores +1 - ideally matched codons would not appear
# too many in a row



# Save (download) updated sequence and 
# randomisation report


