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

input_path = 'final_sequences/'
input_codon_map = 'experimental_data/codon-table.csv'
output_folder = 'results-cleaned/'

np.random.seed(105105)

file_list = [filename for filename in os.listdir(input_path) if 'optimised-randomised_' in filename]

sequences = []
for x in file_list:
    sequence = pd.read_csv(f'{input_path}{x}')
    # sequence['graph_negatives'] = - sequence['final_nucleotide_identity']
    sequence['data_source'] = x
    sequences.append(sequence)


#-------------------------------LOCATE TG REPEATS IN OPTIMISED CODON SEQUENCE------------------------------------
#Locate TG repeats in 'optimised_codon'
T0_list = []
for i in sequences: 
    optimised_codon_list = i['optimised_codon'].tolist()
    optimised_codon_seq = ''.join(optimised_codon_list)
    T_1 = re.sub(r'T(?=G)', '1', optimised_codon_seq)
    T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
    T0 = re.sub(r'[A-Z]', '0', T_1_A)
    i['T0_A_optimised'] = wrap(T0, 3)


#replace '10' repeats with R
min_repeats = 1
for i in sequences:
    repeat_lists = i['T0_A_optimised'].tolist()
    repeat_seq = ''.join(repeat_lists)
    pattern = f"(1010)"
    matches = re.sub(pattern, 'RRRR', repeat_seq)
    i['identified_repeats_optimised'] = wrap(matches, 3)

for i in sequences:
    i['R_found_optimised'] = ['Y' if 'R' in repeats else 'N' for repeats in i['identified_repeats_optimised']]

#-------------------------------LOCATE TG REPEATS IN DE-CLUSTERED CODON SEQUENCE---------------------------------
#Locate TG repeats in declustered codon
T0_list = []
for i in sequences: 
    declustered_codon_list = i['de-clustered_codons'].tolist()
    declustered_codon_seq = ''.join(declustered_codon_list)
    T_1 = re.sub(r'T(?=G)', '1', declustered_codon_seq)
    T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
    T0 = re.sub(r'[A-Z]', '0', T_1_A)
    i['T0_A_declustered'] = wrap(T0, 3)


#replace '10' repeats with R
min_repeats = 1
for i in sequences:
    repeat_lists = i['T0_A_declustered'].tolist()
    repeat_seq = ''.join(repeat_lists)
    pattern = f"(1010)"
    matches = re.sub(pattern, 'RRRR', repeat_seq)
    i['identified_repeats_declustered'] = wrap(matches, 3)

for i in sequences:
    i['R_found_declustered'] = ['Y' if 'R' in repeats else 'N' for repeats in i['identified_repeats_declustered']]

#-------------------------------LOCATE TG REPEATS IN ORIGINAL CODON SEQUENCE---------------------------------
#Locate TG repeats in original codon
T0_list = []
for i in sequences: 
    orig_codon_list = i['original_codon'].tolist()
    orig_codon_seq = ''.join(orig_codon_list)
    T_1 = re.sub(r'T(?=G)', '1', orig_codon_seq)
    T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
    T0 = re.sub(r'[A-Z]', '0', T_1_A)
    i['T0_A_orig'] = wrap(T0, 3)

#replace '10' repeats with R
min_repeats = 1
for i in sequences:
    repeat_lists = i['T0_A_orig'].tolist()
    repeat_seq = ''.join(repeat_lists)
    pattern = f"(1010)"
    matches = re.sub(pattern, 'RRRR', repeat_seq)
    i['identified_repeats_orig'] = wrap(matches, 3)

for i in sequences:
    i['R_found_orig'] = ['Y' if 'R' in repeats else 'N' for repeats in i['identified_repeats_orig']]

#----------------CONCATENATE SEQUENCE DATAFRAMES FOR TG REPEAT REMOVAL IN DE-CLUSTERED SEQUENCE-----------------

sequences = pd.concat(sequences)

#identify 1st codon of RRRR repeats
sequences['R_cluster'] = (sequences['R_found_declustered'] == 'Y') & (sequences['R_found_declustered'].shift(-1) == 'Y')

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
    for codon, cluster in zip(sequences['de-clustered_codons'], sequences['R_cluster']):
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

sequences['pre-syn_codon'] = randomise_codons(sequences)

#----------------------REVERSE TG REPEAT REMOVAL TO ORIGINAL CODON WHERE NEEDED-----------------------------

syn_codon_list = sequences['pre-syn_codon'].tolist()
syn_codon_seq = ''.join(syn_codon_list)
T_1 = re.sub(r'T(?=G)', '1', syn_codon_seq)
T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
T0 = re.sub(r'[A-Z]', '0', T_1_A)
sequences['T0_A_pre-syn'] = wrap(T0, 3)
    

#replace '10' repeats with R
repeat_lists = sequences['T0_A_pre-syn'].tolist()
repeat_seq = ''.join(repeat_lists)
pattern = f"(1010)"
matches = re.sub(pattern, 'RRRR', repeat_seq)
sequences['identified_repeats_pre-syn'] = wrap(matches, 3)
sequences['R_found_pre-syn'] = ['Y' if 'R' in repeats else 'N' for repeats in sequences['identified_repeats_pre-syn']]



# reinsert codon if TG repeat located in RO codon and not in original codon
sequences['syn_codon'] = [x if y == 'N' else z for x,y,z in zip (sequences['pre-syn_codon'], sequences['R_found_pre-syn'], sequences['original_codon'])]



#-------------------------RECALCULATE FINAL NUCLEOTIDE IDENTITY FOR RO_CODONS
sequences['syn_codon_identity'] = [1 if old == new else 0 for old, new in sequences[['optimised_codon', 'syn_codon']].values]
codon_identity = sequences['codon_identity'].sum() / len(sequences) * 100
logger.info(f'The proportion of matched codons is:{codon_identity} %')

nt_identities = []
for old, new in sequences[['original_codon', 'syn_codon']].values:
    if old == new:
        nt_identities.append(3)
    else:
        nt_identity = 0
        for nt_old, nt_new in zip(old, new):
            if nt_old == nt_new:
                nt_identity += 1
    
        nt_identities.append(nt_identity)
sequences['syn_nucleotide_identity'] = nt_identities


#------------------------IDENTIFY REMAINING TG REPEATS IN FINAL SYNONYMOUS CODON SEQUENECE---------------------
syn_codon_list = sequences['syn_codon'].tolist()
syn_codon_seq = ''.join(syn_codon_list)
T_1 = re.sub(r'T(?=G)', '1', syn_codon_seq)
T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
T0 = re.sub(r'[A-Z]', '0', T_1_A)
sequences['T0_A_syn'] = wrap(T0, 3)
    
#replace '10' repeats with R
repeat_lists = sequences['T0_A_syn'].tolist()
repeat_seq = ''.join(repeat_lists)
pattern = f"(1010)"
matches = re.sub(pattern, 'RRRR', repeat_seq)
sequences['identified_repeats_syn'] = wrap(matches, 3)
sequences['R_found_syn'] = ['Y' if 'R' in repeats else 'N' for repeats in sequences['identified_repeats_syn']]


#---------------------------------------SAVE SYN_CODON SEQUENCES----------------------------------------------

output_dir = "final_sequences/final_synonymous_sequences/"
os.makedirs(output_dir, exist_ok=True)

sequences.to_csv(os.path.join(output_dir, 'all_sequences.csv'), index=False)

grouped_sequences = {name: group.copy() for name, group in sequences.groupby('data_source')}

# Convert each group's 'syn_codon' column to a string
syn_codon_strings = {name: ''.join(group['syn_codon'].astype(str)) for name, group in grouped_sequences.items()}

# Save each syn_codon string as a text file

for name, codon_string in syn_codon_strings.items():
    # Create a simple filename based on the group name
    base_name = os.path.splitext(name)[0]  # Remove .csv extension
    base_name = base_name.replace("optimised-randomised", "final-synonymous")
    output_path = os.path.join(output_dir, f"{base_name}.txt")
    with open(output_path, "w") as f:
        f.write(codon_string)
