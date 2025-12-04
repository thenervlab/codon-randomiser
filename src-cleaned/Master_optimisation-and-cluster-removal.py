import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap
from loguru import logger
from itertools import groupby
from operator import itemgetter

np.random.seed(105105)

input_folder = 'sequences'
output_folder = 'results/Optimised-sequences'

# Iterate through all text files in the folder
new_seq = []
old_seq = []
for filename in os.listdir(input_folder):
    if filename.endswith('opt.txt'):  
        # Check if the file is a text file
        file_path = os.path.join(input_folder, filename)
        with open(file_path, 'r') as file:
            content = file.read()
            new_seq.append(content)
    elif filename.endswith('original.txt'):  
        # Check if the file is a text file
        file_path = os.path.join(input_folder, filename)
        with open(file_path, 'r') as file:
            content = file.read()
            old_seq.append(content)

# Break into codons
compared_codons = []
for orig1, opt2 in zip(old_seq, new_seq):
    orig_df = pd.DataFrame(wrap(orig1, 3), columns=['original_codon'])
    opt_df = pd.DataFrame(wrap(opt2, 3), columns=['optimised_codon'])
    combined_df = pd.merge(orig_df, opt_df, left_index=True, right_index=True)
    compared_codons.append(combined_df)

# Calculate the matching condons and nucleotides

def calculate_identity(codons):
    codons['codon_identity'] = [1 if old == new else 0 for old, new in codons[['original_codon', 'optimised_codon']].values]
    codon_identity = codons['codon_identity'].sum() / len(codons) * 100

    nt_identities = []
    for old, new in codons[['original_codon', 'optimised_codon']].values:
        if old == new:
            nt_identities.append(3)
        else:
            nt_identity = 0
            for nt_old, nt_new in zip(old, new):
                if nt_old == nt_new:
                    nt_identity += 1
        
            nt_identities.append(nt_identity)
    codons['nucleotide_identity'] = nt_identities
    
    return codons

codon_nt_identities = []
for i in compared_codons:
    nt_calc = calculate_identity(i)
    codon_nt_identities.append(nt_calc)

# Filter identical nucleotides
identical_nucs = []
for i in codon_nt_identities:
    ident_nucs = i[i['nucleotide_identity'] == 3].copy()
    identical_nucs.append(ident_nucs)

# Check for clusters of identical sequence - modify at least every second codon in a cluster using a random codon

# determine index numbers of identical nucleotide clusters

codon_clusters = []
for g in identical_nucs:
    clusters = [row for row in g.index]
    codon_clusters.append(clusters)

nt_clusters = []
for i in codon_clusters:
    seq_cluster = []
    for k, q in groupby(enumerate(i), lambda ix : ix[0] - ix[1]):
        seq_cluster.append((list(map(itemgetter(1), q))))
    nt_clusters.append(seq_cluster)

#exclude all the single values

filtered_clusters = [[cluster for cluster in sublist if len(cluster) >= 2] for sublist in nt_clusters]

# flatten list
flattened = [[item for sublist in outer_list for item in sublist] for outer_list in filtered_clusters]

# add another column to dodon table for yes/no to indicate when indices are part of clustered list

for i, x in zip(flattened, compared_codons):
    x['clustered?'] = [1 if g in i else 0 for g in x.index.tolist()]


# bring in codon table from randomiser and edit to read table rather then text

# Bring in codon table

input_codon_map = 'experimental_data/codon-table.csv'

codon_map = pd.read_csv(input_codon_map)
stop_codons = codon_map[codon_map['FullName'] == 'Stop'].copy()['Codon'].tolist() 

# utilise to make new randomised+optimised

stringency = 'high'


def randomise_codons(i):
    final_codons = []
    for codon, cluster in zip(codons['optimised_codon'], codons['clustered?']):
        if cluster==1:
            amino_acid = dict(codon_map[['Codon', 'AminoAcid']].values)[codon]

            available_codons = codon_map[codon_map['AminoAcid'] == amino_acid]['Codon'].tolist().copy()
            if stringency == 'high':
                if len(available_codons) > 1:
                    available_codons = [val for val in available_codons if val != codon]        
            final_codons.append(np.random.choice(available_codons))
        else:
            final_codons.append(codon)
    return final_codons

for codons in compared_codons:
    final_codons = randomise_codons(codons)
    codons['de-clustered_codons'] = final_codons

new_output_folder = 'final_sequences/'

sequence_order = []
#------------------------------------------RECALCULATE IDENTITIES------------------------------------

# Calculate the matching condons and nucleotides

def re_calculate_identity(codons):
    codons['de-clustered_codon_identity'] = [1 if old == new else 0 for old, new in codons[['original_codon', 'de-clustered_codons']].values]
    codon_identity = codons['de-clustered_codon_identity'].sum() / len(codons) * 100

    nt_identities = []
    for old, new in codons[['original_codon', 'de-clustered_codons']].values:
        if old == new:
            nt_identities.append(3)
        else:
            nt_identity = 0
            for nt_old, nt_new in zip(old, new):
                if nt_old == nt_new:
                    nt_identity += 1
        
            nt_identities.append(nt_identity)
    codons['de-clustered_nucleotide_identity'] = nt_identities
    
    return codons

final_codon_nt_identities = []
for i in compared_codons:
    nt_calc = re_calculate_identity(i)
    final_codon_nt_identities.append(nt_calc)

# save each dataframe in compared_codons as an individual CSV
for idx, df in enumerate(compared_codons, start=1):
    out_path = os.path.join(output_folder, f'codons_{idx}.csv')
    df.to_csv(out_path, index=False)



