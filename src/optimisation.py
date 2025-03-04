import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap
from loguru import logger


old_seq = 'ATGGGGCAAACTAAAAGTAAAATTAAAAGTAAATATGCCTCTTATCTCAGCTTTATTAAAATTCTTTTAAAAAGAGGGGGAGTTAAAGTATCTACAAAAAATCTAATCAAGCTATTTCAAATAATAGAACAATTTTGCCCATGGTTTCCAGAACAAGGAACTTTAGATCTAAAAGATTGGAAAAGAATTGGTAAGGAACTAAAACAAGCAGGTAGGAAGGGTAATATCATTCCACTTACAGTATGGAATGATTGGGCCATTATTAAAGCAGCTTTAGAACCATTTCAAACAGAAGAAGATAGCGTTTCAGTTTCTGATGCCCCTGGAAGCTGTATAATAGATTGTAATGAAAACACAAGGAAAAAATCCCAGAAAGAAACGGAAGGTTTACATTGCGAATATGTAGCAGAGCCGGTAATGGCTCAGTCAACGCAAAATGTTGACTATAATCAATTACAGGAGGTGATATATCCTGAAACGTTAAAATTAGAAGGAAAAGGTCCAGAATTAGTGGGGCCATCAGAGTCTAAACCACGAGGCACAAGTCCTCTTCCAGCAGGTCAGGTGCCCGTAACATTACAACCTCAAAAGCAGGTTAAAGAAAATAAGACCCAACCGCCAGTAGCCTATCAATACTGGCCTCCGGCTGAACTTCAGTATCGGCCACCCCCAGAAAGTCAGTATGGATATCCAGGAATGCCCCCAGCACCACAGGGCAGGGCGCCATACCCTCAGCCGCCCACTAGGAGACTTAATCCTACGGCACCACCTAGTAGACAGGGTAGTGAATTACATGAAATTATTGATAAATCAAGAAAGGAAGGAGATACTGAGGCATGGCAATTCCCAGTAACGTTAGAACCGATGCCACCTGGAGAAGGAGCCCAAGAGGGAGAGCCTCCCACAGTTGAGGCCAGATACAAGTCTTTTTCGATAAAAATGCTAAAAGATATGAAAGAGGGAGTAAAACAGTATGGACCCAACTCCCCTTATATGAGGACATTATTAGATTCCATTGCTCATGGACATAGACTCATTCCTTATGATTGGGAGATTCTGGCAAAATCGTCTCTCTCACCCTCTCAATTTTTACAATTTAAGACTTGGTGGATTGATGGGGTACAAGAACAGGTCCGAAGAAATAGGGCTGCCAATCCTCCAGTTAACATAGATGCAGATCAACTATTAGGAATAGGTCAAAATTGGAGTACTATTAGTCAACAAGCATTAATGCAAAATGAGGCCATTGAGCAAGTTAGAGCTATCTGCCTTAGAGCCTGGGAAAAAATCCAAGACCCAGGAAGTACCTGCCCCTCATTTAATACAGTAAGACAAGGTTCAAAAGAGCCCTATCCTGATTTTGTGGCAAGGCTCCAAGATGTTGCTCAAAAGTCAATTGCCGATGAAAAAGCCCGTAAGGTCATAGTGGAGTTGATGGCATATGAAAACGCCAATCCTGAGTGTCAATCAGCCATTAAGCCATTAAAAGGAAAGGTTCCTGCAGGATCAGATGTAATCTCAGAATATGTAAAAGCCTGTGATGGAATCGGAGGAGCTATGCATAAAGCTATGCTTATGGCTCAAGCAATAACAGGAGTTGTTTTAGGAGGACAAGTTAGAACATTTGGAGGAAAATGTTATAATTGTGGTCAAATTGGTCACTTAAAAAAGAATTGCCCAGTCTTAAACAAACAGAATATAACTATTCAAGCAACTACAACAGGTAGAGAGCCACCTGACTTATGTCCAAGATGTAAAAAAGGAAAACATTGGGCTAGTCAATGTCGTTCTAAATTTGATAAAAATGGGCAACCATTGTCGGGAAACGAGCAAAGGGGCCAGCCTCAGGCCCCACAACAAACTGGGGCATTCCCAATTCAGCCATTTGTTCCTCAGGGTTTTCAGGGACAACAACCCCCACTGTCCCAAGTGTTTCAGGGAATAAGCCAGTTACCACAATACAACAATTGTCCCCCGCCACAAGCGGCAGTGCAGCAGTAG'
new_seq = 'ATGGGCCAGACCAAGAGCAAGATCAAAAGCAAGTACGCCAGCTACCTGAGCTTCATCAAGATCCTGCTGAAGAGAGGCGGAGTGAAAGTGAGCACCAAGAATCTGATCAAGCTCTTCCAGATCATCGAGCAGTTTTGCCCCTGGTTCCCTGAGCAGGGCACCCTGGATCTGAAGGACTGGAAGCGGATCGGCAAAGAACTGAAGCAGGCCGGCCGGAAGGGCAACATCATCCCTCTTACAGTGTGGAACGACTGGGCCATCATCAAGGCCGCCCTGGAACCTTTCCAGACAGAGGAAGATAGCGTGTCTGTCAGCGATGCCCCTGGAAGCTGCATCATTGATTGCAACGAAAACACCAGAAAGAAGAGTCAGAAGGAAACCGAGGGCCTGCACTGCGAGTACGTGGCAGAACCTGTCATGGCCCAGAGCACCCAGAACGTGGACTACAACCAGCTGCAGGAGGTGATTTACCCCGAGACTCTGAAGCTGGAAGGTAAAGGCCCCGAACTGGTGGGCCCTTCTGAAAGCAAGCCGCGGGGCACCAGCCCCCTGCCTGCTGGCCAAGTGCCCGTGACCCTCCAACCTCAAAAGCAGGTCAAGGAAAACAAGACCCAGCCCCCTGTGGCCTATCAGTACTGGCCTCCCGCCGAGCTACAGTATAGACCTCCACCTGAGAGCCAGTACGGATACCCCGGCATGCCTCCAGCTCCTCAAGGACGCGCCCCTTACCCTCAGCCTCCCACCAGAAGACTGAACCCCACCGCCCCTCCAAGCAGACAGGGCTCTGAACTGCACGAGATCATCGATAAGTCCCGGAAGGAGGGCGACACCGAGGCCTGGCAGTTCCCCGTGACATTGGAACCTATGCCCCCCGGCGAGGGCGCCCAGGAGGGAGAGCCTCCTACAGTTGAGGCCAGATACAAGAGCTTCTCCATCAAGATGCTGAAAGATATGAAAGAAGGCGTGAAACAGTACGGCCCTAACAGCCCTTACATGAGAACACTGCTGGACTCTATCGCCCACGGCCACAGACTGATCCCCTACGACTGGGAAATCCTGGCCAAGTCTTCTCTCTCTCCTAGTCAGTTCCTGCAGTTCAAGACCTGGTGGATCGACGGCGTCCAGGAGCAGGTGCGGCGGAACAGAGCCGCTAATCCTCCTGTGAACATCGACGCCGACCAGCTGCTGGGCATCGGCCAGAATTGGAGCACCATCAGCCAGCAGGCCCTGATGCAGAATGAGGCCATCGAGCAGGTGAGGGCCATCTGCCTGAGAGCATGGGAGAAAATCCAAGACCCCGGCAGCACATGCCCTAGCTTCAACACCGTGCGGCAGGGCAGCAAGGAACCATATCCTGACTTCGTGGCCCGGCTGCAGGACGTGGCCCAGAAGTCCATCGCCGATGAGAAGGCCAGAAAAGTGATCGTGGAACTGATGGCTTACGAGAACGCCAATCCTGAGTGCCAGTCCGCCATCAAGCCCCTGAAAGGCAAGGTGCCAGCTGGCAGCGATGTGATCTCTGAGTACGTGAAGGCCTGCGACGGAATCGGCGGTGCTATGCACAAGGCTATGCTGATGGCCCAGGCCATTACAGGAGTGGTGCTGGGCGGCCAGGTACGGACCTTCGGCGGAAAGTGCTACAACTGTGGCCAGATCGGCCATCTGAAGAAAAACTGTCCCGTGCTGAACAAGCAAAACATCACAATCCAGGCTACAACCACAGGCAGAGAGCCTCCAGACCTGTGCCCAAGATGTAAAAAGGGCAAGCACTGGGCTAGCCAGTGTAGAAGCAAATTCGACAAGAATGGTCAGCCTCTGTCCGGCAACGAGCAAAGAGGCCAGCCACAAGCCCCCCAACAGACCGGCGCTTTTCCTATCCAGCCCTTCGTGCCTCAGGGATTTCAGGGCCAGCAGCCTCCTCTGAGCCAGGTGTTTCAGGGAATTTCTCAGCTGCCCCAATACAACAACTGCCCTCCTCCCCAGGCTGCCGTGCAGCAGTGA'

# Break into codons
old_codons = pd.DataFrame(wrap(old_seq, 3), columns=['original_codon'])
new_codons = pd.DataFrame(wrap(new_seq, 3), columns=['new_codon'])

codons = pd.merge(old_codons, new_codons, left_index=True, right_index=True)

# Calculate the matching condons and nucleotides

def calculate_identity(codons):
    codons['codon_identity'] = [1 if old == new else 0 for old, new in codons[['original_codon', 'new_codon']].values]
    codon_identity = codons['codon_identity'].sum() / len(codons) * 100

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
    
    return codons

codons = calculate_identity(codons)

# Visualise overlap in old and new sequences
def plot_identity():
    fig, ax = plt.subplots(figsize=(20, 5))
    plt.stem(range(len(codons)), codons['nucleotide_identity'])
    plt.xlabel('Position')
    
    return fig

plot_identity()

# Filter identical nucleotides
identical_nucs = codons[codons['nucleotide_identity'] == 3].copy()

# Check for clusters of identical sequence - modify at least every second codon in a cluster using a random codon

# determine index numbers of identical nucleotide clusters

cluster_nucs = []
for row in identical_nucs.index:
    cluster_nucs.append(row)

from itertools import groupby
from operator import itemgetter
    

nucleotide_clus = []
for k, g in groupby(enumerate(cluster_nucs), lambda ix : ix[0] - ix[1]):
    nucleotide_clus.append((list(map(itemgetter(1), g))))


#exclude all the single values

clusters = [x for x in nucleotide_clus if len(x)>=2]

# flatten list 

flattened_nucs = [item for sublist in clusters for item in sublist]

# add another column to dodon table for yes/no to indicate when indices are part of clustered list

codons['clustered?']= [1 if x in flattened_nucs else 0 for x in codons.index.tolist()]



# bring in codon table from randomiser and edit to read table rather then text

# Bring in codon table

input_codon_map = 'experimental_data/codon-table.csv'

codon_map = pd.read_csv(input_codon_map)
stop_codons = codon_map[codon_map['FullName'] == 'Stop'].copy()['Codon'].tolist() 

# utilise to make new randomised+optimised

stringency = 'high'

final_codons = []
for codon, cluster in zip(codons['new_codon'], codons['clustered?']):
    if cluster==1:
        amino_acid = dict(codon_map[['Codon', 'AminoAcid']].values)[codon]

        available_codons = codon_map[codon_map['AminoAcid'] == amino_acid]['Codon'].tolist().copy()
        if stringency == 'high':
            if len(available_codons) > 1:
                available_codons = [val for val in available_codons if val != codon]        
        final_codons.append(np.random.choice(available_codons))
    else:
        final_codons.append(codon)
codons['final_codon'] = final_codons  


# Recalculate identity with updated codons

codons['final_codon_identity'] = [1 if old == new else 0 for old, new in codons[['new_codon', 'final_codon']].values]
codon_identity = codons['codon_identity'].sum() / len(codons) * 100
logger.info(f'The proportion of matched codons is:{codon_identity} %')

nt_identities = []
for old, new in codons[['new_codon', 'final_codon']].values:
    if old == new:
        nt_identities.append(3)
    else:
        nt_identity = 0
        for nt_old, nt_new in zip(old, new):
            if nt_old == nt_new:
                nt_identity += 1
    
        nt_identities.append(nt_identity)
codons['final_nucleotide_identity'] = nt_identities

# Visualise overlap in old and new sequences
def plot_identity():
    fig, ax = plt.subplots(figsize=(20, 5))
    plt.stem(range(len(codons)), codons['final_nucleotide_identity'])
    plt.xlabel('Position')
    
    return fig

plot_identity()
# Compare codon usage to optimal frequency
# -> add optimal frequencies to the codon-table.csv as a new column

# -> calculate frequency usage in both GenScript and randomised sequences
# -> Visualise difference in optimal for each amino acid 