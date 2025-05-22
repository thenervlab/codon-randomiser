import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap
from loguru import logger

logger.info('Import OK')

input_path = 'results/'
input_codon_map = 'experimental_data/codon-table.csv'
output_folder = 'results/'

file_list = [filename for filename in os.listdir(input_path) if 'Optimised-sequences_' in filename]

sequences = []
for x in file_list:
    sequence = pd.read_csv(f'{input_path}{x}')
    sequence['graph_negatives'] = - sequence ['final_nucleotide_identity']
    sequence['data_source'] = x
    sequences.append(sequence)
sequences = pd.concat(sequences)

sequences['protein'] = sequences['data_source'].str.split('_').str[-1].str.split('.').str[0]

codon_map = pd.read_csv(input_codon_map)

codon_map.sort_values(['AminoAcid'])

amino_dic = codon_map.set_index('Codon')['AminoAcid'].to_dict()


sequences['aminoacid'] = sequences['final_codon'].map(amino_dic)






df_props = []
for protein in sequences['protein']:

    proportion = (sequences[['aminoacid', 'final_codon', 'final_nucleotide_identity']].groupby(['aminoacid', 'final_codon']).count() / sequences[['aminoacid', 'final_nucleotide_identity']].groupby(['aminoacid']).count()).reset_index()

    codon_frac_dic = codon_map.set_index('Codon')['Fraction'].to_dict()

    proportion['optimal_frequency'] = proportion['final_codon'].map(codon_frac_dic)

    proportion['optimal_freq-final'] = proportion['optimal_frequency']-proportion['final_nucleotide_identity']
    df_props.append(proportion)