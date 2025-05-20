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
    sequences.append(sequence)



# def plot_identity():
#     fig, ax = plt.subplots( figsize=(25, 20))
#     plt.stem(range(len(codons)), codons['nucleotide_identity'])
#     plt.stem(range(len(codons)), codons['graph_negatives'])
#     plt.xlabel('Position')
#     ax.set_yticks(range(-3,4))
#     ax.set_yticklabels((3, 2, 1, 0, 1, 2, 3))
#     return fig 

# plot_identity()

# plt.savefig('ENV.svg')

# Compare codon usage to optimal frequency
# -> add optimal frequencies to the codon-table.csv as a new column DONE

# -> calculate frequency usage in both GenScript and randomised sequences
    # final sequence
# frequency_count_final = {x: len(codons[codons['final_codon'] == x]) for x in codons['final_codon'].unique()}
# print(frequency_count_final)

# frequency_count_final = codons.groupby(['final_codon']).count()['final_codon_identity'].reset_index()


final_codons = codons['final_codon'].value_counts().reset_index()



# create a new column that has amino acids
    # use dictionary of codon table use map one column using dictionary to another set of information

input_codon_map = 'experimental_data/codon-table.csv'

codon_map = pd.read_csv(input_codon_map)

codon_map.sort_values(['AminoAcid'])


amino_dic = codon_map.set_index('Codon')['AminoAcid'].to_dict()


codons['aminoacid'] = codons['final_codon'].map(amino_dic)

# add column with amino acid
    # group together codon and amino acid
    # count function via groupby

proportion = (codons[['aminoacid', 'final_codon', 'final_nucleotide_identity']].groupby(['aminoacid', 'final_codon']).count() / codons[['aminoacid', 'final_nucleotide_identity']].groupby(['aminoacid']).count()).reset_index()

codon_frac_dic = codon_map.set_index('Codon')['Fraction'].to_dict()

proportion['optimal_frequency'] = proportion['final_codon'].map(codon_frac_dic)

proportion['optimal_freq-final'] = proportion['optimal_frequency']-proportion['final_nucleotide_identity']

# produce stacked bar chart


def plot_freq():
    fig, ax = plt.subplots()
    codon_map.groupby(['AminoAcid', 'Codon']) \
        ['Fraction'].sum() \
        .reset_index() \
        .pivot_table(index='AminoAcid', columns='Codon', values='Fraction') \
        .plot(kind='bar', stacked=True, ax=ax)
    ax.legend(title='Codon', bbox_to_anchor=(1.0, 1), loc='upper left')
    ax.set_ylabel('Fraction')
    ax.set_title('AminoAcid')
    for c in ax.containers:
        # Optional: if the segment is small or 0, customize the labels
        labels = [v.get_height() if v.get_height() > 0 else '' for v in c]
        # remove the labels parameter if it's not needed for customized labels
        ax.bar_label(c, labels=codon_map['Codon'], label_type='center')


plot_freq()

# identifying UG repeats
# Identify TG repeats and single As


proportion.to_csv(f'{output_folder}comparison-opt-ENV.csv')

# bring in codon table from randomiser and edit to read table rather then text

# Bring in codon table

input_codon_map = 'experimental_data/codon-table.csv'

codon_map = pd.read_csv(input_codon_map)
stop_codons = codon_map[codon_map['FullName'] == 'Stop'].copy()['Codon'].tolist() 

# utilise to make new randomised+optimised

stringency = 'high'

TG_free_final = []
for codon, cluster in zip(codons['final_codon'], codons['GTGT_clustered?']):
    if cluster==1:
        amino_acid = dict(codon_map[['Codon', 'AminoAcid']].values)[codon]

        available_codons = codon_map[codon_map['AminoAcid'] == amino_acid]['Codon'].tolist().copy()
        if stringency == 'high':
            if len(available_codons) > 1:
                available_codons = [val for val in available_codons if val != codon]        
        final_codons.append(np.random.choice(available_codons))
    else:
        final_codons.append(codon)
codons['final_codon_GTGT_free'] = final_codons  





# Visualise nucleotide identity

file_list = [filename for filename in os.listdir(input_path) if 'Optimised-sequences_' in filename]

sequences = []
for x in file_list:
    sequence = pd.read_csv(f'{input_path}{x}')
    sequence['graph_negatives'] = - sequence ['final_nucleotide_identity']
    sequence['data_source'] = x
    sequences.append(sequence)
sequences = pd.concat(sequences)

sequences['protein'] = sequences['data_source'].str.split('_').str[-1].str.split('.').str[0]

graph_pos = {'ENV':10,'POL':20,'PRO':30,'GAG':40}

sequences['adjusted_position_ro'] = sequences['protein'].map(graph_pos) + sequences['graph_negatives']

sequences['adjusted_position_o'] = sequences['protein'].map(graph_pos) + sequences['nucleotide_identity']


#Locate TG repeats in final codon

GTG_repeats = sequences[sequences['final_codon'] == 'GTG'].copy()
TGT_repeats = sequences[sequences['final_codon'] == 'TGT'].copy()


cluster_TGT = []
for row in TGT_repeats.index:
    cluster_TGT.append(row)

cluster_GTG = []
for row in GTG_repeats.index:
    cluster_GTG.append(row)

from itertools import groupby
from operator import itemgetter

GTG_cluster = []
for k, g in groupby(enumerate(cluster_GTG), lambda ix : ix[0] - ix[1]):
    GTG_cluster.append((list(map(itemgetter(1), g))))

TGT_cluster = []
for k, g in groupby(enumerate(cluster_TGT), lambda ix : ix[0] - ix[1]):
    TGT_cluster.append((list(map(itemgetter(1), g))))

for x in TGT_cluster:
  GTG_cluster.append(x)

total_cluster = GTG_cluster

clusters_3 = [x for x in total_cluster if len(x)>=3]

flattened_clusters = [item for sublist in clusters_3 for item in sublist]

sequences['GTGT_clustered?']= [1 if x in flattened_clusters else 0 for x in sequences.index.tolist()]


#TG repeats new_codon
GTG_repeats = sequences[sequences['new_codon'] == 'GTG'].copy()
TGT_repeats = sequences[sequences['new_codon'] == 'TGT'].copy()

cluster_TGT = []
for row in TGT_repeats.index:
    cluster_TGT.append(row)

cluster_GTG = []
for row in GTG_repeats.index:
    cluster_GTG.append(row)

from itertools import groupby
from operator import itemgetter

GTG_cluster = []
for k, g in groupby(enumerate(cluster_GTG), lambda ix : ix[0] - ix[1]):
    GTG_cluster.append((list(map(itemgetter(1), g))))

TGT_cluster = []
for k, g in groupby(enumerate(cluster_TGT), lambda ix : ix[0] - ix[1]):
    TGT_cluster.append((list(map(itemgetter(1), g))))

for x in TGT_cluster:
  GTG_cluster.append(x)

total_cluster = GTG_cluster

clusters_3 = [x for x in total_cluster if len(x)>=3]

flattened_clusters = [item for sublist in clusters_3 for item in sublist]

sequences['GTGT_clustered?_new_codon']= [1 if x in flattened_clusters else 0 for x in sequences.index.tolist()]

#Graphing

tick_labels = dict(sequences[['adjusted_position_ro', 'graph_negatives']].values)
tick_labels.update(sequences[['adjusted_position_o', 'nucleotide_identity']].values)

y_labels = ('ENV','POL','PRO','GAG')

col = np.where(sequences['GTGT_clustered?']==1,'k',np.where(sequences['GTGT_clustered?_new_codon'], 'r'))


fig, ax = plt.subplots( figsize=(35, 25))
for protein, df in sequences.groupby('protein'):
    plt.stem(range(len(df)), df['adjusted_position_ro'], bottom=graph_pos[protein])
    plt.stem(range(len(df)), df['adjusted_position_o'], bottom=graph_pos[protein])
plt.xlabel('Position',fontsize=30)
ax.set_yticks(list(tick_labels.keys()),labels=list(tick_labels.values()), fontsize=30)
ax.set_ylim(5,45)
ax.set_xticks(ax.get_xticks(), labels = [int(x) for x in ax.get_xticks()], fontsize=30)
ax.set_xlim(0,1000)


plt.savefig('ENV.svg')














# sum_aa = final_codons.groupby(['aminoacid']).sum()['count']
# proportion = final_codons.set_index('aminoacid')/sum_aa
# final_codons['frequency'] = 
# paired['codon_count'] = paired['final_codon'].map(frequency_count_final)
# paired.replace({"count":frequency_count_final}) 
# convert count for each of codons into the proportion
    # groupby by only amino acid finds how many of amino acid there are 
# logic is same for original barplot 

