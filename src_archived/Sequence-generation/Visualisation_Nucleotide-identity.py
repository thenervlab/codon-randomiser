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

input_path = 'results/'
input_codon_map = 'experimental_data/codon-table.csv'
output_folder = 'results/'

np.random.seed(105105)

file_list = [filename for filename in os.listdir(input_path) if 'Optimised-sequences_' in filename]

sequences = []
for x in file_list:
    sequence = pd.read_csv(f'{input_path}{x}')
    sequence['graph_negatives'] = - sequence['final_nucleotide_identity']
    sequence['data_source'] = x
    sequences.append(sequence)

#NEW CODON TG REPEATS
#Locate TG repeats in optimised codon
T0_list = []
for i in sequences: 
    new_codon_list = i['new_codon'].tolist()
    new_codon_seq = ''.join(new_codon_list)
    T_1 = re.sub(r'T(?=G)', '1', new_codon_seq)
    T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
    T0 = re.sub(r'[A-Z]', '0', T_1_A)
    i['T0_A_new'] = wrap(T0, 3)


#replace '10' repeats with R
min_repeats = 1
for i in sequences:
    repeat_lists = i['T0_A_new'].tolist()
    repeat_seq = ''.join(repeat_lists)
    pattern = f"(1010)"
    matches = re.sub(pattern, 'RRRR', repeat_seq)
    i['identified_repeats_new'] = wrap(matches, 3)

for i in sequences:
    i['R_found_new'] = ['Y' if 'R' in repeats else 'N' for repeats in i['identified_repeats_new']]

#final CODON TG REPEATS
#Locate TG repeats in optimised codon
T0_list = []
for i in sequences: 
    new_codon_list = i['final_codon'].tolist()
    new_codon_seq = ''.join(new_codon_list)
    T_1 = re.sub(r'T(?=G)', '1', new_codon_seq)
    T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
    T0 = re.sub(r'[A-Z]', '0', T_1_A)
    i['T0_A_final'] = wrap(T0, 3)


#replace '10' repeats with R
min_repeats = 1
for i in sequences:
    repeat_lists = i['T0_A_final'].tolist()
    repeat_seq = ''.join(repeat_lists)
    pattern = f"(1010)"
    matches = re.sub(pattern, 'RRRR', repeat_seq)
    i['identified_repeats_final'] = wrap(matches, 3)

for i in sequences:
    i['R_found_final'] = ['Y' if 'R' in repeats else 'N' for repeats in i['identified_repeats_final']]

#Original codon TG repeats
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


#concatenate all dataframes
sequences = pd.concat(sequences)

#identify 1st codon of RRRR repeats

sequences['R_cluster'] = (sequences['R_found_final'] == 'Y') & (sequences['R_found_final'].shift(-1) == 'Y')

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
    for codon, cluster in zip(sequences['final_codon'], sequences['R_cluster']):
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

sequences['pre-ro_codon'] = randomise_codons(sequences)


#PRE-RO_CODON TG REPEATS
ro_codon_list = sequences['pre-ro_codon'].tolist()
ro_codon_seq = ''.join(ro_codon_list)
T_1 = re.sub(r'T(?=G)', '1', ro_codon_seq)
T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
T0 = re.sub(r'[A-Z]', '0', T_1_A)
sequences['T0_A_pre-ro'] = wrap(T0, 3)
    

#replace '10' repeats with R
repeat_lists = sequences['T0_A_pre-ro'].tolist()
repeat_seq = ''.join(repeat_lists)
pattern = f"(1010)"
matches = re.sub(pattern, 'RRRR', repeat_seq)
sequences['identified_repeats_pre-ro'] = wrap(matches, 3)
sequences['R_found_pre-ro'] = ['Y' if 'R' in repeats else 'N' for repeats in sequences['identified_repeats_pre-ro']]



# reinsert codon if TG repeat located in RO codon and not in original codon
sequences['ro_codon'] = [x if y == 'N' else z for x,y,z in zip (sequences['pre-ro_codon'], sequences['R_found_pre-ro'], sequences['original_codon'])]



#NUCLEOTIDE IDENTITY CALCULATION FOR RO_CODONS
sequences['ro_codon_identity'] = [1 if old == new else 0 for old, new in sequences[['new_codon', 'ro_codon']].values]
codon_identity = sequences['codon_identity'].sum() / len(sequences) * 100
logger.info(f'The proportion of matched codons is:{codon_identity} %')

nt_identities = []
for old, new in sequences[['original_codon', 'ro_codon']].values:
    if old == new:
        nt_identities.append(3)
    else:
        nt_identity = 0
        for nt_old, nt_new in zip(old, new):
            if nt_old == nt_new:
                nt_identity += 1
    
        nt_identities.append(nt_identity)
sequences['ro_nucleotide_identity'] = nt_identities

#RO_CODON TG REPEATS
ro_codon_list = sequences['ro_codon'].tolist()
ro_codon_seq = ''.join(ro_codon_list)
T_1 = re.sub(r'T(?=G)', '1', ro_codon_seq)
T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
T0 = re.sub(r'[A-Z]', '0', T_1_A)
sequences['T0_A_ro'] = wrap(T0, 3)
    
#replace '10' repeats with R
repeat_lists = sequences['T0_A_ro'].tolist()
repeat_seq = ''.join(repeat_lists)
pattern = f"(1010)"
matches = re.sub(pattern, 'RRRR', repeat_seq)
sequences['identified_repeats_ro'] = wrap(matches, 3)
sequences['R_found_ro'] = ['Y' if 'R' in repeats else 'N' for repeats in sequences['identified_repeats_ro']]

#Visualisation of nucleotide identity

sequences['protein'] = sequences['data_source'].str.split('_').str[-1].str.split('.').str[0]

graph_pos = {'ENV':10,'POL':18,'PRO':26,'GAG':34}


#Original codon position
sequences['original_codon_position'] = 0
sequences['original_codon_adjusted_positions'] = sequences['protein'].map(graph_pos) + sequences['original_codon_position']

#RO identity position
sequences['ro_graph_negatives'] = - sequences['ro_nucleotide_identity']
sequences['adjusted_position_ro'] = sequences['protein'].map(graph_pos) + sequences['ro_graph_negatives']

#O identity position
sequences['adjusted_position_o'] = sequences['protein'].map(graph_pos) + sequences['nucleotide_identity']

tick_labels = dict(sequences[['adjusted_position_ro', 'ro_graph_negatives']].values)
tick_labels.update(sequences[['adjusted_position_o', 'nucleotide_identity']].values)

abs_tick_labels = {k:abs(v) for k,v in tick_labels.items()}

# Create a color palette 
NRV_palette = [
    '#000000',
    '#141414', 
    '#525252',
    '#520B57', 
    '#915A95', 
    '#80A4DA', 
    '#0BBDA3', 
    '#FF960A',
    '#EA4A07',
    '#A00935'
]

NRV_palette_custom = {
    'POL': '#520B57',  
    'GAG': '#80A4DA', 
    'PRO': '#0BBDA3', 
    'ENV': '#FF960A',
}

stem_circle_size = 1
circle_callout_size = 2

#plot graph
fig, ax = plt.subplots(figsize=(7, 4.375))
for protein, df in sequences.groupby('protein'):
    
    # Use custom palette for this protein
    color = NRV_palette_custom.get(protein, '#897FB8')
    markerline, stemlines, baseline = ax.stem(
        range(len(df)),
        df['adjusted_position_ro'],
        bottom=graph_pos[protein],
        linefmt=color,
        markerfmt=color,
        basefmt='k'
    )
    markerline.set_markersize(stem_circle_size)
    plt.setp(stemlines, 'linewidth', 0.5)

    markerline, stemlines, baseline = ax.stem(range(len(df)), df['adjusted_position_o'], bottom=graph_pos[protein], linefmt='#858585')
    markerline.set_markersize(stem_circle_size)
    plt.setp(stemlines, 'linewidth', 0.5)  
            
    markerline, stemlines, baseline = ax.stem(range(len(df)), df['original_codon_adjusted_positions'], bottom=graph_pos[protein], linefmt='#915A95')
    markerline.set_markersize(stem_circle_size)  # Set marker size
    plt.setp(stemlines, 'linewidth', 0.5)        # Set line width (thinner

    # Clustered points: both R_found_new and R_found_orig are 'Y'
    mask_both = (df['R_found_new'] == 'Y') & (df['R_found_orig'] == 'Y')
    clustered_points = df[mask_both]
    markerline, stemlines, baseline = ax.stem(
        clustered_points.index,
        clustered_points['adjusted_position_o'],
        bottom=graph_pos[protein],
        linefmt='#A00935',  # Black line
        markerfmt='#A00935',  # Black marker
        basefmt='k'
    )
    markerline.set_markersize(circle_callout_size)

    # Clustered points: R_found_new is 'Y' and R_found_orig is not 'Y'
    mask_new_only = (df['R_found_new'] == 'Y') & (df['R_found_orig'] != 'Y')
    clustered_points = df[mask_new_only]
    markerline, stemlines, baseline = ax.stem(
        clustered_points.index,
        clustered_points['adjusted_position_o'],
        bottom=graph_pos[protein],
        linefmt='#A00935',  # Red line
        markerfmt='#A00935',  # Red marker
        basefmt='k'
    )
    markerline.set_markersize(circle_callout_size)
    baseline.set_color('k')

    clustered_points = df[df['R_found_ro'] == 'Y']
    markerline, stemlines, baseline = ax.stem(
        clustered_points.index,
        clustered_points['adjusted_position_ro'],
        bottom=graph_pos[protein],
        linefmt='#A00935',  # Red line
        markerfmt='#A00935',  # Red marker
        basefmt='k'
    )
    markerline.set_markersize(circle_callout_size)

    clustered_points = df[df['R_found_orig'] == 'Y']
    markerline, stemlines, baseline = ax.stem(
        clustered_points.index,
        clustered_points['original_codon_adjusted_positions'],
        bottom=graph_pos[protein],
        linefmt='#292929',  # Red line
        markerfmt= '#292929',  # Marker
        basefmt='k'
    )
    markerline.set_color('#292929')
    markerline.set_markersize(3) 
    

plt.xlabel('Codon',fontsize=8)
plt.ylabel('ΔResidue identity', labelpad= +10)
ax.set_yticks(list(abs_tick_labels.keys()),labels=list(abs_tick_labels.values()), fontsize=8)
ax.set_ylim(5,40)
ax.set_xticks(ax.get_xticks(), labels = [int(x) for x in ax.get_xticks()], fontsize=8)
ax.set_xlim(0,1000)
# Remove the top and right borders
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

#Place labels at the y positions for each protein
for protein, ypos in graph_pos.items():
    trans = mtransforms.ScaledTranslation(-35/72, -11/72, fig.dpi_scale_trans)
    ax.text(
        +38, ypos +0.6,  # x, y position (adjust as needed)
        protein,
        transform=ax.transData + trans,
        fontsize=8, va='bottom',
        rotation=90
    )

plt.savefig('nucleotide_identity.svg', bbox_inches='tight', pad_inches=0.05)


plt.savefig('nucleotide_identity_graph.svg')



#SAVE RO_CODON SEQUENCES

output_dir = "ro-sequences_no_motifs"
os.makedirs(output_dir, exist_ok=True)

sequences.to_csv(os.path.join(output_dir, 'all_sequences.csv'), index=False)

grouped_sequences = {name: group.copy() for name, group in sequences.groupby('data_source')}

# Convert each group's 'ro_codon' column to a string
ro_codon_strings = {name: ''.join(group['ro_codon'].astype(str)) for name, group in grouped_sequences.items()}



# Save each ro_codon string as a text file

for name, codon_string in ro_codon_strings.items():
    # Create a simple filename based on the group name
    base_name = os.path.splitext(name)[0]  # Remove .csv extension
    base_name = base_name.replace("Optimised", "ro")
    output_path = os.path.join(output_dir, f"{base_name}.txt")
    with open(output_path, "w") as f:
        f.write(codon_string)
