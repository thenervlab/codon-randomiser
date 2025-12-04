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


input_path = 'results/Optimised-sequences_ENV.csv'

#checking LTR


sequences = pd.read_csv(input_path)


#original TG REPEATS
ro_codon_list = sequences['original_codon'].tolist()
ro_codon_seq = ''.join(ro_codon_list)
T_1 = re.sub(r'T(?=G)', '1', ro_codon_seq)
T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
T0 = re.sub(r'[A-Z]', '0', T_1_A)
sequences['T0_A_original'] = wrap(T0, 3)
    

#replace '10' repeats with R
repeat_lists = sequences['T0_A_original'].tolist()
repeat_seq = ''.join(repeat_lists)
pattern = f"(1010)"
matches = re.sub(pattern, 'RRRR', repeat_seq)
sequences['identified_repeats_original'] = wrap(matches, 3)
sequences['R_found_original'] = ['Y' if 'R' in repeats else 'N' for repeats in sequences['identified_repeats_original']]

#final TG REPEATS
ro_codon_list = sequences['final_codon'].tolist()
ro_codon_seq = ''.join(ro_codon_list)
T_1 = re.sub(r'T(?=G)', '1', ro_codon_seq)
T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
T0 = re.sub(r'[A-Z]', '0', T_1_A)
sequences['T0_A_final'] = wrap(T0, 3)
    

#replace '10' repeats with R
repeat_lists = sequences['T0_A_final'].tolist()
repeat_seq = ''.join(repeat_lists)
pattern = f"(1010)"
matches = re.sub(pattern, 'RRRR', repeat_seq)
sequences['identified_repeats_final'] = wrap(matches, 3)
sequences['R_found_final'] = ['Y' if 'R' in repeats else 'N' for repeats in sequences['identified_repeats_final']]


#new TG REPEATS
ro_codon_list = sequences['new_codon'].tolist()
ro_codon_seq = ''.join(ro_codon_list)
T_1 = re.sub(r'T(?=G)', '1', ro_codon_seq)
T_1_A = re.sub(r'(?=10)A(?=10)|(?=01)A(?=01)', '2', T_1)
T0 = re.sub(r'[A-Z]', '0', T_1_A)
sequences['T0_A_new'] = wrap(T0, 3)
    

#replace '10' repeats with R
repeat_lists = sequences['T0_A_new'].tolist()
repeat_seq = ''.join(repeat_lists)
pattern = f"(1010)"
matches = re.sub(pattern, 'RRRR', repeat_seq)
sequences['identified_repeats_new'] = wrap(matches, 3)
sequences['R_found_new'] = ['Y' if 'R' in repeats else 'N' for repeats in sequences['identified_repeats_new']]

#Original codon position
sequences['original_codon_position'] = 0

#final nucleotide identity position
sequences['final_graph_negatives'] = - sequences['final_nucleotide_identity']


#plot graph
fig, ax = plt.subplots(figsize=(35, 25))
ax.stem(range(len(sequences)), sequences['final_graph_negatives'])
ax.stem(range(len(sequences)), sequences['nucleotide_identity'])
ax.stem(range(len(sequences)), sequences['original_codon_position'])
clustered_points = sequences[sequences['R_found_new'] == 'Y']
ax.stem(
    clustered_points.index,
    clustered_points['nucleotide_identity'],
    linefmt='r-',  # Red line
    markerfmt='ro',  # Red marker
    basefmt=' '
)
clustered_points = sequences[sequences['R_found_final'] == 'Y']
ax.stem(
    clustered_points.index,
    clustered_points['final_graph_negatives'],
    linefmt='r-',  # Red line
    markerfmt='ro',  # Red marker
    basefmt=' '
)
clustered_points = sequences[sequences['R_found_original'] == 'Y']
markerline, stemlines, baseline = ax.stem(
    clustered_points.index,
    clustered_points['original_codon_position'],
    linefmt='r-',  # Red line
    markerfmt='o',  # Marker
    basefmt=' '
)
markerline.set_color('yellow')

# 
plt.xlabel('Position',fontsize=30)
ax.set_yticks(list(tick_labels.keys()),labels=list(tick_labels.values()), fontsize=30)
ax.set_ylim(5,45)
ax.set_xticks(ax.get_xticks(), labels = [int(x) for x in ax.get_xticks()], fontsize=30)
ax.set_xlim(0,1000)
# Remove the top and right borders
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
