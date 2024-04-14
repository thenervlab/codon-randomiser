import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger

logger.info('Import OK')

input_folder = ''
output_folder = ''

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Optional: Set seed


# Read in codon table


# Read in sequence file -  from txt or fasta?


# Verify sequence is in groups of 3 codons


# Locate start codon (ideally at the start, if not find)



# Locate stop codon (ideally at end, if missing add warning)



# For each codon in the sequence:
## Identify AA coded
## Identify options for other codons
## Randomly select a codon
## Apppend to new sequence


# Calculate % identity with original sequence
## At the codon level
## At the nucleotide level


# Optional: Visualise comparison 
# -> ball & stick plot where each codon is a position,
# and each nucleotide within that codon that is matched
# scores +1 - ideally matched codons would not appear
# too many in a row


# Save (download) updated sequence and 
# randomisation report


