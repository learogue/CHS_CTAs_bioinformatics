#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------------------------
# Script : Merge Unique Mutated Protein FASTA
# Author : Léa ROGUE
# Date : 18/09/2025
# Description : This script scans all mutated protein FASTA files located in ../results/mut_fasta/, extracts their 
# sequences, and produces a merged FASTA file containing only unique protein sequences.
#     1. Read every FASTA file inside ../results/mut_fasta/
#     2. Parse headers and protein sequences
#     3. Deduplicate sequences (sequence = key, header = value)
#     4. Output a single FASTA file containing unique sequences only
#     5. Format sequences to 60 characters per line
# ----------------------------------------------------------------------------------------------------------------------

from collections import defaultdict
import glob

d_fasta = {}  # key = sequence, value = header

# Read fasta files and store unique sequences
for file in glob.glob('../results/mut_fasta/*'):
    with open(file, 'r') as f:
        current_id = ''
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save the previous sequence
                if current_id and seq_lines:
                    seq = ''.join(seq_lines)
                    if seq not in d_fasta:  # unique
                        d_fasta[seq] = current_id
                current_id = line
                seq_lines = []
            else:
                seq_lines.append(line)
        # Save the last sequence in the file
        if current_id and seq_lines:
            seq = ''.join(seq_lines)
            if seq not in d_fasta:
                d_fasta[seq] = current_id

# Write merged unique sequences to output fasta
with open('../results/mut_fasta_merged_unique.fasta', 'w') as out_f:
    for seq, header in d_fasta.items():
        out_f.write(header + '\n')
        for i in range(0, len(seq), 60):
            out_f.write(seq[i:i+60] + '\n')
