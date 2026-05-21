#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------------------------
# Script : Generate mutated protein FASTA from VEP annotations
# Author : Léa ROGUE
# Date : 18/09/2025
# Description : This script parses VEP-annotated variants (VCF) and UniProt human protein FASTA sequences to generate 
# mutated protein sequences for a given sample.
#     1. Load UniProt Swiss-Prot human protein sequences (FASTA)
#     2. Extract gene names and UniProt IDs
#     3. Parse VEP annotations from <sample>.vep.vcf
#     4. Collect HGVS protein changes for each gene
#     5. Convert 3-letter amino acid codes to 1-letter notation
#     6. Apply amino acid substitutions, stops, and frameshifts
#     7. Write mutated protein sequences in FASTA format
#     8. Log mismatches, out-of-range positions, missing genes, and unparsed HGVS entries
# ----------------------------------------------------------------------------------------------------------------------

from collections import defaultdict
import re

sample_id = 'JJ012'  # Sample identifier, used for input/output filenames

# Read UniProt FASTA
d_uniprot = {}        # Dictionary to store human protein sequences and UniProt IDs
current_id = ''       # Track current gene while parsing
human_prot = False    # Flag to track if the current protein is human

with open('../data/uniprot_sprot.fasta') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):                     # Header line
            if 'HUMAN' in line:                      # Only consider human proteins
                uniprot_id = line.split('|')[1]      # Extract UniProt ID
                m = re.search(r'GN=([^ ]+)', line)   # Extract gene name
                if m:
                    gene_name = m.group(1)
                    # Initialize dictionary entry for this gene
                    d_uniprot[gene_name] = {'seq': '', 'uniprot_id': uniprot_id}
                    current_id = gene_name
                    human_prot = True
                else:
                    human_prot = False
            else:
                human_prot = False
        else:
            # Append sequence lines only if it's a human protein
            if human_prot and current_id:
                d_uniprot[current_id]['seq'] += line

# Read VCF file (variants)
d_vcf = defaultdict(list)                                           # Store HGVS protein changes per gene
with open(f'../results/{sample_id}/{sample_id}.vep.vcf') as vcf:
    for line in vcf:
        if not line.startswith('#'):                                # Skip VCF headers
            info_field = line.strip().split('\t')[7]                # INFO column
            hgvsp = info_field.split('|')[11].split(':')[-1]        # HGVS protein notation
            gene_name = info_field.split('|')[3]                    # Gene name
            if hgvsp != '':
                d_vcf[gene_name].append(hgvsp)

# Amino acid conversion dictionary (3-letter to 1-letter)
aa3to1 = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Glu':'E','Gln':'Q',
          'Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F',
          'Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V', 'Ter':'*'}

# Convert HGVS protein changes to 1-letter code and filter invalid entries
d_vcf_new = defaultdict()
for gene, hgvs_list in d_vcf.items():
    new_list = []
    for h in hgvs_list:
        if not h:
            continue

        # Replace URL-encoded characters or frameshift notation
        h = h.replace('%3D', '=')           # Synonymous mutation
        h = h.replace('fsTer', 'fs*')       # Frameshift stop

        # Ignore synonymous mutations
        if h.endswith('='):
            continue

        # Match substitution or stop codon
        m = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)', h)

        # Match frameshift mutation
        m_fs = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})(fs\*\d+)', h)

        if m:
            aa_from3, pos, aa_to3 = m.groups()
            aa_from1 = aa3to1.get(aa_from3, 'X')
            aa_to1   = aa3to1.get(aa_to3, 'X')
            new_h = f'{aa_from1}{pos}{aa_to1}'
            new_list.append(new_h)

        if m_fs:
            aa_from3, pos, aa_to3, fs_stop = m_fs.groups()
            aa_from1 = aa3to1.get(aa_from3, 'X')
            aa_to1   = aa3to1.get(aa_to3, 'X')
            new_h = f'{aa_from1}{pos}{aa_to1}{fs_stop}'
            new_list.append(new_h)

    if new_list != []:
        d_vcf_new[gene] = new_list

# Output mutated FASTA sequences and log warnings
with open(f'../results/{sample_id}/{sample_id}_mutated_gene_mut.fasta', 'w') as out_f, \
     open(f'../results/{sample_id}/{sample_id}_warnings_gene_mut.log', 'w') as log_f:

    for gene, l_mutations in d_vcf_new.items():
        for mutation in l_mutations:
            if gene in d_uniprot:
                seq = d_uniprot[gene]['seq']
                uniprot_id = d_uniprot[gene]['uniprot_id']

                # Match substitution or stop codon
                m = re.match(r'([A-Z])(\d+)([A-Z*])', mutation)
                # Match frameshift
                m_fs = re.match(r'([A-Z])(\d+)([A-Z*])fs\*(\d+)', mutation)

                if m:
                    ref_aa, pos, alt_aa = m.groups()
                    pos = int(pos) - 1
                    if pos < 0 or pos >= len(seq):
                        log_f.write(f'[OUT OF RANGE] {gene} {mutation} (len={len(seq)})\n')
                        continue
                    if seq[pos] != ref_aa:
                        log_f.write(f'[MISMATCH] {gene} {mutation} expected {ref_aa}, found {seq[pos]}\n')
                        continue

                    # Apply mutation: stop codon truncates sequence
                    if alt_aa == '*':  
                        mutated_seq = seq[:pos]
                    else:
                        mutated_seq = seq[:pos] + alt_aa + seq[pos+1:]

                elif m_fs:
                    ref_aa, pos, alt_aa, stop_dist = m_fs.groups()
                    pos = int(pos) - 1
                    stop_dist = int(stop_dist)
                    if pos < 0 or pos >= len(seq):
                        log_f.write(f'[OUT OF RANGE] {gene} {mutation} (len={len(seq)})\n')
                        continue
                    # Frameshift: replace amino acids up to predicted stop
                    mutated_seq = seq[:pos] + alt_aa + seq[pos+1:pos+1+stop_dist]

                else:
                    log_f.write(f'[UNPARSED] {gene} {mutation}\n')
                    continue

                # Write FASTA output: use UniProt ID as header
                out_f.write(f'>{uniprot_id}|{mutation}\n')
                for i in range(0, len(mutated_seq), 60):
                    out_f.write(mutated_seq[i:i+60] + '\n')
            else:
                log_f.write(f'[MISSING] {gene} not in UniProt\n')
