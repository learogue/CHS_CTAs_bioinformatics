#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------------------------
# Script : Compute TPM and FPKM
# Author : Léa ROGUE
# Date : 27-06-2025
# Description : This script converts raw gene-level read counts obtained from featureCounts into normalized expression
# values (FPKM and TPM) for CH2879 and JJ012 cell lines. It integrates gene lengths extracted from the featureCounts 
# output and gene annotations retrieved from the GTF file (GRCh38.114).
# ----------------------------------------------------------------------------------------------------------------------

from collections import defaultdict
import re

# Transform gene ID to gene name
id_to_name = {}
with open('../data/Homo_sapiens.GRCh38.114.gtf', 'r') as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        infos = parts[8]
        match_id = re.search(r'gene_id "([^"]+)"', infos)
        match_name = re.search(r'gene_name "([^"]+)"', infos)
        if match_id and match_name:
            gene_id = match_id.group(1)
            gene_name = match_name.group(1)
            id_to_name[gene_id] = gene_name

# Samples to process
samples = ["CH2879", "JJ"]

# Process
for sample in samples:
    d_counts = defaultdict()    
    d_length = defaultdict() 

    # Read counts and gene lengths
    with open(f'../results/counts/{sample}_fragments_counts.tab', 'r') as f:
        next(f) # skip header line1
        next(f) # skip header line2
        for lig in f:
            lig = lig.strip().split('\t')
            d_length[lig[0]] = int(lig[5])
            d_counts[lig[0]] = int(lig[6])

    # Commpute total
    total_counts = sum(d_counts.values())

    # Compute FPKM
    d_fpkm = defaultdict()
    for gene in d_counts:
        count = d_counts[gene]
        length = d_length[gene]
        d_fpkm[gene] = (count * 1e9) / (length * total_counts)
    sum_fpkm = sum(d_fpkm.values())

    # Compute TPM
    d_tpm = defaultdict()
    for gene in d_counts:
        count = d_counts[gene]
        length = d_length[gene]
        d_tpm[gene] = ((count * 1e3) / length) * (1 / sum_fpkm) * 1e6

    # Write file
    with open(f'../results/{sample}_fpkm_tpm.tsv', 'w') as out:
        out.write("Gene_ID\tGene_Name\tLength\tCount\tFPKM\tTPM\n")
        for gene in d_counts:
            count = d_counts[gene]
            length = d_length[gene]
            fpkm = d_fpkm[gene]
            tpm = d_tpm[gene]
            gene_name = id_to_name.get(gene, "NA")
            out.write(f"{gene}\t{gene_name}\t{length}\t{count}\t{fpkm}\t{tpm}\n")