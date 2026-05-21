#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------------------------
# Script : Clean NetMHC results and create a table
# Author : Léa ROGUE
# Date : 21-10-25
# Description : This script processes the raw output from NetMHC predictions for the CH2879 and JJ012 cell lines. It 
# extracts relevant information such as allele, peptide sequence, binding affinity, and rank, and organizes it into a 
# structured table format. The cleaned data is saved as TSV files for downstream analysis.
# ----------------------------------------------------------------------------------------------------------------------

import re
import pandas as pd
import glob

samples = ['JJ012', 'CH2879', 'all']

for sample in samples:
    data = []  # List to store parsed rows

    # Retrieve the prediction file for the current sample
    files = glob.glob(f'../data/all_predictions_{sample}.txt')

    for filepath in files:
        with open(filepath, 'r') as f:
            for lig in f:
                lig = lig.strip()  # Remove leading/trailing whitespace and newlines

                # Keep only data lines (starting with a digit)
                # Skips headers, blank lines, and NetMHC comment lines
                if re.match(r'^\d', lig):
                    lig = lig.split()

                    # Extract relevant columns: Allele (1), Peptide (2), iCore (9), and trailing columns (11+)
                    # Indices match the tabular output format of NetMHC
                    lig = [lig[1], lig[2], lig[9]] + lig[11:]

                    # Remove the "<=" token used by NetMHC to flag binders (e.g., <= SB, <= WB)
                    lig = [elt for elt in lig if elt != "<="]

                    # Pad with None if the row is incomplete (e.g., no BindLevel for non-binders)
                    while len(lig) < 7:
                        lig.append(None)

                    data.append(lig)

    # Column names matching the extracted fields
    columns = ['Allele', 'Peptide', 'iCore', '1-log50k(aff)', 'Affinity', 'Rank', 'BindLevel']

    df = pd.DataFrame(data, columns=columns)
    df = df.drop_duplicates()  # Remove any duplicate entries

    # Save the cleaned table as a TSV file
    output_path = f'../data/table_res_netmhc_all_{sample}.tsv'
    df.to_csv(output_path, sep='\t', index=False)