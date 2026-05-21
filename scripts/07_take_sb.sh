#!/bin/bash
# -----------------------------------------------------------------------------------------------------------------
# bash script : Extract strong binder peptides from NetMHC results
# Author : Léa ROGUE
# Date : 27-06-2025
# Description : This script extracts peptides classified as strong binders (SB) from the NetMHC prediction results 
# for the JJ012 cell line. It filters the results based on the 'SB' classification and saves the relevant peptides 
# and their associated data into a new TSV file for further analysis.
# -----------------------------------------------------------------------------------------------------------------

# Take first line of the table to get the header, then filter lines with 'SB' and save to new file
head -1 ../data/table_res_netmhc_all_JJ012.tsv > ../data/table_res_netmhc_all_sb_JJ012.tsv   
cat ../data/table_res_netmhc_all_JJ012.tsv | grep 'SB' >> ../data/table_res_netmhc_all_sb_JJ012.tsv

head -1 ../data/table_res_netmhc_all_CH2879.tsv > ../data/table_res_netmhc_all_sb_CH2879.tsv
cat ../data/table_res_netmhc_all_CH2879.tsv | grep 'SB' >> ../data/table_res_netmhc_all_sb_CH2879.tsv