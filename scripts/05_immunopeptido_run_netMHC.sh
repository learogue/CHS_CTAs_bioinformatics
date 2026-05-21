#!/bin/bash
# -----------------------------------------------------------------------------------------------------------------
# bash script : Run netMHC for peptide binding prediction
# Auteur : Léa ROGUE
# Date : 22-09-25
# Description : This script runs the netMHC tool to predict peptide binding affinities for a comprehensive list of 
# HLA alleles. It takes peptide sequences from FASTA files and generates binding predictions for peptides of lengths 
# 8 to 11 amino acids. The results are saved in text files for further analysis.
# -----------------------------------------------------------------------------------------------------------------

# Take all alleles from the list
ALLELES=( HLA-A0101 HLA-A0201 HLA-A0202 HLA-A0203 HLA-A0205 HLA-A0206 HLA-A0207 HLA-A0211 HLA-A0212 HLA-A0216 HLA-A0217 HLA-A0219 HLA-A0250 HLA-A0301 HLA-A0302 HLA-A0319 HLA-A1101 HLA-A2301 HLA-A2402 HLA-A2403 HLA-A2501 HLA-A2601 HLA-A2602 HLA-A2603 HLA-A2902 HLA-A3001 HLA-A3002 HLA-A3101 HLA-A3201 HLA-A3207 HLA-A3215 HLA-A3301 HLA-A6601 HLA-A6801 HLA-A6802 HLA-A6823 HLA-A6901 HLA-A8001 HLA-B0702 HLA-B0801 HLA-B0802 HLA-B0803 HLA-B1401 HLA-B1402 HLA-B1501 HLA-B1502 HLA-B1503 HLA-B1509 HLA-B1517 HLA-B1801 HLA-B2705 HLA-B2720 HLA-B3501 HLA-B3503 HLA-B3701 HLA-B3801 HLA-B3901 HLA-B4001 HLA-B4002 HLA-B4013 HLA-B4201 HLA-B4402 HLA-B4403 HLA-B4501 HLA-B4506 HLA-B4601 HLA-B4801 HLA-B5101 HLA-B5301 HLA-B5401 HLA-B5701 HLA-B5703 HLA-B5801 HLA-B5802 HLA-B7301 HLA-B8101 HLA-B8301 HLA-C0303 HLA-C0401 HLA-C0501 HLA-C0602 HLA-C0701 HLA-C0702 HLA-C0802 HLA-C1203 HLA-C1402 HLA-C1502 HLA-E0101 HLA-E0103 )

# Convert the array of alleles into a comma-separated string
ALLELES_CSV=$(IFS=, ; echo "${ALLELES[*]}")

# Run NetMHC
netMHC -a "$ALLELES_CSV" -f ../results/peptides_CH2879.fa -l 8,9,10,11 > ../data/all_predictions_CH2879.txt
netMHC -a "$ALLELES_CSV" -f ../results/peptides_JJ012.fa -l 8,9,10,11 > ../data/all_predictions_JJ012.txt
netMHC -a "$ALLELES_CSV" -f ../results/peptides_all.fa -l 8,9,10,11 > ../data/all_predictions_all.txt

# Select strong binders (SB) based on the %Rank 
cat ../data/all_predictions_CH2879.txt | grep "SB" > ../data/all_sb_CH2879.txt
cat ../data/all_predictions_JJ012.txt | grep "SB" > ../data/all_sb_JJ012.txt
cat ../data/all_predictions_all.txt | grep "SB" > ../data/all_sb_all.txt
