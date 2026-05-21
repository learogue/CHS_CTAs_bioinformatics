#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------------------------------------------------
# R script : Fasta files creation from immunopeptidomic files
# Author : Léa ROGUE
# Date : 30-01-2026
# Description : This script loads immunopeptidomic tables and create fasta files to run netMHC.
# ------------------------------------------------------------------------------------------------------------------------

library(dplyr)
library(tidyr)

# Read data
df_all <- read.csv("../data/immunopep_all.tsv", sep = "\t", header = T)

# Delete contaminant peptides
df_all_clean <- df_all[df_all$Contaminant != T,]
write.table(df_all_clean, "../data/immunopep_all_clean.tsv", sep = "\t", quote = F, row.names = F)

# CH2879 abundances
df_ch_mhc1 <- df_all_clean[, c("Annotated.Sequence", "Master.Protein.Accessions", "Sequence.Length", "Abundances..Grouped...CH2879..NoINF")]
df_ch_mhc1$Seq <- sub(".*\\.(.*?)\\..*", "\\1", df_ch_mhc1$Annotated.Sequence)

# JJ012 abundances
df_jj_mhc1 <- df_all_clean[, c("Annotated.Sequence", "Master.Protein.Accessions", "Sequence.Length", "Abundances..Grouped...JJ012..NoINF")]
df_jj_mhc1$Seq <- sub(".*\\.(.*?)\\..*", "\\1", df_jj_mhc1$Annotated.Sequence)

# Clean tables
pep_ch_mhc1 <- df_ch_mhc1 %>%
  filter(!is.na(Abundances..Grouped...CH2879..NoINF) & Abundances..Grouped...CH2879..NoINF > 0) %>%
  group_by(Seq) %>%
  slice_max(Abundances..Grouped...CH2879..NoINF, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Seq, Prot_accessions = Master.Protein.Accessions, seq_len = Sequence.Length, abundances_CH2879_noINF = Abundances..Grouped...CH2879..NoINF)

pep_jj_mhc1 <- df_jj_mhc1 %>%
  filter(!is.na(Abundances..Grouped...JJ012..NoINF) & Abundances..Grouped...JJ012..NoINF > 0) %>%
  group_by(Seq) %>%
  slice_max(Abundances..Grouped...JJ012..NoINF, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Seq, Prot_accessions = Master.Protein.Accessions, seq_len = Sequence.Length, abundances_JJ012_noINF = Abundances..Grouped...JJ012..NoINF)

# Save tables
write.table(pep_ch_mhc1, "../data/immunopep_CH2879.tsv", sep = "\t", quote = F, row.names = F)
write.table(pep_jj_mhc1, "../data/immunopep_JJ012.tsv", sep = "\t", quote = F, row.names = F)

# Write fasta to run netMHC
fasta_ch <- paste0(">", seq_len(nrow(pep_ch_mhc1)), "\n", pep_ch_mhc1$Seq)
writeLines(fasta_ch, "../results/peptides_CH2879.fa")

fasta_jj <- paste0(">", seq_len(nrow(pep_jj_mhc1)), "\n", pep_jj_mhc1$Seq)
writeLines(fasta_jj, "../results/peptides_JJ012.fa")

fasta_all <- paste0(">", seq_len(nrow(df_all_clean)), "\n", df_all_clean$Seq)
writeLines(fasta_all, "../results/peptides_all.fa")
