#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------------------------------------------------------------------
# R script : Downloading R packages
# Author  : Léa ROGUE
# Date    : 02-04-2026
# Description : To execute after creating conda env CHS_CTAs_bioinformatics
# ----------------------------------------------------------------------------------------------------------------------------------------

# CRAN packages
install.packages(c(
  "forestplot",   # 3.1.7
  "survminer",    # 0.5.2
  "ggsurvfit",    # 1.2.0
  "ggseqlogo",    # 0.2.2
  "VennDiagram",  # 1.8.2
  "writexl",      # 1.5.4
  "ggpubr",       # 0.6.3
  "ggsignif",     # 0.6.4
  "plotly",       # 4.12.0
  "WGCNA"         # 1.74
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("UniProt.ws")  # 2.46.1

# Create results directory
dir.create("results", showWarnings = FALSE)