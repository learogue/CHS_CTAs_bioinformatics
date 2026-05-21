#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------------------------------------------------
# R script : Preprocessing and processing step to generate expression matrix and immune cells means matrix
# Author  : Léa ROGUE
# Date    : 23-01-2025
# Description : This script use oligo to process microarrays data. 
#   1. Read all the files and create the object with metadata from the .sdrf file to have all the informations. 
#   2. Quality analysis 
#   3. Expression data are normalized with the RMA algorithm (from oligo)
#   4. Probes annotation
#   5. Clean and organize the dataset by consolidating duplicate probes and calculating mean expression values.
#   6. Compute Z-scores to standardize gene expression levels.
# ------------------------------------------------------------------------------------------------------------------------

# If an error occur for the rma fonction, these commands line solve it
#BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE)
#BiocManager::install("oligo", configure.args="--disable-threading", force = TRUE)

# Function to calculate Z scores on rows
calculate_z_scores <- function(df_input, col) {
  # Exclude col PROBEID SYMBOL and CTA
  data_values <- df_input[, -c(col)]
  
  # Calculate Z-scores
  z_scores_row <- t(scale(t(data_values)))
  
  # Add columns
  df_z_scores <- cbind(df_input[, c(col)], z_scores_row)
  
  # Return df
  return(df_z_scores)
}

# Load packages
library(oligo)
library(arrayQualityMetrics)
library(hugene20sttranscriptcluster.db)
library(dplyr)

# Directories paths
data_dir <- "../data/E-MTAB-7264_full"
output_dir <- "../results"
dir.create(output_dir)

# SRDF to store the metadata
SDRF <- read.delim(paste0(data_dir, "/E-MTAB-7264.sdrf.txt"))

# Give raw names the sample name 
rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)

# Read the raw data and associate metadata
raw_data <- read.celfiles(filenames = file.path(data_dir, SDRF$Array.Data.File), phenoData = SDRF)

# Normalization and delete background noise
norm_data <- rma(raw_data, target = "core")

# Array quality control
arrayQualityMetrics(expressionset = raw_data, outdir = "../results/qc/qc_raw", force = TRUE, do.logtransform = TRUE)
arrayQualityMetrics(expressionset = norm_data, outdir = "../results/qc/qc_norm", force = TRUE, do.logtransform = TRUE)

# Save the metadata
metadata_col <- c("Characteristics.exp.subtypes.",
                  "Characteristics.mom.subtypes.",
                  "Characteristics.idh1.aamut.",
                  "Characteristics.idh1.freq.",
                  "Characteristics.idh2.aamut.",
                  "Characteristics.idh2.freq.",
                  "Characteristics.col2a1.",
                  "Characteristics.tp53.",
                  "Characteristics.cdkn2a.copynumber.",
                  "Characteristics.os.delay.",
                  "Unit.time.unit.",
                  "Characteristics.event.death.",
                  "Characteristics.tumor.grading.",
                  "Factor.Value.tumor.grading.")
df <- pData(norm_data)[,metadata_col]
write.table(df, file = "../results/metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Take only the coding genes
mapped_probes <- mappedkeys(hugene20sttranscriptclusterGENENAME)

# Associate genes with probes
anno <- AnnotationDbi::select(hugene20sttranscriptcluster.db, keys = mapped_probes, keytype = "PROBEID", columns = c("SYMBOL", "GENENAME"))

# Take rows that have a gene
anno <- subset(anno, !is.na(SYMBOL))

# Change names of column with patient IDs
colnames(norm_data) <- norm_data$Source.Name
rownames(pData(norm_data)) <- norm_data$Source.Name

# Expression data and rename the rows with probe id
expr_norm_data_df <- as.data.frame(exprs(norm_data))
expr_norm_data_df$PROBEID <- rownames(expr_norm_data_df)

# Merge expression data with SYMBOL
merged_df <- merge(anno, expr_norm_data_df, by = "PROBEID", all.x = TRUE)
merged_df <- merged_df[, c(-1,-3)]

# For genes that have multiple probe, calculate the mean
merged_df_avg <- merged_df %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# Read CTA file
merged_df_avg$CTA <-  NA
df_CTA <- read.table("../data/CTA_list_clean.txt", header = FALSE)

# Rename col
colnames(df_CTA) <- c("SYMBOL")

# Df with CTA and whole genes by adding CTA in col CTA for CTA genes
df_CTA_whole <- merged_df_avg %>%
  left_join(df_CTA, by = "SYMBOL") %>%
  mutate(CTA = ifelse(SYMBOL %in% df_CTA$SYMBOL, "CTA", NA))

# Reorganize columns
df_CTA_whole <- df_CTA_whole %>%
  dplyr::select(SYMBOL, CTA, everything())

# Read immune cells genes
df_immune_sign <- read.table("../data/immune_cells_genes.tsv", header = FALSE, sep = "\t")
colnames(df_immune_sign) <- c("Signature", "Gene")

# Merge df_CTA_whole and df_immune_sign to associate signatures
df_CTA_immune_sign_whole <- merge(df_CTA_whole, df_immune_sign, by.x = "SYMBOL", by.y = "Gene", all.x = TRUE)

# Reorganize the columns if needed
df_CTA_immune_sign_whole <- df_CTA_immune_sign_whole %>%
  select(SYMBOL, CTA, Signature, everything())
write.table(df_CTA_immune_sign_whole, "../results/matrix_complete_intensities.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Combine multiple signatures into one for each gene by concatenating with a comma
df_CTA_immune_sign_whole_clean <- df_CTA_immune_sign_whole %>%
  group_by(SYMBOL) %>%
  summarise(
    Signature = paste(unique(Signature), collapse = ", "),                                              
    across(everything(), ~first(.)),  
    .groups = "drop")
write.table(df_CTA_immune_sign_whole_clean, file = "../results/whole_gene_int_CTA_sign_imm_clean.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
