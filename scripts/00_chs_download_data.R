#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------------------------------------------------------------------
# R script : Downloading CHS  microarray data E-MTAB-7264 from https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7264
# Author  : Léa ROGUE
# Date    : 21-01-2025
# Description : This R script use the package ArrayExpress which ask the data from the EBI database and download the raw data (.CEL). 
# These files are AffyMetrix microarray and the values are intensities. There is also metadata files to have the informations about the 
# experiments and the data.
# ----------------------------------------------------------------------------------------------------------------------------------------

# Loading packages
library(ArrayExpress)

# Data directory path
dir_files <- "../data/E-MTAB-7264_full"

# Create the dir
dir.create(dir_files)

# Download data from EBI
getAE("E-MTAB-7264", path = dir_files, type = "full", extract = FALSE)
