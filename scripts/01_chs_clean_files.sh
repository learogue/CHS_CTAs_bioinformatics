#!/bin/bash
# -----------------------------------------------------------------------------------------------------------------
# bash script : Cleaning files to not use all the raw files
# Author : Léa ROGUE
# Date : 27-01-2025
# Description : This script is useful for cleaning and preprocessing gene expression data. This script filters .CEL
# files by checking if their filenames contain an ID from the list stored in list_files_used_publi.txt because some
# files came from same patients. The 1st step is to read the list of valid IDs from list_files_used_publi.txt. The
# 2nd step is to iterate through all .CEL files in the directory. The 3rd step is to check if the filename contains
# a valid ID, if the filename matches an ID, keep the file, otherwise, delete the file. Finally, the metadata file
# (E-MTAB-7264.sdrf.txt) is updated by keeping only the relevant lines.
# -----------------------------------------------------------------------------------------------------------------

# Path for ID to keep
list_file="../data/list_files_used_publi.txt"

# Read file name of .CEL files
for file in ../data/E-MTAB-7264_full/*.CEL; do
    # Extract
    basename_without_extension=$(basename "$file" .CEL)

    # Initiate a var to know if the file have to be delete
    delete_file=true

    # Verify if the ID is in the file name
    while read -r id; do
        if [[ "$basename_without_extension" == *"$id"* ]]; then
            delete_file=false
            break
        fi
    done < "$list_file"

    # If not, file deleted
    if $delete_file; then
        echo "$file deleted"
        rm "$file"
    fi
done

# Modify the metadata and delete lines
head -n 1 ../data/E-MTAB-7264_full/E-MTAB-7264.sdrf.txt > tmp.txt
awk 'NR==FNR {ids[$1]; next} $1 in ids' ../data/list_files_used_publi.txt ../data/E-MTAB-7264_full/E-MTAB-7264.sdrf.txt >> tmp.txt
cp tmp.txt ../data/E-MTAB-7264_full/E-MTAB-7264.sdrf.txt
rm tmp.txt
