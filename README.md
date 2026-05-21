# CTAs bioinformatics analysis in chondrosarcoma

This repository contains scripts used for the bioinformatic analyses presented in the article:

> **Identification of immunopeptides (pHLA) as candidate therapeutic targets in chondrosarcoma**  
> Léa Rogue, Jean-Marc Monneuse, Céline Béchon, Lola Cepero, Caroline Peyrode, Sandrine Viala, Maud Privat, Yannick Bidet, Adrien Saliou, Paul-Olivier Rouzaire, Elisabeth Miot-Noirault, Florent Cachin and Aurélien Pommier
> *Manuscript in preparation*

The objective of this work is to characterize Cancer-Testis Antigens (CTAs) in chondrosarcoma, integrating transcriptomic and immunopeptidomic data to identify potential peptides presented on the surface of tumor cells via the HLA-I complex.

---

## Data availability

The datasets analyzed in this study are publicly available:
- **Microarray data (chondrosarcoma)**: ArrayExpress accession [E-MTAB-7264](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7264)
- **Protein expression (normal tissues)**: Human Protein Atlas - [https://www.proteinatlas.org/humanproteome/tissue/data#consensus_tissues_rna](https://www.proteinatlas.org/humanproteome/tissue/data#consensus_tissues_rna)
- **Immunopeptidomic data**: Immune Epitope Database - [https://www.iedb.org/database_export_v3.php](https://www.iedb.org/database_export_v3.php)

RNA-seq data from cell lines are available from the corresponding author upon reasonable request.

---

## Usage and requirements

This project uses multiple languages including R, Python, and Bash.

All required packages are listed in the `env/CHS_CTAs_bioinformatics.yml` file and can be installed using Conda.

We also used **NetMHC 4.0** (https://services.healthtech.dtu.dk/services/NetMHC-4.0) for HLA-I binding predictions (authorization is required to use this tool).

The full analysis pipeline can be reproduced by running the scripts in numerical order.

---

## Installation

```bash
git clone https://github.com/learogue/CHS_CTAs_bioinformatics.git
cd CHS_CTAs_bioinformatics
conda env create -f env/CHS_CTAs_bioinformatics.yml
conda activate CHS_CTAs_bioinformatics
Rscript scripts/install_packages.R
```

You can execute the scripts directly or use RStudio to run the RMarkdown files.

---

## Contents

```
CHS_CTAs_bioinformatics
├── data
│   ├── immunopep_all.tsv                       # Immunopeptidomic data (processed data for JJ012 and CH2879)
│   ├── immunopeptido_proteo.tsv                # Immunopeptidomic data (protein data processed data for JJ012 and CH2879)
│   ├── list_files_used_publi.txt               # List of microarray files used (duplicates excluded)
│   └── peptide_table_immg_pred.tsv             # Peptide table with HLA-I binding predictions
├── env
│   └── CHS_CTAs_bioinformatics.yml             # Conda environment file
├── scripts
│   ├── 00_chs_download_data.R                  # Download microarray data (CHS)
│   ├── 01_chs_clean_files.sh                   # Remove duplicate files
│   ├── 02_chs_process.R                        # Preprocessing and normalization (RMA)
│   ├── 03_chs_fig1.Rmd                         # Figure 1 - transcriptomic analysis
│   ├── 04_immunopeptido_fasta.R                # Generate FASTA files for netMHC
│   ├── 05_immunopeptido_run_netMHC.sh          # Run NetMHC for HLA-I binding predictions
│   ├── 06_clean_res_netmhc.py                  # Parse and clean NetMHC results
│   ├── 07_take_sb.sh                           # Filter strong binders
│   ├── 08_rnaseq_processing.sh                 # RNA-seq alignment pipeline (STAR + featureCounts)
│   ├── 09_rnasq_count_tpm.py                   # Compute TPM from counts
│   ├── 10_immunopeptido_fig2.Rmd               # Figure 2 - immunopeptidomic analysis
│   ├── 11_immunopeptido_CTAs_fig3.Rmd          # Figure 3 - CTA immunopeptidomic analysis
│   ├── 12_fig4.Rmd                             # Figure 4 - CTAs associated with OS and affinity
│   ├── 13_fig5.Rmd                             # Figure 5 - Immunophenotype
│   ├── 14_fig6.Rmd                             # Figure 6 - Immunogenicity
│   ├── 15_fig7.Rmd                             # Figure 7 - INFg
│   ├── apdx_0_pipeline_rnaseq_gatk.sh          # Appendix - variant calling pipeline (GATK)
│   ├── apdx_1_vcf_to_fasta.py                  # Appendix - VCF to FASTA conversion
│   ├── apdx_2_merge_fasta.py                   # Appendix - merge FASTA files
│   ├── install_packages.R                      # Install R packages not available via conda
│   ├── notebooks                               # PDF outputs from RMarkdown scripts
│   │   ├── 03_chs_fig1.pdf
│   │   ├── 10_immunopeptido_fig2.pdf
│   │   ├── 11_immunopeptido_CTAs_fig3.pdf
│   │   ├── 12_fig4.pdf
│   │   ├── 13_fig5.pdf
│   │   ├── 14_fig6.pdf
│   │   └── 15_fig7.pdf
│   └── tex_files
│       └── list_of_figures.tex                 # LaTeX file to generate list of figures
└── README.md
```

---

## Code authorship

All bioinformatic analyses, data processing pipelines, and figure generation scripts in this repository were designed and implemented by Léa ROGUE.

This repository accompanies a collaborative research project. Contributions from co-authors relate to study design, data generation, and biological interpretation as part of the associated publication.

---

## Contact

For questions regarding the code, please contact: lea.rogue@uca.fr

---

## License

This repository is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
