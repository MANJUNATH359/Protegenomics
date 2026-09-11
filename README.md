# Proteogenomics: Six-Frame Translation and Peptide Filtering

Python scripts developed as part of the proteogenomic analysis used to refine
the *Entamoeba histolytica* HM-1:IMSS reference genome.

The workflow was used in the following publication:

Malshetty, M.B., Chowdhury, S., Pawar, S. et al. Application of Proteogenomic Approaches for Refinement of the Entamoeba histolytica Reference Genome Using High Resolution Mass-Spectrometry Data. Acta Parasit. 71, 166 (2026). https://doi.org/10.1007/s11686-026-01352-8 

## Overview

This repository contains Python scripts used in the computational workflow
for generating a custom six-frame translated protein database and filtering
peptides against the annotated *E. histolytica* reference proteome.

The workflow supported the identification of genome search-specific peptides
(GSSPs) for downstream proteogenomic analysis and genome annotation refinement.

## Workflow

**Reference genome → Six-frame translation → Custom ORF database → LC-MS/MS
peptide identification → Annotated protein filtering → GSSP analysis**

## Scripts

### 1. Six-frame ORF translation

`orf_finder_enthisto_unbiased_final.py`

Generates a custom six-frame translated protein database from an
*E. histolytica* genome FASTA file.

The script:

- Translates the genome in all six reading frames (+1, +2, +3, −1, −2, −3)
- Identifies open reading frames
- Applies a minimum ORF length of 10 amino acids
- Records genomic coordinates and reading-frame information
- Generates a FASTA file containing predicted ORF protein sequences
- Generates an ORF log containing genomic coordinates and ORF information

### 2. Annotated protein peptide filtering

`filter_annotated_peptides.py`

Filters peptide sequences against the annotated *E. histolytica* HM-1:IMSS
protein database.

The script:

- Reads peptide sequences from an Excel file
- Removes modification annotations from peptide sequences
- Searches peptides against annotated protein sequences
- Separates peptides that match the reference proteome from those that do not
- Generates matched and unmatched peptide tables

The unmatched peptide set can subsequently be used for downstream
proteogenomic/GSSP analysis.

## Data Sources

The study used:

- *Entamoeba histolytica* HM-1:IMSS reference genome from AmoebaDB release 68
- Annotated protein sequences from AmoebaDB release 68
- Proteomic datasets from PRIDE:
  - PXD042282
  - PXD051913

The raw proteomic datasets are publicly available through the PRIDE repository.

## Reproducibility

The six-frame translation workflow generated approximately 1.16 million
predicted ORFs using a minimum translated ORF length of 10 amino acids.

In the published analysis, 38,562 unique peptides were identified from the
proteomic datasets. Of these, 480 peptides did not match the annotated
*E. histolytica* protein database and were investigated as genome
search-specific peptides (GSSPs).

These GSSPs contributed to the identification of 41 novel protein-coding
genes and refinement of existing gene models.
