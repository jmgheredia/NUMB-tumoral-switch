# Data instructions

This repository contains code to reproduce analyses and figures for the manuscript:

"Tumoral Switch in NUMB splicing changes essential transcription pathways and induces malignant properties of tumour cells"

Large datasets (e.g. TCGA/CCLE) and third-party supplementary files are not redistributed in this GitHub repository.

------------------------------------------------------------------------

## Required third-party file (external input)

Some scripts require the following precomputed IsoformSwitchAnalyzeR object:

-   File name: `File_S1_switchAnalyzeRlists.Rdata`
-   Place the file here (local path): `data/external/File_S1_switchAnalyzeRlists.Rdata`

Source (persistent reference): - Figshare DOI: <https://doi.org/10.6084/m9.figshare.4924724>

Note: - Direct download URLs may change over time; the DOI above is the stable record to locate the file.

------------------------------------------------------------------------

## TCGA / CCLE

TCGA and CCLE datasets were used in this study but are not redistributed in this repository.

Users must obtain the corresponding TCGA and CCLE files from their official sources:

- TCGA: NIH Genomic Data Commons (GDC)
- CCLE: Firebrowse 

For reproducibility, the exact file names, releases, and preprocessing steps used in this study are documented in the corresponding analysis scripts.


All downstream scripts assume that these files are placed in `data/external/` following the paths specified in each script header.

------------------------------------------------------------------------

## PDMR (Patient-Derived Model Repository) – RSEM archives

Isoform-level expression matrices used for Figure 4 are derived from RSEM output files obtained from the Patient-Derived Model Repository (PDMR).

Due to redistribution restrictions, these RSEM .results files are not included in this repository and must be downloaded manually by the user.
