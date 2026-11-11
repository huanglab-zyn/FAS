# huanglab
R scripts for single-cell RNA sequencing data analysis and figure generation.
## Overview
This repository contains the R scripts used for the data analysis and visualization in our manuscript. The scripts are organized chronologically by figures. Specific panels (e.g., Figure 1A, Figure 2C) are explicitly annotated as comments within the corresponding R scripts.

## Data Availability
NCBI BioProject：PRJNA1522586
NCBI GEO： GSE135779、GSE214974、GSE133028

## Repository Structure
*   `Figure1.R` / `Figure2.R` / `Figure4.R` / `Figure6.R` / `Figure7.R`
*   `Supplementary Figure 6.R` / `Supplementary Figure 7.R`

## Dependencies
*   **Main Analysis Packages:** `Seurat` (v4.3.0), `DESeq2` (v1.50.2), `metacell` (v0.3.7), `tgconfig` (v0.1.2)
*   **Data Manipulation & Visualization:** `dplyr`, `tidyr`, `stringr`, `foreach`, `ggplot2`, `pheatmap`, `ggrepel`, `svglite`, `Matrix`
