# Single-cell RNA-seq Analysis of CCl₄-induced Liver Fibrosis

This repository contains R scripts for analyzing single-cell RNA sequencing data from CCl₄-induced liver fibrosis mouse models, as described in [Manuscript Title].

## Overview

- **Cell type annotation** and clustering (Figure 1)
- **Hepatocyte zonation** analysis and stress response (Figure 2)  
- **HSC activation** status and fibrosis mechanisms (Figure 3)
- **Cell-cell communication** analysis using CellChat (Figure 4)
- **Therapeutic target** prioritization (Figure 5)

## System Requirements

- R ≥ 4.2.0
- RStudio (recommended)

## Dependencies

### Core Packages
- Seurat (v4/v5)
- tidyverse
- ggplot2
- dplyr

### Analysis-specific Packages
| Analysis | Package | Version |
|:---|:---|:---|
| Cell communication | CellChat | ≥ 1.6.0 |
| Pathway enrichment | clusterProfiler | ≥ 4.0 |
| Visualization | ComplexHeatmap, circlize | latest |

## Installation

```r
# Install CRAN packages
install.packages(c("Seurat", "tidyverse", "ggplot2", "dplyr"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "ComplexHeatmap"))

# Install CellChat from GitHub
devtools::install_github("jinworks/CellChat")

Repository Structure
├── data/                  # Input data (not included, see Data Access)
├── scripts/
│   ├── figure1_cell_annotation.R
│   ├── figure2_hepatocyte_zonation.R
│   ├── figure3_hsc_activation.R
│   ├── figure4_cellchat_communication.R
│   └── figure5_therapeutic_targets.R
├── results/               # Output figures and tables
├── README.md
└── LICENSE

