# Transient YAP hyperactivation enables the generation of assembled human blastoids

This repository contains the code associated with the published paper titled: "Transient YAP hyperactivation enables the generation of assembled human blastoids"

## Data Analysis

The data analysis for this project is divided into two main parts:

### Shell Code
- **Purpose**: To generate per-cell gene expression counts and analyze cell receptor–ligand interactions.
- **Contents**: Shell scripts for running CellRanger (v8.0.1) and CellphoneDB (v5.0.0).

### R Code
- **Purpose**: Applied for cell quanlity control, integration analysis, and data visualization.
- **Contents**: R scripts, gene lists, and cell type information that used in all figures.

For more detailed information, please refer to the respective scripts and documentation within this repository.

*To avoid potential information leakage risks, all absolute paths have been removed before uploading the code. Readers are advised to configure the path information based on their local environment.

## Author Contributions
J.W., H.M.W. and H.W. conceptualized the study. J.W. and H.M.W. supervised this study. H.W. developed YAP-induced functional TE-like cells (iYAP-TE cells), and established the protocol for TE spheroids, human/chimeric blastoids and human ICM aggregates based on iYAP-TE cells. H.W. established epithelial cell-based attachment model and conducted related experiments investigating attachment mechanism. S.O. conducted all mouse embryo-related experiments shown in Figure 1F-J and S1H. L.L. performed the scRNA-seq analysis and built the GitHub codebase. W.X. and X.Q, helped with the library construction for scRNA-seq. H.W., J.W. and H.M.W. wrote the draft.
