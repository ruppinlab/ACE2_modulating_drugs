# ACE2_modulating_drugs
Code for reproducing the results of our study titled "Systematic cell line-based identification of drugs modifying ACE2
expression" by Sanju Sinha*, Kuoyuan Cheng*, Alejandro A. Schäffer, Kenneth Aldape, Eyal Schiff, and Eytan Ruppin†.


### Additional data needed

The Level 3 LINCS l1000 data should be downloaded from the GEO database (GEO dataset identifiers given below) and put into the "data" folder with the following file names:
* Phase 1 data (GSE92742): data/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx
* Phase 2 data (GSE70138): data/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx

### Software dependencies

R was used for all analysis. Additional R packages required for analysis include:
* data.table
* cmapR
* limma
* affy
* GEOquery
* stringr
* fgsea

Additional R packages required for plotting the figures include:
* ggplot2
* ggrepel
* readxl
* ggpubr
* gridExtra
* grid
* lattice

### Summary of folders and files

The "data" folder contains processed data used for various analysis; The "analysis" folder contains various R scripts for the analyses in the study; The "plot" folder contains R scripts for generating the figures in the manuscript.

Specifically, key scripts in the "analysis" folder include:
* run.cmap.de.R: for running differential expression (DE) analysis on the CMAP data
* gsea.R: for running drug class enrichment analysis on the DE results
* geo.de.R: an example of DE analysis on a dataset from GEO

The other scripts in this folder are either files containing functions sourced by other scripts, or scripts used for selecting the relevant samples/datasets or preparing the samples for analysis.
