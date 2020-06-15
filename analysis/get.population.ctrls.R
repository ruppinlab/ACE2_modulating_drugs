##### This is the script for extracting the population control expression values from the CMAP data, to be used as control for the drug-treated vs control differential expression analysis
# the output of this script is already saved in ../data/population.ctrls.RDS

library(data.table)
library(cmapR)

# load info on drugs of interests: this loaded file contains information on clinically approved drugs collected from DrugBank
load("../data/drugs.of.interest.RData")
drugs <- d.app$pert_iname # all approved drugs

# load CMAP instance info (i.e. info on all available samples): the data in this loaded file was downloaded from CMAP in plain text format and saved as an RData file
load("../data/instinfo.RData")

# load CMAP gene info (i.e. info on all available samples): the data in this loaded file was downloaded from CMAP in plain text format and saved as an RData file
load("../data/geneinfo.RData")
# use only the best inferred genes (BING, this includes the landmark genes); ACE2 is one of the BING genes
gns <- geneinfo[pr_is_bing==1, .(id=as.character(pr_gene_id), symbol=pr_gene_symbol)]

# filter for the drugs, cell types, and treatment conditions of our interest; this represent all the 989 approved drugs available in CMAP
instinfo <- instinfo[pert_iname %in% drugs & pert_time=="24" & pert_dose_unit=="um" & pert_dose==10] # use only those with concentration 10 uM and treatment time 24 hrs, since these are the most frequent combinations in the dataset

# compute the population control for each assay plate in CMAP
# this is slow to run, mainly limited by the slow speed of reading the .gctx CMAP data file from disk (this file is too big to load into memory)
plts <- unique(instinfo1$rna_plate)
names(plts) <- plts
res <- sapply(plts, function(p) {
  ids <- instinfo[rna_plate==p, inst_id]
  mat <- parse_gctx("/data/Lab_ruppin/resources/l1000/PhaseI/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx", rid=gns$id, cid=ids)@mat
  rownames(mat) <- gns$symbol
  apply(mat, 1, median)
})

# save result: this output file was already generated and in the folder
#saveRDS(res, file="../data/population.ctrls.RDS")

