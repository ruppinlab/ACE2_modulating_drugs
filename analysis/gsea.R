##### This is the script for analyzing the drug class enrichment in ACE2 regulators with the GSEA method
# here we use the drug classes as defined based on the WHO ATC drug classification

library(data.table)
library(fgsea)

# load previously generated ACE2 differential results: our results presented in the manuscript were based on the extended DE results of the 989 approved drugs across 28 CMAP cell lines
de.res <- readRDS("../results/de.res989.RDS")

# load WHO ATC drug indication data obtained from DrugBank
load("../data/drug.atc.RData")
# define drug classes based on WHO ATC
dset <- lapply(atc.drug.set, tolower)

# GSEA analysis: note that the exact numerical P values and adjusted P values may have minor changes due to that they are computed via a permutation test involving random sampling
x <- de.res$log.fc
names(x) <- de.res$pert_iname
gsea.res <- fgsea(dset, x, 1e4)
gsea.res <- cbind(atc.anno[match(gsea.res$pathway,atc.code), .(name=atc.term)], gsea.res)

# focusing on level 2 ATC anntation (i.e. with three letters/numbers in the ATC code) and with drug classes containing more than 10 drugs
gsea.res <- gsea.res[nchar(pathway)==3 & size>10][, padj:=p.adjust(pval, "BH")][order(padj,pval)]

# save result
#saveRDS(gsea.res, file="../results/gsea.res.atc.RDS")


