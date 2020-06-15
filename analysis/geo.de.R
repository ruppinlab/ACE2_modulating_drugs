##### This is an example script for performing differential expression analysis of drug-treated vs control samples from selected GEO datasets
# in this example, the following GEO dataset is used:
geo.id <- "GSE33561"

# source the functions required for doing DE
source("geo.de.funcs.R")

# get data from GEO by ID
geo <- getGEO(GEO=geo.id, GSEMatrix=TRUE, AnnotGPL=TRUE)
# usually there is only 1 dataset within the GEO ID, but if there are more than one then need to check each manually
if (length(geo)==1) geo <- geo[[1]] else message(length(geo), " datasets available.")

# examine the phenotypic data
pData(geo)
## based on the phenotypic data:
# 1. get the drug used in the experiment
drug <- "losartan"
# 2. obtain the subset of relevant samples by their indices
tmp <- pData(geo)[["agent:ch1"]]
tmp1 <- pData(geo)[["batch:ch1"]]
idx <- c(which(tmp=="Room Air plus Losartan"), which(tmp=="Room Air" & tmp1=="batch 1"))
# 3. create a data.table of control/treated labels
tmp[idx] # check the selected samples
phe <- data.table(group=factor(rep(c("trt","ctrl"),each=3), levels=c("ctrl","trt")))

# pre-process data (selecting the relevant subset of samples by their indicies); this workflow is for microarray
dat <- prep.data(geo[,idx], norm.method="loess")
# DE analysis
de.res <- de(dat, phe, coef="grouptrt")
# check ACE2 (or Ace2 if mice)
( de.res.ace2 <- de.res[id=="Ace2"] )

# save result
#save(geo, idx, dat, phe, de.res, de.res.ace2, file=paste0("../results/",geo.id,"_",drug,".RData"))


