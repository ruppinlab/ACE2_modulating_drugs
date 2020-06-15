##### This is the file for running the ACE2 differential expression analysis between drug-treated and control samples

# source the R script file containing necessary functions
source("cmap.de.funcs.R")

# load the previously generated data file containing samples for our selected drugs of interests
# in this example, here we use the samples corresponding to the 672 drugs tested on the four carcinoma cell lines, used in the main part of our study
instinfo <- readRDS("../results/selected.samples672.RDS")
# alternatively, for example, can use the samples corresponding to the 989 drugs from our extended analysis
#instinfo <- readRDS("../results/selected.samples989.RDS")

# run DE for all selected drugs
# this is slow, with speed mainly limited by reading the CMAP .gctx data file from disk (the file is too big to load into memory)
# in our study, instead of using the lapply() loop as below, we submitted an array job to a computational cluster
# you may also use mclapply(..., mc.cores=`number of cpu cores`) to parallelize the loop
de.res <- rbindlist(lapply(1:nrow(instinfo), function(i) {
  tryCatch(get.de(instinfo$inst[[i]], instinfo$plate[[i]], instinfo$cell[[i]], instinfo$phase[[i]]),
           error=function(e) data.table(id=NA,log.fc=NA,ave.expr=NA,t=NA,pval=NA,padj=NA,B=NA))
}))
de.res <- cbind(instinfo[,1], de.res)
de.res[, padj:=p.adjust(pval, "BH")]
de.res <- de.res[order(padj,pval)]

# save results
#saveRDS(de.res, file="../results/de.res672.RDS")

