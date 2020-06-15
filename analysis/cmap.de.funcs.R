##### this file contains functions for running differential expression (DE) analysis of drug-treated vs population control samples using Level 3 LINCS l1000 data
# The Level 3 LINCS l1000 data should be downloaded from the GEO database to the "data" folder under the working directory, the GSE dataset IDs and the file names are as follows:
# Phase 1 data: data/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx
# Phase 2 data: data/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx

library(data.table)
library(cmapR)
library(limma)

# get best inferred genes (BING, this includes the landmark genes)
load("../data/geneinfo.RData")
gns <- geneinfo[pr_is_bing==1, .(id=as.character(pr_gene_id), symbol=pr_gene_symbol)] # use only the sets of BING (landmark and best-inferred) genes

# load previously generated population controls, i.e. median expression for each plate
pop1 <- readRDS("../data/population.ctrls.RDS")

get.mat <- function(phase, insts) {
  # a function to read data from the l1000 gctx file
  if (phase==1) {
    fn <- "../data/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"
  } else if (phase==2) {
    fn <- "../data/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx"
  }
  mat <- parse_gctx(fn, rid=gns$id, cid=insts)@mat
  rownames(mat) <- gns$symbol
  mat
}

de <- function(dat, pheno, model="~.", coef, robust=FALSE, trend=FALSE) {
  # a function to perform DE analysis with limma
  mat <- dat
  design <- model.matrix(as.formula(model), pheno)
  fit <- lmFit(mat, design)
  fit <- eBayes(fit, robust=robust, trend=trend)
  res <- tryCatch({
    tt <- as.data.table(topTable(fit, coef=coef, number=Inf, genelist=rownames(mat)))
    setnames(tt, c("id","log.fc","ave.expr","t","pval","padj","B"))
    tt
  }, error=function(e) {
    tt <- as.data.table(topTable(fit, coef=coef, number=Inf))
    setnames(tt, c("id","log.fc","ave.expr","t","pval","padj","B"))
    tt
  })
  res
}

get.de <- function(insts, plates, cells, phases, gn="ACE2") {
  # a function to perform DE analysis for a given gene (gn), for treated vs population control samples of one particular drug
  # insts: the l1000 instance IDs of the treated samples for the drug
  # plates, cells, and phases are the plates, cell lines, and data phases corresponding to insts (in the same order)
  # the DE analysis will be performed for the drug-treated vs control, while controlling for plates, cells and phases as covariates in a linear model with limma
  mat1 <- NULL
  phe1 <- NULL
  mat2 <- NULL
  phe2 <- NULL
  p1 <- phases=="i"
  if (any(p1)) {
    mat.trt <- get.mat(1, insts[p1])
    mat.ctl <- pop1[, plates[p1]]
    mat1 <- cbind(mat.trt, mat.ctl)
    phe1 <- data.table(grp=factor(rep(c("trt","ctl"), c(ncol(mat.trt), ncol(mat.ctl))), levels=c("ctl","trt")), repl=rep(plates[p1],2), cell=rep(cells[p1],2), phase="i")
  }
  p2 <- phases=="ii"
  if (any(p2)) {
    mat.trt <- get.mat(2, insts[p2])
    mat.ctl <- pop2[, plates[p2]]
    mat2 <- cbind(mat.trt, mat.ctl)
    phe2 <- data.table(grp=factor(rep(c("trt","ctl"), c(ncol(mat.trt), ncol(mat.ctl))), levels=c("ctl","trt")), repl=rep(plates[p2], 2), cell=rep(cells[p2],2), phase="ii")
  }
  mat <- cbind(mat1, mat2)
  phe <- rbind(phe1, phe2)
  if (uniqueN(phe$cell)==1) phe[, cell:=NULL]
  if (uniqueN(phe$phase)==1) phe[, phase:=NULL]
  de.res <- de(mat, phe, coef="grptrt")
  de.res[id==gn]
}


