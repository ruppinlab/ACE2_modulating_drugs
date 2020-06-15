##### this file contains functions for running differential expression (DE) analysis for datasets from GEO

library(data.table)
library(stringr)
library(GEOquery)
library(affy)
library(limma)

prep.data <- function(dat, log="default", norm.method="loess") {

  if (class(dat)=="ExpressionSet") {
    # for featureData, fix potential improper column names so that later limma::topTable can use them
    fvarLabels(dat) <- make.names(fvarLabels(dat))
    mat <- exprs(dat)
  } else if (is.matrix(dat)) mat <- dat

  # log2 transform
  if (log=="default") {
    # following the method as in GEO2R, for microarray data
    qx <- as.numeric(quantile(mat, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm=TRUE))
    log <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  }
  if (log) {
    nlt0 <- sum(mat<0, na.rm=TRUE)
    if (nlt0>0) {
      warning(sprintf("Log-transformation error: there are %d negative values in the data, the data may be already on log-scale.\nHere is a summary of the data:\n", nlt0))
      print(summary(as.vector(mat)))
      stop()
    }
    mat <- log2(mat+1)
    cat("log2-transformation performed.\n")
  } else cat("log2-transformation NOT performed.\n")

  # normalization
  if (norm.method=="loess") {
    nna <- sum(is.na(mat))
    if (nna>0) {
      stop(sprintf("Loess normalization error: there are %d NA/NaN's in the data.\n", nna))
    } else mat <- affy::normalize.loess(mat, log.it=FALSE)
  } else if (norm.method=="quantile") {
    mat <- limma::normalizeQuantiles(mat)
  } else cat("Normalization NOT performed.\n")

  # return
  if (class(dat)=="ExpressionSet") {
    exprs(dat) <- mat
    return(dat)
  } else if (is.matrix(dat)) return(mat)
}


de <- function(dat, pheno, model="~.", coef, rn=NULL, robust=FALSE, trend=FALSE) {

  if (class(dat)=="ExpressionSet") {
    mat <- exprs(dat)
    idx <- which(tolower(names(fData(dat))) %in% c("gene.symbol","gene_symbol","gene symbol","symbol","orf",rn))
    if (length(idx)!=1) stop("Issue with gene symbols, please check.")
    rownames(mat) <- fData(dat)[[idx]]
  } else if (is.matrix(dat)) mat <- dat

  design <- model.matrix(as.formula(model), pheno)
  fit <- limma::lmFit(mat, design)
  fit <- limma::eBayes(fit, robust=robust, trend=trend)
  res <- tryCatch({
    tt <- as.data.table(limma::topTable(fit, coef=coef, number=Inf, genelist=rownames(mat)))
    setnames(tt, c("id","log.fc","ave.expr","t","pval","padj","B"))
    tt
  }, error=function(e) {
    tt <- as.data.table(limma::topTable(fit, coef=coef, number=Inf))
    setnames(tt, c("id","log.fc","ave.expr","t","pval","padj","B"))
    tt
  })

  res
}
 
