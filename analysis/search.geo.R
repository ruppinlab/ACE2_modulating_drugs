##### This is the script used for a programmatic search of GEO database for candidate datasets involving drug treatment
# the resulting candidate datasets from this script were then manually curated and analyzed

library(data.table)
library(GEOmetadb)
library(stringr)


dbfile <- getSQLiteFile() # download database to current dir, last downloaded on 2020.05.16, database timestamp 2020.05.10
con <- dbConnect(SQLite(), dbfile)
initRegExp(con) # allow regexp in sqlite database

make.sql.string.gse <- function(...) {
  x <- list(...)
  res <- sapply(x, function(i) sprintf("(title regexp '%s' or summary regexp '%s' or overall_design regexp '%s')", i,i,i))
  res <- paste(res, collapse=" and ")
  paste("select gse,title,summary,overall_design from gse where", res, "and type like '%expression%'") # I checked that type like %expression% can select for all gene expression profiling datasets
}

make.sql.string.gds <- function(...) {
  x <- list(...)
  res <- sapply(x, function(i) sprintf("(title regexp '%s' or description regexp '%s')", i,i))
  res <- paste(res, collapse=" and ")
  paste("select gds,title,description from gds where", res, "and sample_organism regexp 'Homo sapiens|Mus musculus'")
}

qgeo <- function(...) {
  rs <- tryCatch(as.data.table(dbGetQuery(con, make.sql.string.gse(...))), error=function(e) NULL)
  if (!is.null(rs) && nrow(rs)!=0) rs <- rs[, .(id=gse, detail=sprintf("# Title: %s\n\n# Summary:\n  %s\n\n# Design:\n  %s\n\n",title,summary,overall_design))] else rs <- NULL
  rs1 <- tryCatch(as.data.table(dbGetQuery(con, make.sql.string.gds(...))), error=function(e) NULL)
  if (!is.null(rs1) && nrow(rs1)!=0) rs1 <- rs1[, .(id=gds, detail=sprintf("# Title: %s\n\n# Description:\n  %s\n\n",title,description))] else rs1 <- NULL
  if (!(is.null(rs) && is.null(rs1))) res <- rbind(rs, rs1) else res <- NULL
  res
}


# all approved drugs from drugbank

load("../data/drugs.of.interest.RData")
drugs <- d.app$pert_iname
drugs.rgx <- paste0("[",toupper(str_sub(drugs,1,1)),tolower(str_sub(drugs,1,1)),"]",str_sub(drugs,2,-1))

# search database

pb <- round(seq(0.01, 0.99, by = 0.01) * length(drugs)) # a progress bar for monitoring
message("0%")
tmp <- 1:length(drugs.rgx)
names(tmp) <- drugs
result <- lapply(tmp, function(i) {
  a <- match(i, pb)
  if (!is.na(a)) message(a, "%")
  qgeo(drugs.rgx[i])
}) # seems that cannot do it in parallel -- will lead to error querying the sql database
result <- result[!sapply(result, is.null)]
result <- rbindlist(result, idcol="drug")

# result is the raw output of the GEO database search for gene expression datasets involving treatment by approved drugs
result <- result[, .(drugs=list(unique(drug)), detail=detail[1]), by=id]
nrow(result) # 9636

# further filter for datasets of lung-derived samples
lung <- result[grep("[ (-][Ll]ung[ -]|[Pp]neumo|[Bb]ronch(ial|us)", detail)][grep("[Tt]reat", detail)][-grep("[ (-][Pp]lant[ -]|[Ss]wine|[Rr]abbit|[Bb]acteria|[Ff]etal", detail)]
nrow(lung) # [1] 248

# further filter for datasets of kidney-derived samples
kidney <- result[grep("[ (-][Kk]idney[ -]|[ (-][Nn]ephr|[ (-][Rr]enal", detail)][grep("[Tt]reat", detail)][-grep("[ (-][Pp]lant[ -]|[Ss]wine|[Rr]abbit|[Bb]acteria|[Ff]etal", detail)]
nrow(kidney) # [1] 148

# these resulting datasets were further manually curated to select the really relevant ones,
# then ACE2 differential expression was performed for each of the finally selected datasets manually, an example is given in geo.de.R


