source("de.funcs.R")


# drugs
load("data/drugs.of.interest.RData")
drugs <- d.app$pert_iname # all approved drugs

# cell info
load("data/cellinfo.RData")
cells <- cellinfo[grep("carcinoma",subtype), cell_id] # use only carcinoma cell lines

# try to get a common set of cell lines that a decent fraction of drugs are being tested on
load("data/instinfo.RData") # p1
instinfo <- instinfo[pert_iname %in% drugs & pert_time=="24" & pert_dose_unit=="um" & pert_dose==10 & cell_id %in% cells] # use only those with concentration 10 uM and treatment time 24 hrs, since these are the most frequent combinations in the dataset
# check the number of drugs tested on each cell linear
tmp <- instinfo[, .(n=uniqueN(pert_iname)), by=cell_id][order(-n)]
# get the drugs that are commonly tested on the top 4 cell lines
cells <- tmp[1:4, cell_id]
drugs <- Reduce(intersect, instinfo[cell_id %in% cells, .(x=list(unique(pert_iname))), by=cell_id]$x)
length(drugs) # 672 drugs
# get the drug-treated instance IDs for the selected drugs in the selected cells
instinfo <- instinfo[pert_iname %in% drugs & cell_id %in% cells]
instinfo <- instinfo[, .(phase="i", inst_id, plate=rna_plate, cell_id, pert_iname)]
instinfo <- instinfo[, .(inst=list(inst_id), plate=list(plate), cell=list(cell_id), phase=list(phase)), by=.(pert_iname)]

# run DE for all selected drugs
de.res <- rbindlist(lapply(1:nrow(instinfo), function(i) {
  tryCatch(get.de(instinfo$inst[[i]], instinfo$plate[[i]], instinfo$cell[[i]], instinfo$phase[[i]]),
           error=function(e) data.table(id=NA,log.fc=NA,ave.expr=NA,t=NA,pval=NA,padj=NA,B=NA))
}))
de.res <- cbind(instinfo[,1], de.res)
de.res[, padj:=p.adjust(pval, "BH")]
de.res <- de.res[order(padj,pval)]

saveRDS(de.res, file="de.res.RDS")

