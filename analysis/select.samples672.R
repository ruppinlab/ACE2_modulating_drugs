##### This is the script for obtaining the CMAP sample info for a selected set of drugs and cells of interests, to be used for subsequent treated vs control differential expression analysis
# this output of this file is the CMAP sample info on the 672 drugs on the 4 carcinoma cell lines, used in the main part of our analysis

library(data.table)

# load info on drugs of interests: this loaded file contains information on clinically approved drugs collected from DrugBank
load("../data/drugs.of.interest.RData")
drugs <- d.app$pert_iname # all approved drugs

# load CMAP cell info: the data in this loaded file was downloaded from CMAP in plain text format and saved as an RData file
load("../data/cellinfo.RData")
cells <- cellinfo[grep("carcinoma",subtype), cell_id] # use only carcinoma cell lines, for the main part of our analysis

# load CMAP instance info (i.e. info on all available samples): the data in this loaded file was downloaded from CMAP in plain text format and saved as an RData file
load("../data/instinfo.RData")

# filter for the drugs, cell types, and treatment conditions of our interest
instinfo <- instinfo[pert_iname %in% drugs & pert_time=="24" & pert_dose_unit=="um" & pert_dose==10 & cell_id %in% cells] # use only those with concentration 10 uM and treatment time 24 hrs, since these are the most frequent combinations in the dataset

# check the number of drugs tested on each cell type -- we want to find a set of cell types that a reasonable number of drugs were tested on all of these cells
( tmp <- instinfo[, .(n=uniqueN(pert_iname)), by=cell_id][order(-n)] )
# get the drugs that are commonly tested on the top 4 cell lines: this is the 672 drugs we focused on in the main part of our study
cells <- tmp[1:4, cell_id]
drugs <- Reduce(intersect, instinfo[cell_id %in% cells, .(x=list(unique(pert_iname))), by=cell_id]$x)
length(drugs) # 672 drugs

# get the drug-treated instance IDs for the selected drugs in the selected cells
instinfo <- instinfo[pert_iname %in% drugs & cell_id %in% cells]
instinfo <- instinfo[, .(phase="i", inst_id, plate=rna_plate, cell_id, pert_iname)]
instinfo <- instinfo[, .(inst=list(inst_id), plate=list(plate), cell=list(cell_id), phase=list(phase)), by=.(pert_iname)]

# save selected sample info
#saveRDS(instinfo, file="../results/selected.samples672.RDS")


