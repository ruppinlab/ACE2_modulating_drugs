##### This is the script for obtaining the CMAP sample info for a selected set of drugs and cells of interests, to be used for subsequent treated vs control differential expression analysis
# this output of this file is the CMAP sample info on the 989 approved drugs on the 28 CMAP cell types, used in our extended analysis

library(data.table)

# load info on drugs of interests: this loaded file contains information on clinically approved drugs collected from DrugBank
load("../data/drugs.of.interest.RData")
drugs <- d.app$pert_iname # all approved drugs

# load CMAP cell info: the data in this loaded file was downloaded from CMAP in plain text format and saved as an RData file
load("../data/cellinfo.RData")
cells <- cellinfo[grep("carcinoma",subtype), cell_id] # use only carcinoma cell lines, for the main part of our analysis

# load CMAP instance info (i.e. info on all available samples): the data in this loaded file was downloaded from CMAP in plain text format and saved as an RData file
load("../data/instinfo.RData")

# filter for the drugs, cell types, and treatment conditions of our interest; this represent all the 989 approved drugs available in CMAP
instinfo <- instinfo[pert_iname %in% drugs & pert_time=="24" & pert_dose_unit=="um" & pert_dose==10] # use only those with concentration 10 uM and treatment time 24 hrs, since these are the most frequent combinations in the dataset
instinfo <- instinfo[, .(phase="i", inst_id, plate=rna_plate, cell_id, pert_iname)]
instinfo <- instinfo[, .(inst=list(inst_id), plate=list(plate), cell=list(cell_id), phase=list(phase)), by=.(pert_iname)]
instinfo[, uniqueN(pert_iname)] # 989 drugs
instinfo[, uniqueN(unlist(cell))] # 28 cell types

# save selected sample info
#saveRDS(instinfo, file="../results/selected.samples989.RDS")


