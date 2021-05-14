#----------------------------------------------------------
# Notes or References:
#
#----------------------------------------------------------
########
# Preliminary Items
########
# load utilities and functions
source(file.path(getwd(), "R/functions&utilities.R"))

# prep environment and get data
sapply(
  c("dada2", "phyloseq"), 
  require, 
  character.only=TRUE
)

# set seed
set.seed(20190611)

#----------------------------------------------------------
########
# load data & prep paths
########
taxtab<-readRDS(file.path(rds_path_16S,"taxtab.RDS"))
sd<-readRDS(file.path(rds_path_16S,"sd.RDS"))
seqtab<-readRDS(file.path(rds_path_16S,"seqtab.RDS"))

#----------------------------------------------------------
########
# ensure that samples are in the same/correct order
########
dim(sd)
dim(seqtab)

# these are sorted
rownames(seqtab)
# these are in numerical/counting order, not the same as above. 
names(sd)
sd$SRA_File_Name_1

# sort sd to match same order as in seqtab
sum(
  rownames(seqtab) ==
  sd[order(sd$SRA_File_Name_1),]$SRA_File_Name_1
)

sd <- sd[order(sd$SRA_File_Name_1),]

rownames(sd) <- sd[order(sd$SRA_File_Name_1),]$SRA_File_Name_1

#----------------------------------------------------------
########
# create phyloseq object
########
ps <- phyloseq(
  tax_table(taxtab),
  sample_data(sd),
  otu_table(seqtab, taxa_are_rows = FALSE)
)

ps 

# cleanup names
sample_names(ps)

cbind(
  sample_names(ps),
  sd$SRA_Sample_Name,
  sd$isNeg
)

sample_names(ps) <- sd$SRA_Sample_Name
sample_names(ps)

#----------------------------------------------------------
########
# save
########

# save in RDS directory
saveRDS(ps, file.path(rds_path_16S, "ps-orig_16S.RDS"))
# also save in data directory 
saveRDS(ps, file.path(data_path, "ps-orig_16S.RDS"))

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
#
########