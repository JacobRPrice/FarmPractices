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
taxtab<-readRDS(file.path(rds_path_ITS,"taxtab.RDS"))
sd<-readRDS(file.path(rds_path_ITS,"sd.RDS"))
seqtab<-readRDS(file.path(rds_path_ITS,"seqtab.RDS"))

dim(sd)
dim(seqtab)
sum(rownames(seqtab) == rownames(sd))


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

# save in RDS directory
saveRDS(ps, file.path(rds_path_ITS, "ps-orig_ITS.RDS"))
# also save in data directory 
saveRDS(ps, file.path(data_path, "ps-orig_ITS.RDS"))

#----------------------------------------------------------
########
# decontam results
########
# https://astrobiomike.github.io/amplicon/dada2_workflow_ex
# https://github.com/benjjneb/decontam
# https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0605-2
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html