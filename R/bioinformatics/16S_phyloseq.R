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
saveRDS(ps, file.path(rds_path_16S, "ps-orig_16S.RDS"))
# also save in data directory 
saveRDS(ps, file.path(data_path, "ps-orig_16S.RDS"))

#----------------------------------------------------------
########
# decontam results
########
# https://astrobiomike.github.io/amplicon/dada2_workflow_ex
# https://github.com/benjjneb/decontam
# https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0605-2
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html