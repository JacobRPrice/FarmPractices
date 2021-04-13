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
  c("dada2", "phyloseq", "decontam"), 
  require, 
  character.only=TRUE
)

# set seed
set.seed(20190611)

#----------------------------------------------------------
########
# load data & prep paths
########
# seqtab <- readRDS(file.path(rds_path_16S, "seqtab.RDS"))
# sampdat <- readRDS(file.path(rds_path_16S, "sampdat.RDS"))

taxtab<-readRDS(file.path(rds_path_16S,"taxtab.RDS"))
samdf<-readRDS(file.path(rds_path_16S,"samdf.RDS"))
seqtab<-readRDS(file.path(rds_path_16S,"seqtab.RDS"))
# fitGTR<-readRDS(file.path(rds_path_16S,"fitGTR.RDS"))





#----------------------------------------------------------
########
# create phyloseq object
########
ps <- phyloseq(
  tax_table(taxtab),
  sample_data(samdf),
  otu_table(seqtab, taxa_are_rows = FALSE)#,
  # phy_tree(fitGTR$tree)
)

saveRDS(ps, file.path(rds_path_16S, "ps-orig.RDS"))

#----------------------------------------------------------
########
# decontam results
########
# https://astrobiomike.github.io/amplicon/dada2_workflow_ex
# https://github.com/benjjneb/decontam
# https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0605-2
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html