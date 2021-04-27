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
  c("dada2"), 
  require, 
  character.only=TRUE
)

# set seed
set.seed(20190611)

#----------------------------------------------------------
########
# load data & prep paths
########
seqtab <- readRDS(file.path(rds_path_16S, "seqtab.RDS"))

#----------------------------------------------------------
########
# assign taxonomy 
########
taxtab.1to6 <- assignTaxonomy(
  seqs = seqtab,
  refFasta = file.path(
    data_path, "/SILVA/silva_nr99_v138.1_train_set.fa.gz"),
  minBoot = 80,
  multithread = TRUE,
  verbose = TRUE
)

saveRDS(
  taxtab.1to6,
  file.path(rds_path_16S, "taxtab.1to6.RDS")
)

dim(taxtab.1to6)
head(unname(taxtab.1to6))
colnames(taxtab.1to6)

#----------------------------------------------------------
########
# add species
########
taxtab <- addSpecies(
  taxtab = taxtab.1to6,
  refFasta = file.path(
    data_path, "/SILVA/silva_species_assignment_v138.1.fa.gz"),
  allowMultiple = TRUE,
  tryRC = FALSE,
  verbose = TRUE
)

saveRDS(
  taxtab,
  file.path(rds_path_16S, "taxtab.RDS")
)

head(unname(taxtab))
colnames(taxtab)
head(rownames(taxtab))
table(unname(taxtab)[,1],useNA="always")
table(unname(taxtab)[,2],useNA="always")
table(unname(taxtab)[,3],useNA="always")
table(unname(taxtab)[,4],useNA="always")
table(unname(taxtab)[,5],useNA="always")
table(unname(taxtab)[,6],useNA="always")
table(unname(taxtab)[,7],useNA="always")

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########