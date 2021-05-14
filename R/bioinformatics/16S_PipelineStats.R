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
# load tracking objects
########

track_QC <-  readRDS(file.path(rds_path_16S, "track_QC.RDS"))
dim(track_QC)
str(track_QC)
track_QC

#----------------------------------------------------------
########
# Import dada2 pipeline results and collect statistics
########
dadaFs <- readRDS(file.path(rds_path_16S, "dadaFs.RDS"))
dadaRs <- readRDS(file.path(rds_path_16S, "dadaRs.RDS"))
mergers <- readRDS(file.path(rds_path_16S, "mergers.RDS"))
seqtab <- readRDS(file.path(rds_path_16S, "seqtab.RDS"))

getN <- function(x) sum(getUniques(x))

track_dada <- cbind(
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab)
)
colnames(track_dada) <- c(
  "denoisedF",
  "denoisedR",
  "merged",
  "nonchim"
)
track_dada <- as.data.frame(track_dada)
track_dada$file <- rownames(track_dada)

dim(track_dada)
str(track_dada)
track_dada

#----------------------------------------------------------
########
# Merge Tracking Objects
########
names(track_QC)
names(track_dada)

sum(track_dada$file %in% track_QC$file)
# we can merge using $file in each data frame
sum(track_dada$file == track_QC$file)
# the entries are in the same order

track <- merge(
  track_QC, track_dada, 
  by = "file"
)

# visually check that everything worked out the way we expected it to. 
head(track)
head(track_QC)
head(track_dada)

# assign rownames
rownames(track) <- track$file

# save results as RDS
saveRDS(track, file.path(rds_path_16S, "track.RDS"))
# also export as csv for easy reading outside of R
write.csv(track, file.path(output_path_16S, "track.csv"))

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
#
########