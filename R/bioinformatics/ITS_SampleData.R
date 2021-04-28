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
# load sample data 
########
dat <- datBAK <- read.csv(file.path(data_path,"metadata.csv"),sep=",",header=TRUE)
dim(dat)
names(dat)
str(dat)
# are there any N/A left over from exporting the csv?
# which(dat == "N/A")

# ensure variables are coded correctly
dat$Date_Collected <- as.Date(dat$Date_Collected, format = "%m/%d/%y")
dat$Sample_Depth <- factor(dat$Sample_Depth)
dat$Field_Rep <- factor(dat$Field_Rep)
dat$Field_Sys <- factor(dat$Field_Sys)
dat$Field_EP <- factor(dat$Field_EP)
dat$Field <- factor(dat$Field)
dat$System <- factor(dat$System)
dat$Till_Type <- factor(dat$Till_Type)
dat$System2 <- factor(dat$System2)
dat$Cover_Crop <- factor(dat$Cover_Crop)

str(dat)

#----------------------------------------------------------
########
# subset sample data to include only the relevant entries
########
dat <- subset(dat, Gene == "ITS2")
dim(dat)

#----------------------------------------------------------
########
# load tracking objects
########

track_QC <-  readRDS(file.path(rds_path_ITS, "track_QC.RDS"))
dim(track_QC)
str(track_QC)
track_QC

#----------------------------------------------------------
########
# Import dada2 pipeline results and collect statistics
########
dadaFs <- readRDS(file.path(rds_path_ITS, "dadaFs.RDS"))
dadaRs <- readRDS(file.path(rds_path_ITS, "dadaRs.RDS"))
mergers <- readRDS(file.path(rds_path_ITS, "mergers.RDS"))
seqtab <- readRDS(file.path(rds_path_ITS, "seqtab.RDS"))

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
saveRDS(track, file.path(rds_path_ITS, "track.RDS"))
# also export as csv for easy reading outside of R
write.csv(track, file.path(output_path_ITS, "track.csv"))


#----------------------------------------------------------
########
# Merge (the merged) tracking object with sample data
########
names(dat)
names(track)

# we can merge by the file names
sum(track$file %in% dat$SRA_File_Name_1)
# but the entries are not in the same order. 
# The entries in the sample data file need to be in the same order that the dada results are in!
dat$file <- dat$SRA_File_Name_1

sd <- merge(
  track,
  dat, 
  by = "file"
)
rownames(sd) <- sd$file

dim(dat)
dim(track)
dim(sd)
names(sd)

head(sd,11)
# everything looks good to go. 

# save results as RDS
saveRDS(sd, file.path(rds_path_ITS, "sd.RDS"))
# also export as csv for easy reading outside of R
write.csv(sd, file.path(output_path_ITS, "sd.csv"))

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
#
########