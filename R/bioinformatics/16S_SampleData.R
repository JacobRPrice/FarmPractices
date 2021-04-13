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
sd <- read.csv(file.path(data_path,"FST_metadata.csv"),sep=",",header=TRUE)


#----------------------------------------------------------
########
# format sample data
########
dim(sd)
names(sd)
str(sd)

# sd$Sample_ID
# sd$Seq.tube.ID
sd$SWRC_.ID <- as.numeric(sd$SWRC_.ID)
sd$Sample_or_Control <- as.factor(sd$Sample_or_Control)

sd$Sample.Site.ID
sd$Sample.Site.ID[c(97:99)] <- NA
sd$Sample.Site.ID

sd$Date.Collected
sd$Date.Collected[c(97:99)] <- NA
sd$Date.Collected

sd$Sample_depth_cm
sd$Sample_depth_cm[c(97:99)] <- NA
sd$Sample_depth_cm
sd$Sample_depth_cm <- as.factor(sd$Sample_depth_cm)
sd$Sample_depth_cm

# sd$Site
# sd$Sample_ID.1
# sd$Site_PE
# sd$Point_of_Entry
# sd$Replicate

sd$System
sd$System[c(97:99)] <- NA
sd$System
sd$System <- as.factor(sd$System)
sd$System


sd$Till_Type
sd$Till_Type[c(97:99)] <- NA
sd$Till_Type
sd$Till_Type <- as.factor(sd$Till_Type)
sd$Till_Type

sd$CoverCrop
sd$CoverCrop[c(97:99)] <- NA
sd$CoverCrop
sd$CoverCrop <- as.factor(sd$CoverCrop)
sd$CoverCrop

sd$CoverCrop_Type



#----------------------------------------------------------
########
# load track object
########
trackobj <- readRDS(file.path(rds_path_16S, "trackobj.RDS"))
# trackobj <- readRDS("/Users/jprice/Desktop/trackobj.RDS")
# trackobj
class(trackobj)
dim(trackobj)
trackobj <- as.data.frame(trackobj)
trackobj
sum(trackobj$reads.in)
sum(trackobj$reads.out)
trackobj$fastq <- rownames(trackobj)
rownames(trackobj) <- NULL
trackobj

trackobj$fastq

?strsplit
strsplit(
  x = trackobj$fastq,
  split = "_R"
)

#----------------------------------------------------------
########
# merge track obj and sample data 
########

###
# create vector to link entries in both dataframes
###
# metadata file contains samples ordered by (lab) tube number used during library prep. 
# fastq files were processed according to their (sorted)  filename order
# The result is that they're in different orders. 
# Since we are stitching together the phyloseq object according to the order in which the fastq files were processed through dada, we should make the metadata file conform to the order of the fastq files. 

paste0(
  sd$Site,
  "_",
  sd
)
names(sd)










