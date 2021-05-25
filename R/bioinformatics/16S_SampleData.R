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
# sapply(
#   c("dada2"),
#   require,
#   character.only=TRUE
# )

# set seed
set.seed(20190611)

#----------------------------------------------------------
########
# load sample data 
########
dat <- read.csv(file.path(data_path,"metadata.csv"),sep=",",header=TRUE)
dim(dat)
names(dat)
str(dat)
# are there any N/A left over from exporting the csv?
# which(dat == "N/A")
# which(dat == "NA")

# ensure variables are coded correctly
dat$Date_Collected <- as.Date(dat$Date_Collected, format = "%m/%d/%y")
dat$Sample_Depth <- factor(dat$Sample_Depth)
dat$ID_Rep <- factor(dat$ID_Rep)
dat$ID_Sys <- factor(dat$ID_Sys)
dat$ID_EP <- factor(dat$ID_EP)
dat$ID_subplot <- factor(dat$ID_subplot)
dat$Plot <- factor(dat$Plot)
dat$Fertility_Source <- factor(dat$Fertility_Source)
dat$Management_System <- factor(dat$Management_System)
dat$Pesticide_Application <- factor(dat$Pesticide_Application)
dat$Tillage <- factor(dat$Tillage)
dat$Cover_Crop <- factor(dat$Cover_Crop)
dat$K <- as.numeric(dat$K)
dat$Mg <- as.numeric(dat$Mg)
dat$P <- as.numeric(dat$P)

str(dat)

#----------------------------------------------------------
########
# subset sample data to include only the relevant entries
########
dim(dat)
dat <- subset(dat, Gene == "16S")
dim(dat)

#----------------------------------------------------------
########
# save sample data
########

# save results as RDS
saveRDS(dat, file.path(rds_path_16S, "sd.RDS"))
# also export as csv for easy reading outside of R
write.csv(dat, file.path(output_path_16S, "sd.csv"))

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
#
########
