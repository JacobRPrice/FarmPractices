#----------------------------------------------------------
#----------------------------------------------------------
# Notes or References:
# 
#----------------------------------------------------------
########
# define filepaths and create directory structure
########

### 
# project-wide paths
proj_path<-getwd()
data_path<-file.path(proj_path,"data")

### 
# 16S-specific paths
figs_path_16S<-file.path(proj_path,"figs_16S")
if(!file_test("-d",figs_path_16S)) dir.create(figs_path_16S)

output_path_16S <- file.path(proj_path, "output_16S")
if(!file_test("-d", output_path_16S)) dir.create(output_path_16S)

rds_path_16S<-file.path(output_path_16S,"RDS_16S")
if(!file_test("-d", rds_path_16S)) dir.create(rds_path_16S)


# ### 
# # ITS2-specific paths
# figs_path_ITS<-file.path(proj_path,"figs_ITS")
# if(!file_test("-d",figs_path_ITS)) dir.create(figs_path_ITS)
# 
# output_path_ITS <- file.path(proj_path, "output_ITS")
# if(!file_test("-d", output_path_ITS)) dir.create(output_path_ITS)
# 
# rds_path_ITS<-file.path(output_path_ITS,"RDS_ITS")
# if(!file_test("-d", rds_path_ITS)) dir.create(rds_path_ITS)


# figs_path<-file.path(proj_path,"figs")
# if(!file_test("-d",figs_path)) dir.create(figs_path)
# output_path<-file.path(proj_path,"output")
# if(!file_test("-d",output_path)) dir.create(output_path)
# rds_path<-file.path(output_path,"RDS")
# if(!file_test("-d",rds_path)) dir.create(rds_path)

#----------------------------------------------------------
########
# Set Seed and Options
########
set.seed(20190611)

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########