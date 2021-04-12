#!/bin/bash

###
# copy raw data and other contextual information to project data directory 
###

###
# 16S dada2 pipeline
###
date
Rscript ./R/bioinformatics/16S_QC.R

date
Rscript ./R/bioinformatics/16S_dada.R

date
Rscript ./R/bioinformatics/16S_Taxonomy.R

# date
#Rscript ./R/bioinformatics/16S_SampleData.R

# date
#Rscript ./R/bioinformatics/16S_Phyloseq.R

###
# ITS dada2 pipeline
###

# date
# Rscript ./R/bioinformatics/ITS_QC.R

# date
# Rscript ./R/bioinformatics/ITS_dada.R

# date
# Rscript ./R/bioinformatics/ITS_Taxonomy.R

# date
#Rscript ./R/bioinformatics/ITS_SampleData.R

# date
#Rscript ./R/bioinformatics/ITS_Phyloseq.R
