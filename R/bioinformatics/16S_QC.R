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
  c("dada2", "ShortRead", "Biostrings", "ggplot2"), 
  require, 
  character.only=TRUE
)

# set seed
set.seed(20190611)

#----------------------------------------------------------
########
# raw reads
########
rawread_path<-file.path(data_path, "16S")

fns <- sort(list.files(rawread_path, full.names = TRUE))
length(fns)

fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]
head(fns)
head(fnFs)
head(fnRs)

#----------------------------------------------------------
########
# pre-filter read files to remove any with "N" in their sequences
########
filt.filtN_path <- file.path(data_path, "16S_filtN")
if(!file_test("-d", filt.filtN_path)) dir.create(filt.filtN_path)

fnFs.filtN <- file.path(filt.filtN_path, basename(fnFs))
fnRs.filtN <- file.path(filt.filtN_path, basename(fnRs))

out.filtN <- filterAndTrim(
  fwd = fnFs,
  filt = fnFs.filtN, 
  rev = fnRs,
  filt.rev = fnRs.filtN,
  maxN = 0,
  multithread = parallel::detectCores() - 1
)
out.filtN

saveRDS(out.filtN, file.path(rds_path_16S, "out.filtN.RDS"))

#----------------------------------------------------------
########
# identify location and orientation of primers in the reads
########
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "GGACTACNVGGGTWTCTAAT"
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# identify/count primer appearances within forward and reverse reads
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
#                   Forward Complement Reverse RevComp
# FWD.ForwardReads  317252          0       0       0
# FWD.ReverseReads       0          0       0     350
# REV.ForwardReads       3          0       0     858
# REV.ReverseReads  313855          0       0       2

# orientation of reads are normal, with only a few readthrough


#----------------------------------------------------------
########
# use cutadapt to remove library prep primers
########
cutadapt <- "/Users/jprice/Library/Python/3.9/bin/cutadapt"
system2(cutadapt, args = "--version")

# create output filenames/paths for cutadapt-ed fastq files
cut_path <- file.path(data_path, "16S_filtN_cut")
if (!file_test("-d", cut_path)) dir.create(cut_path)

fnFs.cut <- file.path(cut_path, basename(fnFs.filtN))
fnRs.cut <- file.path(cut_path, basename(fnRs.filtN))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(
    cutadapt, 
    args = c(
      R1.flags, 
      R2.flags,
      # number of CPU cores to use
      "-j", 7,
      # -n 2 required to remove FWD and REV from reads
      "-n", 2, 
      # set minimum length so empty reads are not printed. This can interfere with other steps down the road (at least with visualizing the quality profiles.). 
      # The zero length reads to appear to be removed by the trim and filter step but it's better to be cautious and do this task explicitly here.
      # https://github.com/benjjneb/dada2/issues/159
      # https://cutadapt.readthedocs.io/en/stable/guide.html#filtering-reads
      "--minimum-length", 1,
      # output files
      "-o", fnFs.cut[i], 
      "-p", fnRs.cut[i], 
      # input files
      fnFs.filtN[i], 
      fnRs.filtN[i]
    )
  )
}

# sanity check
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
# damn-near perfect. 
# Forward Complement Reverse RevComp
# FWD.ForwardReads       0          0       0       0
# FWD.ReverseReads       0          0       0       0
# REV.ForwardReads       3          0       0       0
# REV.ReverseReads       0          0       0       2

#----------------------------------------------------------
########
# check read quality post-cutadapt
########
cutFs <- sort(list.files(cut_path, pattern = "_R1.fq.gz", full.names = TRUE))
cutRs <- sort(list.files(cut_path, pattern = "_R2.fq.gz", full.names = TRUE))

pF <- plotQualityProfile(cutFs[1:6]) +
  geom_vline(xintercept = 201, color = "red")
ggsave(
  file.path(figs_path_16S,"ReadQuality_post-cutadapt_F.png")
  )

pR <- plotQualityProfile(cutRs[1:6]) +
  geom_vline(xintercept = 201, color = "red")
ggsave(
  file.path(figs_path_16S,"ReadQuality_post-cutadapt_R.png")
)

#----------------------------------------------------------
########
# trim and filter 
########
# filtered read directory
filt_path<-file.path(data_path,"16S_filt")
if(!file_test("-d",filt_path)) dir.create(filt_path)

filtFs<-file.path(filt_path,basename(cutFs))
filtRs<-file.path(filt_path,basename(cutRs))

# carry out trimming and filtering of reads
trackobj <- filterAndTrim(
  fwd=cutFs,
  filt=filtFs,
  rev=cutRs,
  filt.rev=filtRs,
  # trimLeft=10,
  truncLen=c(200,200),
  maxN=0,
  maxEE=2,
  truncQ=2, # default 2
  # rm.lowcomplex = 4, # remove low complexity sequences. 
  compress=TRUE,
  verbose=TRUE,
  multithread = parallel::detectCores() - 1
)

trackobj

saveRDS(trackobj, file.path(rds_path_16S, "trackobj.RDS"))

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########