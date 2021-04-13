#----------------------------------------------------------
# Notes or References:
#
#----------------------------------------------------------
########
# Preliminary Items
########
getwd()
# load utilities and functions
source(file.path(getwd(), "R/functions&utilities.R"))

# prep environment and get data
sapply(
  c("dada2","ggplot2", "magrittr", "dplyr"), 
  require, 
  character.only=TRUE
)

# set seed
set.seed(20190611)

#----------------------------------------------------------
########
# import sequences
########
filt_path <- file.path(data_path,"ITS_filt")

filts <- sort(list.files(filt_path, full.names = TRUE))
filtFs <- filts[grepl("R1", filts)]
filtRs <- filts[grepl("R2", filts)]

head(filtFs)
head(filtRs)

#----------------------------------------------------------
########
# specify error function decision
########
# We will be using an error function that: 
#   loess() : set span = 2 and implement changes to weights
#   enforce monotonicity

loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution 
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity 
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>% 
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

#----------------------------------------------------------
########
# dada - FORWARD
########
dadaFs <- dada(
  derep = filtFs,
  err = NULL,
  selfConsist = TRUE,
  pool = FALSE,
  errorEstimationFunction = loessErrfun_mod,
  multithread = TRUE, 
  verbose = 2
)
dada2:::checkConvergence(dadaFs[[1]])
saveRDS(dadaFs, file.path(rds_path_ITS, "dadaFs.RDS"))
ggsave(
  plotErrors(dadaFs, nominalQ=TRUE),
  file.path(figs_path_ITS,"dada_errF.png"),
  width = 7, height = 7, units = "in"
)

#----------------------------------------------------------
########
# dada - REVERSE
########
dadaRs <- dada(
  derep = filtRs,
  err = NULL,
  selfConsist = TRUE,
  pool = FALSE,
  errorEstimationFunction = loessErrfun_mod,
  multithread = TRUE, 
  verbose = 2
)
dada2:::checkConvergence(dadaRs[[1]])
saveRDS(dadaRs, file.path(rds_path_ITS, "dadaRs.RDS"))
ggsave(
  plotErrors(dadaRs, nominalQ=TRUE),
  file.path(figs_path_ITS,"dada_errR.png"),
  width = 7, height = 7, units = "in"
)

#----------------------------------------------------------
########
# merge paired reads
########
mergers <- mergePairs(
  dadaFs, 
  filtFs, 
  dadaRs, 
  filtRs,
  # minOverlap=20,
  verbose=TRUE
)
names(mergers)
table(mergers[[1]]$nmatch)
head(mergers[[1]])
saveRDS(mergers,file.path(rds_path_ITS,"mergers.RDS"))

#----------------------------------------------------------
########
# Sequence Table
########
seqtab.all <- makeSequenceTable(mergers)
dim(seqtab.all)
table(nchar(getSequences(seqtab.all)))
saveRDS(seqtab.all, file.path(rds_path_ITS, "seqtab.all.RDS"))


#----------------------------------------------------------
########
# remove chimeras
########
seqtab <- removeBimeraDenovo(
  seqtab.all,
  multithread = TRUE, 
  verbose = TRUE
)
sum(seqtab)/sum(seqtab.all)
saveRDS(seqtab, file.path(rds_path_ITS, "seqtab.RDS"))


#----------------------------------------------------------
########
# collect statistics
########
# trackobj <- readRDS(file.path(rds_path_ITS, "trackobj.RDS"))
# 
# head(trackobj)
# 
# getN <- function(x) sum(getUniques(x))
# track <- cbind(
#   trackobj, 
#   sapply(dadaFs, getN), 
#   sapply(dadaRs, getN), 
#   sapply(mergers, getN), 
#   rowSums(seqtab)
# )
# colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
# head(track)
# 
# saveRDS(track, file.path(rds_path_ITS, "track.RDS"))

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########