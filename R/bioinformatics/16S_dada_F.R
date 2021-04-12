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
  c("dada2","ggplot2"), 
  require, 
  character.only=TRUE
)

# set seed
set.seed(20190611)

#----------------------------------------------------------
########
# import sequences
########
filt_path <- file.path(data_path,"16S_filt")

filts <- sort(list.files(filt_path, full.names = TRUE))
filtFs <- filts[grepl("R1", filts)]
filtRs <- filts[grepl("R2", filts)]

head(filtFs)
head(filtRs)

# #----------------------------------------------------------
# ########
# # specify error model
# ########
# # binned quality scores require alterations to the error model in order to achieve a better (and more realistic) fit. 
# # see dada issue 938 for some context and the function we are using. 
# # https://github.com/benjjneb/dada2/issues/938
# 
# loessErrfun_mod <- function (trans) {
#   qq <- as.numeric(colnames(trans))
#   est <- matrix(0, nrow = 0, ncol = length(qq))
#   for (nti in c("A", "C", "G", "T")) {
#     for (ntj in c("A", "C", "G", "T")) {
#       if (nti != ntj) {
#         errs <- trans[paste0(nti, "2", ntj), ]
#         tot <- colSums(trans[paste0(nti, "2", c("A",
#                                                 "C", "G", "T")), ])
#         rlogp <- log10((errs + 1)/tot)
#         rlogp[is.infinite(rlogp)] <- NA
#         df <- data.frame(q = qq, errs = errs, tot = tot,
#                          rlogp = rlogp)
#         mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
#         pred <- predict(mod.lo, qq)
#         maxrli <- max(which(!is.na(pred)))
#         minrli <- min(which(!is.na(pred)))
#         pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
#         pred[seq_along(pred) < minrli] <- pred[[minrli]]
#         est <- rbind(est, 10^pred)
#       }
#     }
#   }
#   MAX_ERROR_RATE <- 0.25
#   MIN_ERROR_RATE <- 1e-07
#   est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
#   est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
#   err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ], est[4,
#   ], 1 - colSums(est[4:6, ]), est[5:6, ], est[7:8, ], 1 -
#     colSums(est[7:9, ]), est[9, ], est[10:12, ], 1 - colSums(est[10:12,
#     ]))
#   rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4),
#                           "2", c("A", "C", "G", "T"))
#   colnames(err) <- colnames(trans)
#   return(err)
# }



# #----------------------------------------------------------
# ########
# # specify error model
# ########
# # now we'll define a new error model, building off the idea in the modified function above
# # 
# # Try 1
# #   loess() : set span = 2 and implement changes to weights
# #   enforce monotonicity
# 
# library(magrittr)
# library(dplyr)
# 
# loessErrfun_mod <- function(trans) {
#   qq <- as.numeric(colnames(trans))
#   est <- matrix(0, nrow=0, ncol=length(qq))
#   for(nti in c("A","C","G","T")) {
#     for(ntj in c("A","C","G","T")) {
#       if(nti != ntj) {
#         errs <- trans[paste0(nti,"2",ntj),]
#         tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
#         rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
#         rlogp[is.infinite(rlogp)] <- NA
#         df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
#         
#         # original
#         # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
#         # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
#         # #        mod.lo <- loess(rlogp ~ q, df)
#         
#         # Gulliem Salazar's solution 
#         # https://github.com/benjjneb/dada2/issues/938
#         mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
#         
#         pred <- predict(mod.lo, qq)
#         maxrli <- max(which(!is.na(pred)))
#         minrli <- min(which(!is.na(pred)))
#         pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
#         pred[seq_along(pred)<minrli] <- pred[[minrli]]
#         est <- rbind(est, 10^pred)
#       } # if(nti != ntj)
#     } # for(ntj in c("A","C","G","T"))
#   } # for(nti in c("A","C","G","T"))
#   
#   # HACKY
#   MAX_ERROR_RATE <- 0.25
#   MIN_ERROR_RATE <- 1e-7
#   est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
#   est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
#   
#   # enforce monotonicity 
#   # https://github.com/benjjneb/dada2/issues/791
#   estorig <- est
#   est <- est %>% 
#     data.frame() %>%
#     mutate_all(funs(case_when(. < X40 ~ X40,
#                               . >= X40 ~ .))) %>% as.matrix()
#   rownames(est) <- rownames(estorig)
#   colnames(est) <- colnames(estorig)
#   
#   # Expand the err matrix with the self-transition probs
#   err <- rbind(1-colSums(est[1:3,]), est[1:3,],
#                est[4,], 1-colSums(est[4:6,]), est[5:6,],
#                est[7:8,], 1-colSums(est[7:9,]), est[9,],
#                est[10:12,], 1-colSums(est[10:12,]))
#   rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
#   colnames(err) <- colnames(trans)
#   # Return
#   return(err)
# }
# 
# # check what this looks like
# errF <- learnErrors(
#   filtFs, 
#   multithread = TRUE, 
#   errorEstimationFunction = loessErrfun_mod,
#   verbose = TRUE
# )
# dada2:::checkConvergence(errF)
# plotErrors(errF, nominalQ=TRUE)
# ggsave(
#   file.path(figs_path_16S,"learnErrors_errF_span2_enforceMono.pdf"),
#   width = 7, height = 7, units = "in"
# )
# saveRDS(errF, file.path(rds_path_16S, "errF_span2_enforceMono.RDS"))
# 
# # conclusions:
# # this looks to work relatively well. The curve is smooth without the hooks/dips that we see in the default approach. 
# # It also does not have the increasing A2G set like we saw with Guillian's approach of setting span = 2. 



# #----------------------------------------------------------
# ########
# # specify error model
# ########
# # Try 2
# #   loess() : keep original call
# #   enforce monotonicity
# # 
# # trying this because we may not need to change the loess call at all, enforcing monotonicity may be good enough. 
# 
# 
# library(magrittr)
# library(dplyr)
# 
# loessErrfun_mod <- function(trans) {
#   qq <- as.numeric(colnames(trans))
#   est <- matrix(0, nrow=0, ncol=length(qq))
#   for(nti in c("A","C","G","T")) {
#     for(ntj in c("A","C","G","T")) {
#       if(nti != ntj) {
#         errs <- trans[paste0(nti,"2",ntj),]
#         tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
#         rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
#         rlogp[is.infinite(rlogp)] <- NA
#         df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
#         
#         # original
#         # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
#         mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
#         # #        mod.lo <- loess(rlogp ~ q, df)
#         
#         # Gulliem Salazar's solution 
#         # https://github.com/benjjneb/dada2/issues/938
#         # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
#         
#         pred <- predict(mod.lo, qq)
#         maxrli <- max(which(!is.na(pred)))
#         minrli <- min(which(!is.na(pred)))
#         pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
#         pred[seq_along(pred)<minrli] <- pred[[minrli]]
#         est <- rbind(est, 10^pred)
#       } # if(nti != ntj)
#     } # for(ntj in c("A","C","G","T"))
#   } # for(nti in c("A","C","G","T"))
#   
#   # HACKY
#   MAX_ERROR_RATE <- 0.25
#   MIN_ERROR_RATE <- 1e-7
#   est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
#   est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
#   
#   # enforce monotonicity 
#   # https://github.com/benjjneb/dada2/issues/791
#   estorig <- est
#   est <- est %>% 
#     data.frame() %>%
#     mutate_all(funs(case_when(. < X40 ~ X40,
#                               . >= X40 ~ .))) %>% as.matrix()
#   rownames(est) <- rownames(estorig)
#   colnames(est) <- colnames(estorig)
#   
#   # Expand the err matrix with the self-transition probs
#   err <- rbind(1-colSums(est[1:3,]), est[1:3,],
#                est[4,], 1-colSums(est[4:6,]), est[5:6,],
#                est[7:8,], 1-colSums(est[7:9,]), est[9,],
#                est[10:12,], 1-colSums(est[10:12,]))
#   rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
#   colnames(err) <- colnames(trans)
#   # Return
#   return(err)
# }
# 
# 
# # check what this looks like
# errF <- learnErrors(
#   filtFs, 
#   multithread = TRUE, 
#   errorEstimationFunction = loessErrfun_mod,
#   verbose = TRUE
# )
# dada2:::checkConvergence(errF)
# plotErrors(errF, nominalQ=TRUE)
# ggsave(
#   file.path(figs_path_16S,"learnErrors_errF_enforceMono.pdf"),
#   width = 7, height = 7, units = "in"
# )
# saveRDS(errF, file.path(rds_path_16S, "errF_enforceMono.RDS"))
# 
# # conclusions:
# # without altering the span settings in loess() we still have sharp peaks and dips. 
# # enforcing decreasing rates the way that we do does appear to clean up that issue quite nicely. 


# #----------------------------------------------------------
# ########
# # specify error model
# ########
# # now we'll define a new error model, building off the idea in the modified function above
# # 
# # Try 3
# #   loess() : implement changes to weights, don't change span
# #   enforce monotonicity
# 
# library(magrittr)
# library(dplyr)
# 
# loessErrfun_mod <- function(trans) {
#   qq <- as.numeric(colnames(trans))
#   est <- matrix(0, nrow=0, ncol=length(qq))
#   for(nti in c("A","C","G","T")) {
#     for(ntj in c("A","C","G","T")) {
#       if(nti != ntj) {
#         errs <- trans[paste0(nti,"2",ntj),]
#         tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
#         rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
#         rlogp[is.infinite(rlogp)] <- NA
#         df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
#         
#         # original
#         # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
#         # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
#         # #        mod.lo <- loess(rlogp ~ q, df)
#         
#         # Gulliem Salazar's solution 
#         # https://github.com/benjjneb/dada2/issues/938
#         # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
#         
#         # only change the weights
#         mod.lo <- loess(rlogp ~ q, df, weights = log10(tot))
#         
#         pred <- predict(mod.lo, qq)
#         maxrli <- max(which(!is.na(pred)))
#         minrli <- min(which(!is.na(pred)))
#         pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
#         pred[seq_along(pred)<minrli] <- pred[[minrli]]
#         est <- rbind(est, 10^pred)
#       } # if(nti != ntj)
#     } # for(ntj in c("A","C","G","T"))
#   } # for(nti in c("A","C","G","T"))
#   
#   # HACKY
#   MAX_ERROR_RATE <- 0.25
#   MIN_ERROR_RATE <- 1e-7
#   est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
#   est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
#   
#   # enforce monotonicity 
#   # https://github.com/benjjneb/dada2/issues/791
#   estorig <- est
#   est <- est %>% 
#     data.frame() %>%
#     mutate_all(funs(case_when(. < X40 ~ X40,
#                               . >= X40 ~ .))) %>% as.matrix()
#   rownames(est) <- rownames(estorig)
#   colnames(est) <- colnames(estorig)
#   
#   # Expand the err matrix with the self-transition probs
#   err <- rbind(1-colSums(est[1:3,]), est[1:3,],
#                est[4,], 1-colSums(est[4:6,]), est[5:6,],
#                est[7:8,], 1-colSums(est[7:9,]), est[9,],
#                est[10:12,], 1-colSums(est[10:12,]))
#   rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
#   colnames(err) <- colnames(trans)
#   # Return
#   return(err)
# }
# 
# # check what this looks like
# errF <- learnErrors(
#   filtFs, 
#   multithread = TRUE, 
#   errorEstimationFunction = loessErrfun_mod,
#   verbose = TRUE
# )
# dada2:::checkConvergence(errF)
# plotErrors(errF, nominalQ=TRUE)
# ggsave(
#   file.path(figs_path_16S,"learnErrors_errF_weights_enforceMono.pdf"),
#   width = 7, height = 7, units = "in"
# )
# saveRDS(errF, file.path(rds_path_16S, "errF_weights_enforceMono.RDS"))
# 
# # conclusions:
# # without adjusting span, we sill have dips and peaks, although they are much less distinct. This is most likely due to the effects of altering the argument passed to weight. 
# # enforcing decreasing error rates did it's job as expected. 


#----------------------------------------------------------
########
# error function decision
########
# We will be using the error function that 
# Try 1
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
# learn errors
########
# errF <- learnErrors(
#   filtFs, 
#   multithread = TRUE, 
#   errorEstimationFunction = loessErrfun_mod,
#   verbose = TRUE
# )
# dada2:::checkConvergence(errF)
# plotErrors(errF, nominalQ=TRUE)
# ggsave(
#   file.path(figs_path_16S,"learnErrors_errF.pdf"),
#   width = 7, height = 7, units = "in"
# )
# saveRDS(errF, file.path(rds_path_16S, "errF.RDS"))

#----------------------------------------------------------
########
# dada
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
saveRDS(dadaFs, file.path(rds_path_16S, "dadaFs.RDS"))
ggsave(
  plotErrors(dadaFs, nominalQ=TRUE),
  file.path(figs_path_16S,"dada_errF.pdf"),
  width = 7, height = 7, units = "in"
)

#----------------------------------------------------------
########
# are the error objects identical? 
########
# identical(dada2::getErrors(errF), dada2::getErrors(dadaFs))

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########