# MARKNEW next 4 lines
setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/2-validate-innsbruck/ImputedSNPs/")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/5-prediction-and-validation.R")
load("../../Data/dat.outcome.rda")
load("../../Data/gabriel.1707.Rdata")
# MARKNEW
load("../../Data/ind.part.rda")
# do for all  for imputations
# MARKNEW
load("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/ImputedSNPs/Results_genomewide/aucs/aucs_dda_genomewide.rda")
aucs.b.i[1:nrow(aucs.b.i), 2:ncol(aucs.b.i)] <- NA
colnames(aucs.b.i)[2:7] <- c("snps_only", "snps_only.lo", "snps_only.up", "env_only", "env_only.lo","env_only.up" )
all.aucs <- list()
all.aucs$dda <- aucs.b.i
# colnames(dat.outcome)[colnames(dat.outcome)=="dda.wh12_notany"] <- "wh12_notany"

# boostrap samples
B=10000
cores=20
set.seed(13412342)
for(outcome in c("dda")){
  # MARKNEW: next 3 lines
    ind <- ind.part$cvfolds.without.inn$all
    y <- lapply(ind, function(x) gabriel.1707$dd_asthma[setdiff(1:1707, ind.part$all.inn) ][x])
    w <- lapply(ind, function(x) gabriel.1707$weight_dna[setdiff(1:1707, ind.part$all.inn) ][x])
    # w.test <- gabriel.1707$weight_dna[fin.ind$train.fin.ind ]
  # }
  for(mod in c( "snps_only")){
    for(alg in c("ranger", "glm_gwas", "elnet_weights", "glmnet_weights")){
      if( ( alg=="ipflasso_weights" | alg=="ipflasso"| alg=="glm_mv_gwas") & mod != "both"){
        break
      } else{
        aucs <- list()
        for(i in 1:5){
          load(paste0("Results_genomewide/imp",i,"/predictions/pred.",alg,"_",outcome,"_genomewide_",mod,".rda"))
          aucs[[i]] <- my.bs.2samples.roc(yhat1=unlist(pred),yhat2=NULL, y=unlist(y), prob=unlist(w), B=B, cores = cores)# NOTE: it is not 1/weights !!! (these would be inverse-inverse-probality weights which is wrong)
          lo <- mean(unlist(lapply(aucs, function(x) x[1] )))
          up <- mean(unlist(lapply(aucs, function(x) x[3] )))
          med <- mean(unlist(lapply(aucs, function(x) x[2] )))
          all.aucs[[outcome]][aucs.b.i$method==alg, which(colnames(aucs.b.i) == mod):(which(colnames(aucs.b.i) == mod) +2) ] <- c(med, lo, up)
        }
      }
    }
  }
  print(paste("methods done."))

}
save(all.aucs, file=paste0("Results_genomewide/aucs/all.aucs.with.tests.",B ,"bootstrap_samples.rda"))
# performances
all.aucs$dda[-c(8,9), ]
# # CI dda both vs snps_only
# c(mean(unlist(lapply(all.aucs$tests$dda$both.snps_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda$both.snps_only, function(x) x$diffaucs[3]))))
# # CI dda both vs env_only
# c(mean(unlist(lapply(all.aucs$tests$dda$both.env_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda$both.env_only, function(x) x$diffaucs[3]))))
# # CI dda both vs fh_only
# c(mean(unlist(lapply(all.aucs$tests$dda$both.fh_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda$both.fh_only, function(x) x$diffaucs[3]))))
# 
# 
# # performances
# all.aucs$dda_notany[-c(8,9), ]
# # CI dda_notany both vs snps_only
# c(mean(unlist(lapply(all.aucs$tests$dda_notany$both.snps_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda_notany$both.snps_only, function(x) x$diffaucs[3]))))
# # CI dda_notany both vs env_only
# c(mean(unlist(lapply(all.aucs$tests$dda_notany$both.env_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda_notany$both.env_only, function(x) x$diffaucs[3]))))
# # CI dda_notany both vs fh_only
# c(mean(unlist(lapply(all.aucs$tests$dda_notany$both.fh_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda_notany$both.fh_only, function(x) x$diffaucs[3]))))
# 
# # performances
# all.aucs$dda.wh12_notany[-c(8,9), ]
# # CI dda.wh12_notany both vs snps_only
# c(mean(unlist(lapply(all.aucs$tests$dda.wh12_notany$both.snps_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda.wh12_notany$both.snps_only, function(x) x$diffaucs[3]))))
# # CI dda.wh12_notany both vs env_only
# c(mean(unlist(lapply(all.aucs$tests$dda.wh12_notany$both.env_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda.wh12_notany$both.env_only, function(x) x$diffaucs[3]))))
# # CI dda.wh12_notany both vs fh_only
# c(mean(unlist(lapply(all.aucs$tests$dda.wh12_notany$both.fh_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda.wh12_notany$both.fh_only, function(x) x$diffaucs[3]))))
# 
# ### specific tests:
# dda - interactions - elnet weighted vs. fh best alg. (weighted)
# aucs.test <- list()
# for(i in 1:5){
#   load(paste0("2-LearningMethods/Results19/imp",i,"/predictions/pred.ranger_dda_estdos_both.rda"))
#   pred.fh_only <- pred
#   load(paste0("2-LearningMethods/Results19/imp",i,"/predictions/pred.elnet_weights_dda_estdos_interactions_only.rda"))
#   pred.interactions <- pred
#   # NOTE: it is not 1/weights !!! (these would be inverse-inverse-probality weights which is wrong)
#   aucs.test[[i]] <-   my.bs.2samples.roc(unlist(pred.interactions), unlist(pred.fh_only), y=unlist(y), prob=unlist(w), B=B, alpha.2samples = .05/3, cores = cores)
# }
# c(mean(unlist(lapply(aucs.test, function(x) x$diffaucs[1]))), mean(unlist(lapply(aucs.test, function(x) x$diffaucs[3]))))
# 

# 
