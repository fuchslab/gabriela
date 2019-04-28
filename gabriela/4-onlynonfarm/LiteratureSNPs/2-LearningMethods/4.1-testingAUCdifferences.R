# MARKNEW 4 lines
setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/4-onlynonfarm/LiteratureSNPs/")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/5-prediction-and-validation.R")
load("../../Data/dat.outcome.rda")
load("../../Data/gabriel.1707.Rdata")
# MARKNEW
load("../../Data/ind.part.rda")
# do for all  for imputations
# MARKNEW
load("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/LiteratureSNPs/2-LearningMethods/Results19/aucs/aucs_template.rda")
aucs.b.i[1:nrow(aucs.b.i), 2:ncol(aucs.b.i)] <- NA

## MARKNEW add glm-weights in template already here
aucs.b.i <- rbind(aucs.b.i, aucs.b.i[1, ])
aucs.b.i$method[nrow(aucs.b.i)] <- "glm_weights"
colnames(aucs.b.i)[2:7] <- c("snps_only", "snps_only.lo", "snps_only.up", "env_only", "env_only.lo","env_only.up" )

## MARKNEW add conf_only
aucs.b.i[1:nrow(aucs.b.i), 2:ncol(aucs.b.i)] <- NA
aucs.b.i <- data.frame(method=aucs.b.i$method, conf_only=NA, conf_only.lo=NA, conf_only.up=NA, aucs.b.i[ ,-1])


all.aucs <- list()
all.aucs$dda <- aucs.b.i
all.aucs$dda_notany <- aucs.b.i
all.aucs$dda.wh12_notany <- aucs.b.i
all.aucs$asthma_wo_bronch <- aucs.b.i 
# colnames(dat.outcome)[colnames(dat.outcome)=="dda.wh12_notany"] <- "wh12_notany"

# boostrap samples
B=10000
cores=20
set.seed(13412342)
for(outcome in c( "dda")){
  
  # MARKNEW: next 3 lines
  ind <- ind.part$cvfolds.without.inn$nonfarm
  y <- lapply(ind, function(x) gabriel.1707$dd_asthma[setdiff(ind.part$nonfarm, ind.part$nonfarm.inn) ][x])
  w <- lapply(ind, function(x) gabriel.1707$weight_dna[setdiff(ind.part$nonfarm, ind.part$nonfarm.inn) ][x])
  
  for(mod in c("conf_only" ,"snps_only", "env_only","fh_only", "interactions_only",  "both")){
    # for(alg in c("ranger", "elnet", "elnet_weights", "glmnet", "glmnet_weights", "ipflasso",  "ipflasso_weights")){
      for(alg in c("ranger", "elnet_weights", "glmnet_weights",  "ipflasso_weights")){
      if( ( alg=="ipflasso_weights" | alg=="ipflasso") & mod != "both"){
        break
      } else{
        aucs <- list()
        for(i in 1:5){
          load(paste0("2-LearningMethods/Results19/imp",i,"/predictions/pred.",alg,"_",outcome,"_estdos_",mod,".rda"))
          aucs[[i]] <- my.bs.2samples.roc(yhat1=unlist(pred),yhat2=NULL, y=unlist(y), prob=unlist(w), B=B, cores = cores)# NOTE: it is not 1/weights !!! (these would be inverse-inverse-probality weights which is wrong)
          lo <- mean(unlist(lapply(aucs, function(x) x[1] )))
          up <- mean(unlist(lapply(aucs, function(x) x[3] )))
          med <- mean(unlist(lapply(aucs, function(x) x[2] )))
          all.aucs[[outcome]][aucs.b.i$method==alg, which(colnames(aucs.b.i) == mod):(which(colnames(aucs.b.i) == mod) +2) ] <- c(med, lo, up)
        }
      }
    }
  }
  print(paste("methods done - tests starting..."))
  ## do the testing:
  method.both.max <- aucs.b.i$method[which(all.aucs[[outcome]]$both == max(all.aucs[[outcome]]$both, na.rm = T)) ]
  method.snps_only.max <- aucs.b.i$method[which(all.aucs[[outcome]]$snps_only == max(all.aucs[[outcome]]$snps_only, na.rm = T)) ]
  method.env_only.max <- aucs.b.i$method[which(all.aucs[[outcome]]$env_only == max(all.aucs[[outcome]]$env_only, na.rm = T)) ]
  method.fh_only.max <- aucs.b.i$method[which(all.aucs[[outcome]]$fh_only == max(all.aucs[[outcome]]$fh_only, na.rm = T)) ]
  
  for(i in 1:5){
    load(paste0("2-LearningMethods/Results19/imp",i,"/predictions/pred.", method.both.max,"_", outcome, "_estdos_both.rda"))
    pred.both <- pred
    load(paste0("2-LearningMethods/Results19/imp", i, "/predictions/pred.", method.fh_only.max, "_", outcome, "_estdos_fh_only.rda"))
    pred.fh_only <- pred
    load(paste0("2-LearningMethods/Results19/imp", i, "/predictions/pred.", method.snps_only.max, "_", outcome, "_estdos_snps_only.rda"))
    pred.snps_only <- pred
    load(paste0("2-LearningMethods/Results19/imp", i, "/predictions/pred.", method.env_only.max, "_", outcome, "_estdos_env_only.rda"))
    pred.env_only <- pred
    
    # NOTE: it is not 1/weights !!! (these would be inverse-inverse-probality weights which is wrong)
    all.aucs$tests[[outcome]]$both.snps_only[[i]] <-   my.bs.2samples.roc(unlist(pred.snps_only), unlist(pred.both), y=unlist(y), prob=unlist(w), B=B, alpha.2samples = .05/3, cores = cores)
    all.aucs$tests[[outcome]]$both.env_only[[i]] <-   my.bs.2samples.roc(unlist(pred.env_only), unlist(pred.both), y=unlist(y), prob=unlist(w), B=B, alpha.2samples = .05/3, cores = cores)
    all.aucs$tests[[outcome]]$both.fh_only[[i]] <-  my.bs.2samples.roc(unlist(pred.fh_only), unlist(pred.both), y=unlist(y), prob=unlist(w), B=B, alpha.2samples = .05/3, cores = cores)
  }
  print(paste("outcome", outcome, "done."))
}
save(all.aucs, file=paste0("2-LearningMethods/Results19/aucs/all.aucs.with.tests.",B ,"bootstrap_samples.rda"))
# performances
all.aucs$dda[-c(8,9), ]
# CI dda both vs snps_only
c(mean(unlist(lapply(all.aucs$tests$dda$both.snps_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda$both.snps_only, function(x) x$diffaucs[3]))))
# CI dda both vs env_only
c(mean(unlist(lapply(all.aucs$tests$dda$both.env_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda$both.env_only, function(x) x$diffaucs[3]))))
# CI dda both vs fh_only
c(mean(unlist(lapply(all.aucs$tests$dda$both.fh_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda$both.fh_only, function(x) x$diffaucs[3]))))


# performances
all.aucs$dda_notany[-c(8,9), ]
# CI dda_notany both vs snps_only
c(mean(unlist(lapply(all.aucs$tests$dda_notany$both.snps_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda_notany$both.snps_only, function(x) x$diffaucs[3]))))
# CI dda_notany both vs env_only
c(mean(unlist(lapply(all.aucs$tests$dda_notany$both.env_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda_notany$both.env_only, function(x) x$diffaucs[3]))))
# CI dda_notany both vs fh_only
c(mean(unlist(lapply(all.aucs$tests$dda_notany$both.fh_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda_notany$both.fh_only, function(x) x$diffaucs[3]))))

# performances
all.aucs$dda.wh12_notany[-c(8,9), ]
# CI dda.wh12_notany both vs snps_only
c(mean(unlist(lapply(all.aucs$tests$dda.wh12_notany$both.snps_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda.wh12_notany$both.snps_only, function(x) x$diffaucs[3]))))
# CI dda.wh12_notany both vs env_only
c(mean(unlist(lapply(all.aucs$tests$dda.wh12_notany$both.env_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda.wh12_notany$both.env_only, function(x) x$diffaucs[3]))))
# CI dda.wh12_notany both vs fh_only
c(mean(unlist(lapply(all.aucs$tests$dda.wh12_notany$both.fh_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$dda.wh12_notany$both.fh_only, function(x) x$diffaucs[3]))))


# performances
all.aucs$asthma_wo_bronch[-c(8,9), ]
# CI asthma_wo_bronch both vs snps_only
c(mean(unlist(lapply(all.aucs$tests$asthma_wo_bronch$both.snps_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$asthma_wo_bronch$both.snps_only, function(x) x$diffaucs[3]))))
# CI asthma_wo_bronch both vs env_only
c(mean(unlist(lapply(all.aucs$tests$asthma_wo_bronch$both.env_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$asthma_wo_bronch$both.env_only, function(x) x$diffaucs[3]))))
# CI asthma_wo_bronch both vs fh_only
c(mean(unlist(lapply(all.aucs$tests$asthma_wo_bronch$both.fh_only, function(x) x$diffaucs[1]))), mean(unlist(lapply(all.aucs$tests$asthma_wo_bronch$both.fh_only, function(x) x$diffaucs[3]))))







### specific tests:
# dda - interactions - elnet weighted vs. fh best alg. (weighted)
aucs.test <- list()
for(i in 1:5){
  load(paste0("2-LearningMethods/Results19/imp",i,"/predictions/pred.ranger_dda_estdos_both.rda"))
  pred.fh_only <- pred
  load(paste0("2-LearningMethods/Results19/imp",i,"/predictions/pred.elnet_weights_dda_estdos_interactions_only.rda"))
  pred.interactions <- pred
  # NOTE: it is not 1/weights !!! (these would be inverse-inverse-probality weights which is wrong)
  aucs.test[[i]] <-   my.bs.2samples.roc(unlist(pred.interactions), unlist(pred.fh_only), y=unlist(y), prob=unlist(w), B=B, alpha.2samples = .05/3, cores = cores)
}
c(mean(unlist(lapply(aucs.test, function(x) x$diffaucs[1]))), mean(unlist(lapply(aucs.test, function(x) x$diffaucs[3]))))


# 
# bla <- list()
# for(i in 1:5){
#   load(paste0("2-LearningMethods/Results19/imp",i,"/predictions/pred.",alg,"_",outcome,"_estdos_",mod,".rda"))
#  rocobject <- roc(unlist(y), unlist(pred))
#   bla[[i]] <- ci.auc.bootstrap.w(rocobject, conf.level = .95, boot.n = B, weights = unlist(w.test), parallel = F)
# }
# c(mean(unlist(lapply(bla, function(x) x[1]))), 
#   mean(unlist(lapply(bla, function(x) x[2]))), 
#   mean(unlist(lapply(bla, function(x) x[3]))))
#### bootstrap test for comparing two samples...


set.seed(234)
y <- rnorm(200)
x <- y + rnorm(100, mean=.1, sd=1)

##
t.test(x,y, paired = T)

## package boot
# library(boot)
# dat <- data.frame(y,x)
# m <- function(d,x) mean(d$x - x)
# b1 <- boot(dat, statistic=m, R=1000)
# boot.ci(b1, conf=.95)
# ?boot.ci
# 
# ratio <- function(d, w) mean(d$x * w)/sum(d$u * w)
# boot(dat, ratio, R = 999)


#### 
my.bs.2samples <- function(x,y,B, prob=NULL){
  ind <- lapply(1:B, function(z) sample(1:length(x), size=length(x), replace=T, prob=prob))
  diffs <- unlist(lapply(ind, function(z) mean(x[z]) - mean(y[z]) ))
  quantile(diffs, probs = c(.025, 0.5, .95))
}

my.bs.2samples(x,y,B=10000, prob=NULL)
t.test(x,y, paired=T)

###### now on our AUC differences

## modify function
library(pROC)


##### version 2
# library(pROC)
# 
# my.bs.2samples.roc2 <- function(yhat1,yhat2,y,B, prob=NULL){
#   roc(y[z], yhat1[z], direction = "<")$auc
#   auc1 <- ci.auc.bootstrap.w(roc=, boot.n=1, parallel = FALSE, conf.level=.95, weights=prob)
#   
#   ind <- lapply(1:B, function(z) sample(1:length(y), size=length(y), replace=T, prob=prob))
#   aucs <- lapply(ind, function(z){
#     auc1 <- roc(y[z], yhat1[z], direction = "<")$auc
#     auc2 <- roc(y[z], yhat2[z], direction = "<")$auc 
#     # print(auc1, auc2)
#     c(auc1, auc2)
#   }                            )
#   diffs <- unlist(lapply(aucs, function(x) x[2] - x[1] ))
#   results <- list()
#   
#   results$aucs1 <- quantile(unlist(lapply(aucs, function(x) x[1])), probs=c(.025, .5, .975))
#   results$aucs2 <- quantile(unlist(lapply(aucs, function(x) x[2])), probs=c(.025, .5, .975))
#   results$diffaucs <- quantile(diffs, probs = c(.025, 0.5, 1-0.025))
#   results
# }

source("~/1-Projects/1-R-functions/5-prediction-and-validation.R")


setwd("~/1-Projects/2-Asthma/LiteratureSNPs/")
# do for all  for imputations
diffs.imps <- list()
for(i in 1:5){
  load(paste0("2-LearningMethods/Results19/imp",i,"/predictions/pred.ranger_dda_estdos_fh_only.rda"))
  pred.fh <- pred
  load(paste0("2-LearningMethods/Results19/imp", i, "/predictions/pred.ipflasso_weights_dda_estdos_both.rda"))
  pred.both <- pred
  load("../Data/fin.ind.rda")
  load("../Data/gabriel.1707.Rdata")
  ind <- fin.ind$ind.within$ind.test
  y <- lapply(ind, function(x) gabriel.1707$dd_asthma[fin.ind$train.fin.ind ][x])
  w <- lapply(ind, function(x) gabriel.1707$weight_dna[fin.ind$train.fin.ind ][x])
  diffs.imps[[i]] <- my.bs.2samples.roc(unlist(pred.fh), unlist(pred.both), y=unlist(y), prob=unlist(w), B=10)
}
diffs.imps

# means seperately:
mean(unlist(lapply(diffs.imps, function(x) x$aucs1[2])))
mean(unlist(lapply(diffs.imps, function(x) x$aucs2[2])))
### all overlap with 0 --> unambiguously not significant


############################################################
#### now same for dda vs not any
############################################################
load("../Data/fin.ind.rda")
# load("../Data/gabriel.1707.Rdata")
load("../Data/dat.outcome.rda")


y.with.nas <- dat.outcome$dda_notany[fin.ind$train.fin.ind]
y.unlist <- na.omit(y.with.nas)
w.with.nas <- gabriel.1707$weight_dna[fin.ind$train.fin.ind ]
w.unlist <- w.with.nas[!is.na(y.with.nas)]
#get right innner indices back 
######### CV (or general arbitrary partitioning)  
  k <- 5
  n <- length(y.unlist)
  set.seed(3534)
  ind <- sample(1:n)
  
  ind.test <- list()
  ind.train <- list()
  for(i in 1:k){
    if(i < k ) ind.test[[i]] <- ind[(1+(i-1)*floor(n/k)):(i*(floor(n/k))) ]
    else ind.test[[i]] <- ind[(1+(i-1)*floor(n/k)):n]
    ind.train[[i]] <- setdiff(ind, ind.test[[i]])
  }
  ind <- ind.test 
  y <- lapply(ind, function(x) y.unlist[x])
  w <- lapply(ind, function(x) w.unlist[x])

setwd("~/1-Projects/2-Asthma/LiteratureSNPs/")
# do for all  for imputations
diffs.imps <- list()
for(i in 1:5){
  load(paste0("2-LearningMethods/Results19/imp",i,"/predictions/pred.ranger_dda_notany_estdos_fh_only.rda"))
  pred.fh <- pred
  load(paste0("2-LearningMethods/Results2ipflong/imp", i, "/predictions/pred.ipflasso_weights_dda_notany_estdos_both.rda"))
  pred.both <- pred
  diffs.imps[[i]] <- my.bs.2samples.roc(unlist(pred.fh), unlist(pred.both), y=unlist(y), prob=unlist(w), B=100, alpha.2samples = .05/3)
}
diffs.imps
mean(unlist(lapply(diffs.imps, function(x) x$aucs1[2])))
mean(unlist(lapply(diffs.imps, function(x) x$aucs2[2])))
str(diffs.imps)


mean(unlist(lapply(diffs.imps, function(x) x$aucs1[2])))

bla <- roc(unlist(y), unlist(pred.both))





# ### check for no significant difference: compare pred.elnet with pred.ranger by fh
# diffs.imps <- list()
# for(i in 1:5){
#   load(paste0("/Users/norbertkrautenbacher/1-Projects/2-Asthma/LiteratureSNPs/2-LearningMethods/Results19/imp",i,"/predictions/pred.ranger_dda_estdos_fh_only.rda"))
#   pred.fh <- pred
#   load(paste0("/Users/norbertkrautenbacher/1-Projects/2-Asthma/LiteratureSNPs/2-LearningMethods/Results19/imp", i, "/predictions/pred.elnet_dda_estdos_fh_only.rda"))
#   pred.both <- pred
#   load("/Users/norbertkrautenbacher/1-Projects/2-Asthma/Data/fin.ind.rda")
#   load("/Users/norbertkrautenbacher/1-Projects/2-Asthma/Data/gabriel.1707.Rdata")
#   ind <- fin.ind$ind.within$ind.test
#   y <- lapply(ind, function(x) gabriel.1707$dd_asthma[fin.ind$train.fin.ind ][x])
#   w <- lapply(ind, function(x) gabriel.1707$weight_dna[fin.ind$train.fin.ind ][x])
#   diffs.imps[[i]] <- my.bs.2samples.roc(yhat1=unlist(pred.fh),yhat2= unlist(pred.both), y=unlist(y), prob=1/unlist(w), B=100)
# }
# diffs.imps
# 
# roc(unlist(y), unlist(pred.fh))
# roc(unlist(y), unlist(pred.both))

######## same for dda vs not_any


