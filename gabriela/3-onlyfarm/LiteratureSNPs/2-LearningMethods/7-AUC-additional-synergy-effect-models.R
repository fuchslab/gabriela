# MARKNEW 4 lines
setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/3-onlyfarm/LiteratureSNPs/")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/5-prediction-and-validation.R")
load("../../Data/dat.outcome.rda")
load("../../Data/gabriel.1707.Rdata")
# MARKNEW
load("../../Data/ind.part.rda")
# do for all  for imputations


env.ifs <- c("fh_only", "fh_conf", "fh_conf_env", "fh_conf_snps", "fh_conf_env_snps",  "both")
methods <- c("glmnet_weights", "ranger", "ipflasso_weights")


names.df <- unlist(lapply(env.ifs, function(x) paste0(x, c("", ".lo", ".up"))))
names.df

df <- data.frame(matrix(NA, 3,length(names.df)))
colnames(df) <-  names.df
df$method <- methods




all.aucs <- list()
all.aucs$dda <- df
# colnames(dat.outcome)[colnames(dat.outcome)=="dda.wh12_notany"] <- "wh12_notany"

# boostrap samples
B=10000
cores=50
set.seed(13412342)
for(outcome in c( "dda")){
  
  # MARKNEW: next 3 lines
  ind <- ind.part$cvfolds.without.inn$farm
  y <- lapply(ind, function(x) gabriel.1707$dd_asthma[setdiff(ind.part$farm, ind.part$farm.inn) ][x])
  w <- lapply(ind, function(x) gabriel.1707$weight_dna[setdiff(ind.part$farm, ind.part$farm.inn) ][x])
  
  for(mod in c( "fh_conf", "fh_conf_env", "fh_conf_snps", "fh_conf_env_snps")){
    # for(alg in c("ranger", "elnet", "elnet_weights", "glmnet", "glmnet_weights", "ipflasso",  "ipflasso_weights")){
    for(alg in methods){
      if( ( alg=="ipflasso_weights" | alg=="ipflasso") & mod == "fh"){
        break
      } else{
        aucs <- list()
        for(i in 1:5){
          load(paste0("2-LearningMethods/Results19/imp",i,"/predictions/pred.",alg,"_",outcome,"_estdos_",mod,".rda"))
          aucs[[i]] <- my.bs.2samples.roc(yhat1=unlist(pred),yhat2=NULL, y=unlist(y), prob=unlist(w), B=B, cores = cores)# NOTE: it is not 1/weights !!! (these would be inverse-inverse-probality weights which is wrong)
          lo <- mean(unlist(lapply(aucs, function(x) x[1] )))
          up <- mean(unlist(lapply(aucs, function(x) x[3] )))
          med <- mean(unlist(lapply(aucs, function(x) x[2] )))
          all.aucs[[outcome]][df$method==alg, which(colnames(df) == mod):(which(colnames(df) == mod) +2) ] <- c(med, lo, up)
        }
      }
    }
  }
  print(paste("methods done - tests starting..."))
  ## do the testing:
  # method.both.max <- df$method[which(all.aucs[[outcome]]$both == max(all.aucs[[outcome]]$both, na.rm = T)) ]
  # method.snps_only.max <- df$method[which(all.aucs[[outcome]]$snps_only == max(all.aucs[[outcome]]$snps_only, na.rm = T)) ]
  # method.env_only.max <- df$method[which(all.aucs[[outcome]]$env_only == max(all.aucs[[outcome]]$env_only, na.rm = T)) ]
  # method.fh_only.max <- df$method[which(all.aucs[[outcome]]$fh_only == max(all.aucs[[outcome]]$fh_only, na.rm = T)) ]
  # 
  # for(i in 1:5){
  #   load(paste0("2-LearningMethods/Results19/imp",i,"/predictions/pred.", method.both.max,"_", outcome, "_estdos_both.rda"))
  #   pred.both <- pred
  #   load(paste0("2-LearningMethods/Results19/imp", i, "/predictions/pred.", method.fh_only.max, "_", outcome, "_estdos_fh_only.rda"))
  #   pred.fh_only <- pred
  #   load(paste0("2-LearningMethods/Results19/imp", i, "/predictions/pred.", method.snps_only.max, "_", outcome, "_estdos_snps_only.rda"))
  #   pred.snps_only <- pred
  #   load(paste0("2-LearningMethods/Results19/imp", i, "/predictions/pred.", method.env_only.max, "_", outcome, "_estdos_env_only.rda"))
  #   pred.env_only <- pred
  #   
  #   # NOTE: it is not 1/weights !!! (these would be inverse-inverse-probality weights which is wrong)
  #   all.aucs$tests[[outcome]]$both.snps_only[[i]] <-   my.bs.2samples.roc(unlist(pred.snps_only), unlist(pred.both), y=unlist(y), prob=unlist(w), B=B, alpha.2samples = .05/3, cores = cores)
  #   all.aucs$tests[[outcome]]$both.env_only[[i]] <-   my.bs.2samples.roc(unlist(pred.env_only), unlist(pred.both), y=unlist(y), prob=unlist(w), B=B, alpha.2samples = .05/3, cores = cores)
  #   all.aucs$tests[[outcome]]$both.fh_only[[i]] <-  my.bs.2samples.roc(unlist(pred.fh_only), unlist(pred.both), y=unlist(y), prob=unlist(w), B=B, alpha.2samples = .05/3, cores = cores)
  # }
  print(paste("outcome", outcome, "done."))
}

all.aucs.temp <- all.aucs
rm(all.aucs)
# performances

## add fh and both
load("2-LearningMethods/Results19/aucs/all.aucs.with.tests.10000bootstrap_samples.rda")
to.be.replaced <- c("fh_only", "fh_only.lo", "fh_only.up", "both", "both.lo", "both.up")
all.aucs.temp$dda[ , to.be.replaced] <- subset(all.aucs$dda,  select=to.be.replaced)[match(methods, all.aucs$dda$method), ]

rm(all.aucs)

all.aucs.syn <- all.aucs.temp

# for farm we had good prediction for snps genomewide on standalone basis, thus we add the results for the synergy effect in addition

load(file=paste0("../ImputedSNPs/Results_genomewide/aucs/all.aucs_for_synergy_models_with.tests.",B ,"bootstrap_samples_genomewide.rda"))
(all.aucs.syn.plus.genomewide <- base::merge(all.aucs.syn$dda, all.aucs.syn.genomewide$dda, by="method", all.x=T))


## lasso itself only interesting for fh_only --> make one line for ipf/lasso

all.aucs.syn.present <- all.aucs.syn.plus.genomewide[ all.aucs.syn.plus.genomewide$method %in% c("ranger", "ipflasso_weights") , ]
all.aucs.syn.present[ all.aucs.syn.present$method=="ipflasso_weights" , c("fh_only", "fh_only.lo", "fh_only.up")] <- all.aucs.syn$dda[ all.aucs.syn$dda$method=="glmnet_weights" , c("fh_only", "fh_only.lo", "fh_only.up")]
all.aucs.syn.present$method[all.aucs.syn.present$method == "ipflasso_weights"] <- "Log.Reg. (IPF-LASSO)"
all.aucs.syn.present$method[all.aucs.syn.present$method == "ranger"] <- "Random Forest"
all.aucs.syn.present


save(all.aucs.syn, all.aucs.syn.plus.genomewide, all.aucs.syn.present, file=paste0("2-LearningMethods/Results19/aucs/all.aucs_for_synergy_models_with.tests.",B ,"bootstrap_samples.rda"))


###  visualize
B=10000
# setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/3-onlyfarm/LiteratureSNPs/")
load(paste0("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/3-onlyfarm/LiteratureSNPs/2-LearningMethods/Results19/aucs/all.aucs_for_synergy_models_with.tests.",B ,"bootstrap_samples.rda"))
env.ifs <- c("fh_only", "fh_conf", "fh_conf_env", "fh_conf_snps", "fh_conf_env_snps",  "both")
## check it out
all.aucs.syn$dda[ , c("method", env.ifs)]

# pdf("onlyfarmsynergytemp.pdf")
barplot(as.matrix(all.aucs.syn$dda[ , env.ifs]), beside = T, ylim=c(0.5, 0.75), xpd=FALSE, main = "onlyfarm", legend.text = all.aucs.syn$dda$method)
# dev.off()