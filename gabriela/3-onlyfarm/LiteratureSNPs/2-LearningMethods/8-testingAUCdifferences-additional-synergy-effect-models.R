# MARKNEW 4 lines

setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/3-onlyfarm/LiteratureSNPs/")
best.syn.model <- "fh_conf_snps"


source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/5-prediction-and-validation.R")
load("../../Data/dat.outcome.rda")
load("../../Data/gabriel.1707.Rdata")
# MARKNEW
load("../../Data/ind.part.rda")
# do for all  for imputations



cores=1
B=10000
load(file=paste0("2-LearningMethods/Results19/aucs/all.aucs_for_synergy_models_with.tests.",B ,"bootstrap_samples.rda"))



set.seed(13412342)
outcome = "dda"

# MARKNEW: next 3 lines
ind <- ind.part$cvfolds.without.inn$farm
y <- lapply(ind, function(x) gabriel.1707$dd_asthma[setdiff(ind.part$farm, ind.part$farm.inn) ][x])
w <- lapply(ind, function(x) gabriel.1707$weight_dna[setdiff(ind.part$farm, ind.part$farm.inn) ][x])


method.syn.best.max <- all.aucs.syn[[outcome]]$method[which(all.aucs.syn[[outcome]][[best.syn.model]] == max(all.aucs.syn[[outcome]][[best.syn.model]], na.rm = T)) ]
method.fh_only.max <- all.aucs.syn[[outcome]]$method[which(all.aucs.syn[[outcome]][["fh_only"]] == max(all.aucs.syn[[outcome]][["fh_only"]], na.rm = T)) ]


tests <- list()
for(i in 1:5){
  load(paste0("2-LearningMethods/Results19/imp",i,"/predictions/pred.", method.syn.best.max,"_", outcome, "_estdos_", best.syn.model , ".rda"))
  pred.syn.best <- pred
  load(paste0("2-LearningMethods/Results19/imp", i, "/predictions/pred.", method.fh_only.max, "_", outcome, "_estdos_fh_only.rda"))
  pred.fh_only <- pred
  
  ### caclulate Bayesfactor instead
  my.seed <- 123
  set.seed(my.seed)
  syn.boot = boot_AUC( unlist(y), unlist(pred.syn.best), nboot=B, weights=unlist(w), cores=cores)
  set.seed(my.seed)
  fh.boot = boot_AUC(unlist(y), unlist(pred.fh_only), nboot=B, weights=unlist(w), cores=cores)
  
  # save allAUCs
  tests$syn.boot[[i]] <- syn.boot
  tests$fh.boot[[i]] <- fh.boot
  
  # Bayes factors
  tests$bayesfactor[[i]] <- (sum(syn.boot$AUC > fh.boot$AUC)/B)/(sum(syn.boot$AUC <= fh.boot$AUC)/B)
  
  # p-value
  tests$pval[[i]] <- 1-(sum(syn.boot$AUC > fh.boot$AUC)/B)
  
  # percentile CIs
  tests$percCI[[i]] <- quantile((syn.boot$AUC - fh.boot$AUC), probs = c(.025, .975))
  tests
}
save(tests, file=paste0("2-LearningMethods/Results19/aucs/aucs_difference_test_synergy_vs_standalone.",B ,"bootstrap_samples.rda") )  

