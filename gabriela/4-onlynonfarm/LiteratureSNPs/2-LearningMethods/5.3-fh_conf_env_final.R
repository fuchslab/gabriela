# NEWMARK
### lit snps
# NEWMARK
### lit snps
setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/4-onlynonfarm/LiteratureSNPs/2-LearningMethods/")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/3-wrp-learn-on-snp-env.R")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/0-functions.R")

library("ranger")

# NEWMARK

### Best model here was ranger


# NEWMARK
### data
load("../../../LiteratureSNPs/Data/dat.lit.snps.19.RData")
load("../../../Data/gabriel.1707.Rdata")
load("../../../Data/dat.outcome.rda")
# load("../../Data/x.env.imp.rda")
load("../../../Data/x.wheeze_onset.rda")

# NEWMARK
# load imputed data sets:
load(paste0("../../../Data/5imps.rda"))

# Due to interpretation and literature we recode month of birth to season of birth
for(i in 1:5){
  eval(parse(text=paste0("season <- imp",i, "$child_b_month")))
  levels(season) <- c("Winter", "Winter", rep("Spring", 3), rep("Summer", 3), rep("Autumn", 3), "Winter")
  # summer as reference category for dummy coding:
  season <- factor(season, levels = c( "Summer", "Spring","Autumn", "Winter"))
  eval(parse(text=paste0("imp",i, "$child_b_month <- NULL" ))) # drop month of birth variable
  eval(parse(text=paste0("imp",i, "$birth_b_season <- season" ))) # add season variable
}

# NEWMARK
# new: leave out a final test set
load("../../../Data/ind.part.rda")
test.fin.ind <- ind.part$nonfarm.inn
train.fin.ind <- setdiff(ind.part$nonfarm, ind.part$nonfarm.inn)
library(parallel)

do.for.imp.set <- function(i){
  
  # NEWMARK

  
 
  
  
  eval(parse(text=paste0("x.env.imp <- imp",i)))
  
  # NEWMARK
  # x.env.imp$farm <- NULL
  x.env.imp$center <- NULL
  
  # NEWMARK
  # x.env <- model.matrix(~. # + farm:Sex_female
  # ,x.env.imp)[,-1]
  
  # NEWMARK: for RF we can let factors be factors
  x.env <- data.frame(x.env.imp#, farm.Sex_female = as.factor(model.matrix(~.+ farm:Sex_female, x.env.imp)[ ,"farm2:Sex_female2"])
                      )
  #                       
  
  
  
  ################# adjust snp.data, env.data and outcome (same number of observations etc)
  # snp.data  <- dat.lit.snps.19
  
  ######### change: we define several data sets replacing x.env
  
  ### confounders: always to include - center, age, month of birth
  cnames.env <- colnames(x.env)
  confs <- cnames.env[c(grep("center", cnames.env),  #grep("month", cnames.env), (let month of birth stay environmental)
                        which(cnames.env=="age"), 
                        which(cnames.env=="child_age"), which(cnames.env=="child_BMI"), which(cnames.env=="Sex_female"))]
  
  x.conf <- x.env[ , confs]
  x.env <- x.env[ , -which(cnames.env %in% confs)]
  
  
  ## 1. x.env.nofh (no family anamnese)
  fam.an <- c("fhasthma", "fhhayfev", "fheczema", "FHx_Atopy" )
  x.env.nofh <-  x.env[ , - which(colnames(x.env) %in% fam.an)]
  
  ## 2. x.fh (family anamnese)
  x.fh <- x.env[ , fam.an]
  
 
  
  
  rm(x.env)
  
  # outcome
  get.y <- function(y.label){
    if(y.label == "dda") y <- dat.outcome$dd_asthma
    if(y.label == "dda_notany") y <- dat.outcome$dda_notany
    if(y.label == "dda.wh12_notany") y <- dat.outcome$dda.wh12_notany
    if(y.label == "asthma_wo_bronch") y <- dat.outcome$asthma_wo_bronch
    
    eval(parse(text=paste0("data.frame(",y.label,"=y)"))) 
  }
  
  # env-type
  # get.env <- function(env.label){
  #   if(env.label == "env_only") x.env <- x.env.nofh
  #   if(env.label == "fh_only") x.env <- x.fh
  #   if(env.label == "interactions_only") x.env <- x.interact
  #   if(env.label == "both") x.env <- x.env.3blocks
  #   if(env.label == "snps_only") x.env <- NULL
  #   x.env 
  # }
  
  
  # weights:
  
  
  ###### application:
  
  # NEWMARK
  # test.fin.ind <- test.farm.inn
  
  
  
  
  y.label="dda"
  obs.weights=gabriel.1707$weight_dna
  
  ######################################################################################################
  ################################## dda: ##############################################################
  ######################################################################################################
  
  y <- get.y(y.label)
  
  y.which <- colnames(y)
  y <- y[ ,1]
  
  
  #### note: we don't do this here since our final indices refere to the WHOLE data sets: we subset later s. below    
  #     ind.y <- which(!is.na(y)) 
  #     y <- y[ind.y]
  
  # 
  # x.env <- x.env
  #### s. 3-wrp-learn-on-snps-env.R   
  ########## create x.snp
  # 
   # x.snp <-  model.matrix(~.+0, data.frame(snp.data)[,-1])
  
  
  
  #### remove observations where outcome is NA ## should be done before
  # obs weights
  if(length(obs.weights)==1){
    obs.weights <- rep(1, length(y))
  } else{
    obs.weights <- obs.weights
  }
  
 
    x.all <- cbind(x.conf, x.env.nofh, x.fh)
  
  
  ################# set-up for algorithms
  # form <- formula( ~. +0)
  dat <- data.frame(y=factor(y), x.all)
  
  # if x argument has to be set, make it formula specific
  # x.pred <- model.matrix(form, data.frame(x.all))
  
  
  
  
  options(na.action='na.omit')
  
  ind.non.na <- which(!is.na(y))
  
  ### if missing values were removed (bec. of outcome)--> we have to adjust final training and test blocks
  train.fin.ind2 <- intersect(train.fin.ind, ind.non.na)
  test.fin.ind2 <- intersect(test.fin.ind, ind.non.na)
  
  
  ntrees=2000
  set.seed(23423)
  train.ranger <- ranger(dependent.variable.name = "y",  data = dat[ train.fin.ind2, ], write.forest = T, 
                         save.memory = F, probability = TRUE, num.trees =ntrees, num.threads = 8, verbose = FALSE, importance = "permutation") 
  # add altmann pvalue calculation for variance importance
  importance.altmann <- importance_pvalues(train.ranger, method="altmann", formula = formula(y~.), data = dat[ train.fin.ind2, ])
  
  # add prediction for innsbruck
  pred.rf <- predict(train.ranger,  data= data.frame(dat[ test.fin.ind2, ]) ) $predictions[ ,2]
  l <- list()
  l$train.ranger <- train.ranger
  l$importance.altmann <- importance.altmann
  l$pred.rf <- pred.rf
  l
}
obj.fin <- mclapply(1:5, function(x) do.for.imp.set(x), mc.cores=5) 

########################################## validation Innsbruck ############################################################################
####### validation
response.inn = dat.outcome$dd_asthma[test.fin.ind]



# unweighted
print(roc.a.rf <- lapply(1:5, function(z) loss.roc(y=response.inn, pred=obj.fin[[z]]$pred.rf, ci=TRUE)))


print(roc.b.i.rf <- mclapply(1:5, function(z) ci.auc.bootstrap.w(roc=roc.a.rf[[z]], boot.n=10000, parallel = FALSE, conf.level=.95, weights=gabriel.1707$weight_dna[test.fin.ind]) , mc.cores = 5 ))

# just for comparison
mult.preds <- lapply(1:5,  function(z) multiply.obs(y=response.inn, prediction=obj.fin[[z]]$pred.rf, weight.evals = gabriel.1707$weight_dna[test.fin.ind]))
print (roc.b.ii.rf <- lapply(1:5,  function(z) loss.roc(y= mult.preds[[z]][ , "y"], pred=mult.preds[[z]][ , "predicted"], ci=TRUE)))



save(obj.fin, 
     roc.a.rf, roc.b.i.rf, roc.b.ii.rf,
     # roc.a.ipf, roc.b.i.ipf, roc.b.ii.ipf,
     file=paste0("Results19/finalmodel/fin.models.", "dda" , "_", "fh_conf_env" ,"_", "ipf_lassoANDranger",".rda" )) 
