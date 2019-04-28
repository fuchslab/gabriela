# outcome:

#### pipeline: 
# 1. get AUCs aggregated from imputations for different outcomes
# 2. get betas/importances from imputations for different outcomes


# NEWMARK
### lit snps
setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/3-onlyfarm/LiteratureSNPs/2-LearningMethods/")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/3-wrp-learn-on-snp-env.R")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/0-functions.R")

library("ranger")

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
library(parallel)
# NEWMARK
train.fin.ind <- setdiff(ind.part$farm, ind.part$farm.inn)
test.fin.ind <- ind.part$farm.inn

do.for.imp.set <- function(i){
  
  
  
  
  
  eval(parse(text=paste0("x.env.imp <- imp",i)))
  
  # NEWMARK
  x.env.imp$farm <- NULL
  x.env.imp$center <- NULL
  
  # NEWMARK
  x.env <- model.matrix(~. # + farm:Sex_female
                        ,x.env.imp)[,-1]
  

  
  ################# adjust snp.data, env.data and outcome (same number of observations etc)
  snp.data  <- dat.lit.snps.19
  
  ######### change: we define several data sets replacing x.env
  
  ### confounders: always to include - center, age, month of birth
  cnames.env <- colnames(x.env)
  confs <- cnames.env[c(grep("center", cnames.env),  #grep("month", cnames.env), (let month of birth stay environmental)
                        which(cnames.env=="age"), 
                        which(cnames.env=="child_age"), which(cnames.env=="child_BMI"), which(cnames.env=="Sex_female2"))]
  
  x.conf <- x.env[ , confs]
  x.env <- x.env[ , -which(cnames.env %in% confs)]
  
  
  ## 1. x.env.nofh (no family anamnese)
  fam.an <- c("fhasthma2", "fhhayfev2", "fheczema2", "FHx_Atopy2" )
  x.env.nofh <-  x.env[ , - which(colnames(x.env) %in% fam.an)]
  
  ## 2. x.fh (family anamnese)
  x.fh <- x.env[ , fam.an]
  
  ## 3. x.interact (several interactions of snps with fam-hist and env)
  # formula:
  
  # NEWMARK
  form.interact.char <- paste( "~(", paste(colnames(snp.data[,-1]), collapse="+"), ")", ":", 
                               "(", 
                               paste(c("Num_Sibs_12", 
                                       # "farm2", "farm2:Sex_female2",
                                       fam.an ), collapse="+"),
                               ")")
  form.interact <- formula(form.interact.char)
  
  # NEWMARK
  x.interact <- model.matrix(  form.interact,
                               data= data.frame(snp.data[ ,-1], x.env[ , c("Num_Sibs_12", 
                                                                           # "farm2", 
                                                                           fam.an )], x.conf[ , "Sex_female2", drop=FALSE] ))
  
  
  x.env.3blocks <- cbind(x.env.nofh, x.fh, x.interact)
  
  blocks.ipf <- list(1:ncol(x.fh),
                     (ncol(x.fh)+1):(ncol(x.fh) + ncol(x.conf)),
                     (ncol(x.fh) + ncol(x.conf) +1): (ncol(x.fh) + ncol(x.conf) + ncol(snp.data[,-1]) )
  )
  
  ## old (age of onset wheeze-variable)
  # x.wo <- as.matrix(x.wheeze_onset)[train.fin.ind,,drop=FALSE ]
  
  
  
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
  get.env <- function(env.label){
    if(env.label == "conf_only") x.env <- x.conf
    if(env.label == "env_only") x.env <- x.env.nofh
    if(env.label == "fh_only") x.env <- x.fh
    if(env.label == "interactions_only") x.env <- x.interact
    if(env.label == "both") x.env <- x.env.3blocks
    if(env.label == "snps_only") x.env <- NULL
    if(env.label == "fh_conf" | env.label == "fh_conf_snps") x.env <- cbind(x.fh, x.conf) 
    if(env.label == "fh_conf_env" | env.label=="fh_conf_env_snps") x.env <- cbind(x.fh, x.conf, x.env.nofh)
    x.env 
  }
  
  # weights:
  
  
  ###### application:
  
  # NEWMARK
  # test.fin.ind <- ind.part$farm.inn
  
  
  y.label="dda"
  method="ipflasso_weights"
  env.if = "fh_conf_snps"
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
  
  
  x.env <- get.env(env.if)
  #### s. 3-wrp-learn-on-snps-env.R   
  ########## create x.snp
  
  x.snp <-  model.matrix(~.+0, data.frame(snp.data)[,-1])
  
  
  if(env.if=="fh_conf_snps"){
    x.all = cbind(get.env(env.if), x.snp)
    pflist = list(c(1,1,1),  # first block: confounders and  weight always set to 1 
                  c(2,1,1), c(1,1,2), # weight one block more
                  c(2,1,2), # weight two blocks more
                  c(3,1,1), c(1,1,3), # weight one block three times
                  c(3,1,3) # weight two blocks three times
    )
  }
  
  #### remove observations where outcome is NA ## should be done before
  # obs weights
  if(length(obs.weights)==1){
    obs.weights <- rep(1, length(y))
  } else{
    obs.weights <- obs.weights
  }
  
  
  # create x.all
  if(env.if == "both"){
    x.all <- cbind(x.conf, x.env, x.snp)
    x.snp <- cbind(x.conf, x.snp) # since we need it for some pca-approaches
    x.env <- cbind(x.conf, x.env)
    
  } 
  if(env.if == "env_only" | env.if == "fh_only" |env.if ==  "interactions_only"){
    x.all <- cbind(x.conf, x.env)
    rm(x.snp)
  } 
  
  if(env.if == "snps_only"){
    x.all <- cbind(x.conf, x.snp)
    rm(x.env)
  } 
  ################# set-up for algorithms
  form <- formula( ~. +0)
  dat <- data.frame(y=factor(y), x.all)
  
  # if x argument has to be set, make it formula specific
  x.pred <- model.matrix(form, data.frame(x.all))
  
  
  
  
  options(na.action='na.omit')
  
  ind.non.na <- which(!is.na(y))
  
  ### if missing values were removed (bec. of outcome)--> we have to adjust final training and test blocks
  train.fin.ind2 <- intersect(train.fin.ind, ind.non.na)
  test.fin.ind2 <- intersect(test.fin.ind, ind.non.na)
  
  l <- list()
  ## the trained model:
   train.glmnet <- cvr2.ipflasso.ed(Y=y[train.fin.ind2], X=x.pred[ train.fin.ind2, ], family="binomial", type.measure="auc", alpha=1, 
   blocks=blocks.ipf, pflist=pflist,
   nfolds=5,ncv=5, weights=obs.weights[ train.fin.ind2])
   pred.ipf <- ipflasso.predict(train.glmnet, Xtest = x.pred[ test.fin.ind2, ])$probabilitiestest
   
   ################ compare this to ranger
   ntrees=2000
   set.seed(23423)
   dat <- data.frame(y=factor(y), x.pred)
   train.ranger <- ranger(dependent.variable.name = "y",  data = dat[ train.fin.ind2, ], write.forest = T, 
                          save.memory = F, probability = TRUE, num.trees =ntrees, num.threads = 8, verbose = FALSE, importance = "permutation"
                          ) 
   # add altmann pvalue calculation for variance importance
   importance.altmann <- importance_pvalues(train.ranger, method="altmann", formula = formula(y~.), data = dat[ train.fin.ind2, ])
   
   # add prediction for innsbruck
   pred.rf <- predict(train.ranger,  data= data.frame(dat[ test.fin.ind2, ]) ) $predictions[ ,2]
   l <- list()

  
  
  l$train.ranger <- train.ranger
  l$importance.altmann <- importance.altmann
  l$pred.rf <- pred.rf
  
  l$train.glmnet <- train.glmnet
  l$pred.ipf <- pred.ipf
  l
  
}
obj.fin <- mclapply(1:5, do.for.imp.set, mc.cores = 5)




 



########################################## validation Innsbruck ############################################################################
####### validation
response.inn = dat.outcome$dd_asthma[test.fin.ind]


##### random forest
# unweighted
print(roc.a.rf <- lapply(1:5, function(z) loss.roc(y=response.inn, pred=obj.fin[[z]]$pred.rf, ci=TRUE)))
# same for ipflasso
# print(roc.a.rf <- lapply(1:5, function(z) loss.roc(y=response.inn, pred=obj.fin[[z]]$pred.ipf, ci=TRUE)))

set.seed(23434)
print(roc.b.i.rf <- mclapply(1:5, function(z) ci.auc.bootstrap.w(roc=roc.a.rf[[z]], boot.n=10000, parallel = FALSE, conf.level=.95, 
                                                              n.cases.min = .1,
                                                              weights= gabriel.1707$weight_dna[test.fin.ind])
                          , mc.cores = 5 ))

# just for comparison
mult.preds <- lapply(1:5,  function(z) multiply.obs(y=response.inn, prediction=obj.fin[[z]]$pred.rf, weight.evals = gabriel.1707$weight_dna[test.fin.ind]))
print (roc.b.ii.rf <- lapply(1:5,  function(z) loss.roc(y= mult.preds[[z]][ , "y"], pred=mult.preds[[z]][ , "predicted"], ci=TRUE)))




##### ipflasso
# unweighted
print(roc.a.ipf <- lapply(1:5, function(z) loss.roc(y=response.inn, pred=obj.fin[[z]]$pred.ipf, ci=TRUE)))
# same for ipflasso
# print(roc.a.ipf <- lapply(1:5, function(z) loss.roc(y=response.inn, pred=obj.fin[[z]]$pred.ipf, ci=TRUE)))

set.seed(23434)
print(roc.b.i.ipf <- mclapply(1:5, function(z) ci.auc.bootstrap.w(roc=roc.a.ipf[[z]], boot.n=10000, parallel = FALSE, conf.level=.95, 
                                                                 n.cases.min = .1,
                                                                 weights= gabriel.1707$weight_dna[test.fin.ind])
                             , mc.cores = 5 ))

# just for comparison
mult.preds <- lapply(1:5,  function(z) multiply.obs(y=response.inn, prediction=obj.fin[[z]]$pred.ipf, weight.evals = gabriel.1707$weight_dna[test.fin.ind]))
print (roc.b.ii.ipf <- lapply(1:5,  function(z) loss.roc(y= mult.preds[[z]][ , "y"], pred=mult.preds[[z]][ , "predicted"], ci=TRUE)))




save(obj.fin, 
     roc.a.rf, roc.b.i.rf, roc.b.ii.rf,
     roc.a.ipf, roc.b.i.ipf, roc.b.ii.ipf,
     file=paste0("Results19/finalmodel/fin.models.", "dda" , "_", "fh_conf_snp" ,"_", "ipf_lassoANDranger",".rda" )) 



