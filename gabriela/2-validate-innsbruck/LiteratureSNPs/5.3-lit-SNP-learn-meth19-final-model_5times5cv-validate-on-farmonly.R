# outcome:

#### pipeline: 
# 1. get AUCs aggregated from imputations for different outcomes
# 2. get betas/importances from imputations for different outcomes


# NEWMARK
### lit snps
setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/2-validate-innsbruck/LiteratureSNPs/2-LearningMethods/")
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

do.for.imp.set <- function(i, env.if="both", method2=NULL){
  
  # NEWMARK
  train.fin.ind <- setdiff(1:1707, ind.part$all.inn)
  
  test.farm.inn <- ind.part$farm.inn
  
  
  eval(parse(text=paste0("x.env.imp <- imp",i)))
  
  # NEWMARK
  # x.env.imp$farm <- NULL
  x.env.imp$center <- NULL
  
  # NEWMARK
  x.env <- model.matrix(~.  + farm:Sex_female
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
                                        "farm2", "farm2:Sex_female2",
                                       fam.an ), collapse="+"),
                               ")")
  form.interact <- formula(form.interact.char)
  
  # NEWMARK
  x.interact <- model.matrix(  form.interact,
                               data= data.frame(snp.data[ ,-1], x.env[ , c("Num_Sibs_12", 
                                                                            "farm2", 
                                                                           fam.an )], x.conf[ , "Sex_female2", drop=FALSE] ))
  
  
  x.env.3blocks <- cbind(x.env.nofh, x.fh, x.interact)
  
  blocks.ipf <- list(1:ncol(x.conf), # confounders
                     (ncol(x.conf)+1):(ncol(x.conf) + ncol(x.env.nofh)), # environmnet
                     (ncol(x.conf) + ncol(x.env.nofh) + 1): (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh)), # family history
                     (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh)+1): (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh) + ncol(x.interact)), # interactions
                     (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh) + ncol(x.interact)+1):
                       (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh) + ncol(x.interact)+ncol(snp.data[,-1])) # snps
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
    if(env.label == "env_only") x.env <- x.env.nofh
    if(env.label == "fh_only") x.env <- x.fh
    if(env.label == "interactions_only") x.env <- x.interact
    if(env.label == "both") x.env <- x.env.3blocks
    if(env.label == "snps_only") x.env <- NULL
    x.env 
  }
  
  
  # weights:
  
  
  ###### application:
  
  # NEWMARK
  test.fin.ind <- test.farm.inn
  
  
  

    y.label="dda"
    method="ipflasso_weights"
    env.if = "both"
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
    
    
      ## the trained model:
      # train.glmnet <- cvr2.ipflasso.ed(Y=y[train.fin.ind2], X=x.pred[ train.fin.ind2, ], family="binomial", type.measure="auc", alpha=1, 
                                       # blocks=blocks, pflist=pflist,
                                       # nfolds=5,ncv=5, weights=obs.weights[ train.fin.ind2])
    #... can just be loaded:
    load(paste0("Results19/finalmodel/fin.models.", y.label, ".imp",i,"_weightscv5times.rda" ))
    
    pred <- ipflasso.predict(obj.fin.both$model.obj, Xtest = x.pred[ test.fin.ind2, ])$probabilitiestest
    
    }
    
    
    
    
    ##################################
    ####### validation
    response = y[test.fin.ind2]
    # unweighted
    print(roc.a <- loss.roc(y=response, pred=pred, ci=TRUE))
    
    print(roc.b.i <- ci.auc.bootstrap.w(roc=roc.a, boot.n=10000, parallel = FALSE, conf.level=.95, weights=obs.weights[test.fin.ind2]))
    
    mult.preds <- multiply.obs(y=response, prediction=pred, weight.evals = obs.weights[test.fin.ind2])
    print (roc.b.ii <- loss.roc(y= mult.preds[ , "y"], pred=mult.preds[ , "predicted"], ci=TRUE))
    
    class(roc.b.ii)
    prediction.object <- prediction(mult.preds[ , "predicted"], mult.preds[ , "y"])
    perf <- performance(prediction.object,"tpr","fpr")
     plot(perf)
    # list(model.obj=obj.fin.both$model.obj, pred = pred, 
         # perf.values = list(roc.a=roc.a, roc.b.i = roc.b.i, roc.b.ii, mult.preds=mult.preds ))
  
  
  

      obj.fin.both <- model.19.final(env.if="both", y.label=y.label, method = "ipflasso_weights")



method2 = NULL
mclapply(1:5, function(x) do.for.imp.set(x, env.if = "both", method2 = method2), mc.cores=5) 
# mclapply(1:5, function(x) do.for.imp.set(x, env.if = "fh_only", method2 = method2), mc.cores=5)
# 
# method2 = "elnet_weights"
# mclapply(1:5, function(x) do.for.imp.set(x, env.if = "snps_only", method2 = method2), mc.cores=5)
# mclapply(1:5, function(x) do.for.imp.set(x, env.if = "fh_only", method2 = method2), mc.cores=5)
# mclapply(1:5, function(x) do.for.imp.set(x, env.if = "env_only", method2 = method2), mc.cores=5)

# run time (on 3*5 (outcomes * imputations) kernels):
# started at 16:45



