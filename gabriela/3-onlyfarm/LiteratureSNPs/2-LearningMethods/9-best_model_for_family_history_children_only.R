

### to be careful about here:

# each time another env block is used as env.only --> use another 'both' so as not to overwrite sth (or solve this in any way)
# for combination approach:
### just cbind the env in order: fam, env.new, interactions 
### bild list with the right.names and give it to function as argument

# be aware of: interactions have the highest dimension!! if a big selection of snps is used then something has to be changed
# however!!: we never can take into account so many interactions so we should stay with the 19 snps here anyway!!


### MARKNEW
### lit snps
setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/3-onlyfarm/LiteratureSNPs/2-LearningMethods/")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/3-wrp-learn-on-snp-env.R")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/0-functions.R")


library(parallel)



# load imputed data sets:
load(paste0("../../../Data/5imps.rda"))

### data
load("../../../LiteratureSNPs/Data/dat.lit.snps.19.RData")
load("../../../Data/gabriel.1707.Rdata")
load("../../../Data/dat.outcome.rda")
# load("../../Data/x.env.imp.rda")
load("../../../Data/x.wheeze_onset.rda")

# new: leave out a final test set
# load("../../Data/fin.ind.rda")

# NEWMARK
# take only farm and no innsbruck
load("../../../Data/ind.part.rda")
train.fin.ind <- setdiff(ind.part$farm, ind.part$farm.inn)

per.imp <- list()
for(i in 1:5){
  eval(parse(text=paste0("x.env.imp <- imp",i)))
  #### MARKNEW
  ## remove variables farm and center for the corresponding subgroupanalysis
  x.env.imp$center <- NULL
  x.env.imp$farm <- NULL
  
  #### MARKNEW
  # thus no interaction with farm either:
  x.env <- model.matrix(~. # + farm:Sex_female
                        ,x.env.imp)[,-1][ train.fin.ind,]
  
  # outcome:
  
  
  ################# adjust snp.data, env.data and outcome (same number of observations etc)
  snp.data  <- dat.lit.snps.19[ train.fin.ind,]
  
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
  form.interact.char <- paste( "~(", paste(colnames(snp.data[,-1]), collapse="+"), ")", ":", 
                               "(", 
                               paste(c("Num_Sibs_12", # "farm2", "farm2:Sex_female2",
                                       fam.an ), collapse="+"),
                               ")")
  form.interact <- formula(form.interact.char)
  
  x.interact <- model.matrix(  form.interact,
                               data= data.frame(snp.data[ ,-1], x.env[ , c("Num_Sibs_12", # "farm2", 
                                                                           fam.an )], x.conf[ , "Sex_female2", drop=FALSE] ))
  
  x.env.3blocks <- cbind(x.env.nofh, x.fh, x.interact)
  
  
  
  
  ## old (age of onset wheeze-variable)
  # x.wo <- as.matrix(x.wheeze_onset)[train.fin.ind,,drop=FALSE ]
  
  
  
  rm(x.env)
  
  # outcome
  get.y <- function(y.label){
    if(y.label == "dda") y <- dat.outcome$dd_asthma
    if(y.label == "dda_notany") y <- dat.outcome$dda_notany
    if(y.label == "dda.wh12_notany") y <- dat.outcome$dda.wh12_notany
    
    y <- y[ train.fin.ind]
    
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
  
  get.blocks.ipf <- function(env.label){
    if(env.label == "both"){
      blocks.ipf <- list(1:ncol(x.conf), 
                         (ncol(x.conf)+1):(ncol(x.conf) + ncol(x.env.nofh)),
                         (ncol(x.conf) + ncol(x.env.nofh) + 1): (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh)),
                         (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh)+1): (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh) + ncol(x.interact)),
                         (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh) + ncol(x.interact)+1):
                           (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh) + ncol(x.interact)+ncol(snp.data[,-1])))
    } else if(env.label == "fh_conf"){
      blocks.ipf <- list(1:ncol(x.fh), 
                         (ncol(x.fh)+1):(ncol(x.fh) + ncol(x.conf)) )
    } else if(env.label == "fh_conf_env"){
      blocks.ipf <- list(1:ncol(x.fh),
                         (ncol(x.fh)+1):(ncol(x.fh) + ncol(x.conf)),
                         (ncol(x.fh) + ncol(x.conf) +1): (ncol(x.fh) + ncol(x.conf) + ncol(x.env.nofh) )
      )
    } else if(env.label == "fh_conf_snps"){ 
      blocks.ipf <- list(1:ncol(x.fh),
                         (ncol(x.fh)+1):(ncol(x.fh) + ncol(x.conf)),
                         (ncol(x.fh) + ncol(x.conf) +1): (ncol(x.fh) + ncol(x.conf) + ncol(snp.data[,-1]) )
      )
    }else if(env.label == "fh_conf_env_snps"){
      blocks.ipf <- list(1:ncol(x.fh),
                         (ncol(x.fh)+1):(ncol(x.fh) + ncol(x.conf)),
                         (ncol(x.fh) + ncol(x.conf) +1): (ncol(x.fh) + ncol(x.conf) + ncol(x.env.nofh)), 
                         (ncol(x.fh) + ncol(x.conf) + ncol(x.env.nofh) +1): (ncol(x.fh) + ncol(x.conf) + ncol(x.env.nofh) + ncol(snp.data[,-1]))                
      )
    }
    blocks.ipf
  }
  
  # weights:
  obs.weights <- gabriel.1707$weight_dna[ train.fin.ind]
  
  ###### application:
  
  y.label = "dda"
  y=get.y(y.label)
  env.ifs <- c("fh_conf_snps")
  
  
  methods <-   "ranger"
  
  x.snp <- snp.data[ ,-1]
  
  x.all = cbind(get.env(env.ifs), x.snp)
  
  dat.temp <- data.frame(y=factor(y$dda), x.all)
  
  weights.fh1 <- obs.weights[dat.temp$fhasthma2==1]
  
  dat.fh1 <- dat.temp[dat.temp$fhasthma2==1, ]
  dat.fh1$fhasthma2 <- NULL
  dat.fh1$fhhayfev2 <- NULL
  dat.fh1$FHx_Atopy2 <- NULL
  
  head(dat.fh1)
  
  k=5
  # 5-Fold-CV:
  k <- k
  n <- nrow(dat.fh1)
  set.seed(3534)
  ind <- sample(1:n)
  
  ind.test <- list()
  ind.train <- list()
  for(j in 1:k){
    if(j < k ) ind.test[[j]] <- ind[(1+(j-1)*floor(n/k)):(j*(floor(n/k))) ]
    else ind.test[[j]] <- ind[(1+(j-1)*floor(n/k)):n]
    ind.train[[j]] <- setdiff(ind, ind.test[[j]])
  }
  
  l <- list()
  for(j in 1:k){
    ntrees=2000
    l$train.ranger[[j]] <- train.ranger <- ranger(dependent.variable.name = "y", data = dat.fh1[ ind.train[[j]], ], write.forest = T, 
                           save.memory = FALSE, probability = TRUE, num.trees =ntrees, num.threads = 8, verbose = FALSE
    ) 
    l$pred[[j]] <- predict(train.ranger, data= dat.fh1[ind.test[[j]], ])$predictions[ ,2]
    l$y[[j]] <- dat.fh1$y[ind.test[[j]]]
    l$w[[j]] <- weights.fh1[ ind.test[[j]] ]
    print(roc(  l$y[[j]] , l$pred[[j]] , direction = "<" ))
  }
  per.imp[[i]] <- l
}  
  
  
  


aucs.imps <- list()
for(i in 1:5){
  pred.imp_i <- unlist(per.imp[[i]]$pred)
  y.imp_i <- unlist(per.imp[[i]]$y)
  w.imp_i <- unlist(per.imp[[i]]$w)
  aucs.imps[[i]] <-  boot_AUC(y.imp_i, pred.imp_i, nboot = 10000, weights = w.imp_i, cores = 1  )
}

quantile(aucs.imps[[3]]$AUC, prob = c(0.025, .5, 0.975)  )

save(aucs.imps, file="Results19/aucs/aucs_fh1_ranger_fh_conf_snps.rda")

