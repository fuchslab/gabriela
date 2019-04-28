

### to be careful about here:

# each time another env block is used as env.only --> use another 'both' so as not to overwrite sth (or solve this in any way)
# for combination approach:
### just cbind the env in order: fam, env.new, interactions 
### bild list with the right.names and give it to function as argument

# be aware of: interactions have the highest dimension!! if a big selection of snps is used then something has to be changed
# however!!: we never can take into account so many interactions so we should stay with the 19 snps here anyway!!


# MARKNEW: change paths
### lit snps
setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/3-wrp-learn-on-snp-env.R")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/0-functions.R")


library(parallel)



# load imputed data sets:
load(paste0("Data/5imps.rda"))

### data
 load("LiteratureSNPs/Data/dat.lit.snps.19.RData")



load("ImputedSNPs/Data/dt.imp.pruned.rda")

load("Data/gabriel.1707.Rdata")
load("Data/dat.outcome.rda")
# load("../../Data/x.env.imp.rda")
load("Data/x.wheeze_onset.rda")

# # new: leave out a final test set
# load("Data/fin.ind.rda")


#### MARKNEW
# take only farm and no innsbruck
load("Data/ind.part.rda")
train.fin.ind <- setdiff(1:1707, ind.part$all.inn)






do.for.imp.set <- function(i,  env.ifs, snp.as){
  eval(parse(text=paste0("x.env.imp <- imp",i)))
  #### MARKNEW
  ## remove variables farm and center for the corresponding subgroupanalysis
  x.env.imp$center <- NULL

  
  #### MARKNEW
  # thus no interaction with farm either:
  x.env <- model.matrix(~.  + farm:Sex_female
                        ,x.env.imp)[,-1][ train.fin.ind,]
  
  # outcome:
  
  
################# adjust snp.data, env.data and outcome (same number of observations etc)
  snp.data  <- dt.imp.pruned[ train.fin.ind,]
  # snp.data for interactions
  snp.for.int <- dat.lit.snps.19[train.fin.ind, -1]
  
  
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
  
  #### MARKNEW
  # leave out farm-interactions
  ## 3. x.interact (several interactions of snps with fam-hist and env)
  # formula:
  form.interact.char <- paste( "~(", paste(colnames(snp.for.int), collapse="+"), ")", ":", 
                               "(", 
                               paste(c("Num_Sibs_12", "farm2", "farm2:Sex_female2", fam.an ), collapse="+"),
                               ")")
  form.interact <- formula(form.interact.char)
  
  x.interact <- model.matrix(  form.interact,
                               data= data.frame(snp.for.int, x.env[ , c("Num_Sibs_12", "farm2", fam.an )], x.conf[ , "Sex_female2", drop=FALSE] ))
  
  
  x.env.3blocks <- cbind(x.env.nofh, x.fh, x.interact)
  
  blocks.ipf <- list(1:ncol(x.conf), 
                     (ncol(x.conf)+1):(ncol(x.conf) + ncol(x.env.nofh)),
                     (ncol(x.conf) + ncol(x.env.nofh) + 1): (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh)),
                     (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh)+1): (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh) + ncol(x.interact)),
                     (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh) + ncol(x.interact)+1):
                       (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh) + ncol(x.interact)+ncol(snp.data))
  )
  
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
    if(env.label == "env_only") x.env <- x.env.nofh
    if(env.label == "fh_only") x.env <- x.fh
    if(env.label == "interactions_only") x.env <- x.interact
    if(env.label == "both") x.env <- x.env.3blocks
    if(env.label == "snps_only") x.env <- NULL
    x.env 
  }
  
  
  # weights:
  obs.weights <- gabriel.1707$weight_dna[ train.fin.ind]
  
  ###### application:
  
  y.label = "dda"
  
 
  # methods <- c("glmnet", "glmnet_weights", "elnet", "elnet_weights", "ipflasso", "ipflasso_weights",
  #              "ranger", "pca_snps_rf",  "pca_rf")
  methods <- "ranger"
  
  
  ## test:
  # 
  # y.label="dda.wh12_notany"
  # env.if <- "both"
  # method = "ipflasso_weights"
  # 
  # 
  # method="glmnet_weights"
  # perform.asthma(y=get.y(y.label), snp.as="estdos", env.if=env.if, 
  #                method=method, x.env=get.env(env.if), snp.data=snp.data, 
  #                obs.weights=obs.weights, x.conf=x.conf, blocks.ipf=blocks.ipf,
  #                path.results="Results19/" )
  # 
  
  
  ### Perform:

  for(method in methods){
    for(env.if in env.ifs){
      # for(y.label in y.labels){
      perform.asthma(y=get.y(y.label), snp.as=snp.as, env.if=env.if, 
                     method=method, x.env=get.env(env.if), snp.data=snp.data, 
                     obs.weights=obs.weights, x.conf=x.conf, blocks.ipf=blocks.ipf,
                     #### MARKNEW
                     ## change path!
                     path.results=paste0("2-validate-innsbruck/ImputedSNPs/Results_genomewide/imp",i,"/"), 
                     needed.as.matrix=FALSE, ntrees=20000 )
    } 
    # print(paste(method, env.if, ": for all outcomes performed"))
    # } 
    print(paste(method, "for all outcomes performed and SNP, env and both performed"))
  }
}


 mclapply(1:5, function(x) do.for.imp.set(x,  env.ifs = "snps_only", snp.as = "genomewide"), mc.cores=5)
 # mclapply(1:5, function(x) do.for.imp.set(x,  env.ifs = "both", snp.as = "genomewide"), mc.cores=5)


