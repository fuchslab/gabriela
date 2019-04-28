# NEWMARK
### lit snps
setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/2-validate-innsbruck/LiteratureSNPs/2-LearningMethods/")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/3-wrp-learn-on-snp-env.R")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/0-functions.R")

library("ranger")

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
test.fin.ind <- ind.part$all.inn
train.fin.ind <- setdiff(1:1707, ind.part$all.inn)
library(parallel)


env.if = "fh_conf_env_snps_farm.interact1"

do.for.imp.set <- function(i){
  
  # NEWMARK

  

  
  
  eval(parse(text=paste0("x.env.imp <- imp",i)))
  
  # NEWMARK
  # x.env.imp$farm <- NULL
  x.env.imp$center <- NULL
  
  # NEWMARK
  x.env <- model.matrix(~. # + farm:Sex_female
  ,x.env.imp)[,-1]
  
  # NEWMARK: for RF we can let factors be factors
  # x.env <- data.frame(x.env.imp, farm.Sex_female = as.factor(model.matrix(~.+ farm:Sex_female, x.env.imp)[ ,"farm2:Sex_female2"]))
  #                       
  ################# adjust snp.data, env.data and outcome (same number of observations etc)
  snp.data  <- dat.lit.snps.19
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
  fam.an <- c("fhasthma2", "fhhayfev2", "fheczema2", "FHx_Atopy2" )
  x.env.nofh <-  x.env[ , - which(colnames(x.env) %in% fam.an)]
  
  ## 2. x.fh (family anamnese)
  x.fh <- x.env[ , fam.an]
  
  
  
  # NEWMARK
  ##################################### formula farm interactions 1
  # since fh_conf_snps best for onlyfarm (and after decrease)
  # for version 1: just include all predictors as farm interactions
  
  
  form.farm.interact1.char <- paste( "~(", paste(colnames(cbind( x.fh,x.conf, x.env.nofh[,- which(colnames(x.env.nofh)=="farm2")],  snp.data[,-1])), collapse="+"), ")", ":", 
                                     "farm2"
  )
  
  form.farm.interact1 <- formula(form.farm.interact1.char)
  
  # NEWMARK
  x.farm.interact1 <- model.matrix(   form.farm.interact1,
                                      data= data.frame(cbind( x.fh,x.conf, x.env.nofh,  snp.data[,-1])))
  
  
  x.env.3blocks.farm.interact1 <- cbind(x.env.nofh, x.fh, x.farm.interact1)
  
  blocks.farm.interact1.ipf <- list(1:ncol(x.conf), # confounders
                                    (ncol(x.conf)+1):(ncol(x.conf) + ncol(x.env.nofh)), # environmnet
                                    (ncol(x.conf) + ncol(x.env.nofh) + 1): (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh)), # family history
                                    (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh)+1): (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh) + ncol(x.farm.interact1)), # interactions
                                    (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh) + ncol(x.farm.interact1)+1):
                                      (ncol(x.conf) + ncol(x.env.nofh) + ncol(x.fh) + ncol(x.farm.interact1)+ncol(snp.data[,-1])) # snps
  )
  
  
  
  
  
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
    if(env.label=="fh_conf_env_snps_farm.interact1") x.env <- x.env.3blocks.farm.interact1
    if(env.label=="fh_conf_env_snps_farm.interact2") x.env <- x.env.3blocks.farm.interact2
    x.env 
  }
  
  
  
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
  
  # x.snp <-  model.matrix(~.+0, data.frame(snp.data)[,-1])
  x.snp <-  model.matrix(~.+0, data.frame(snp.data)[,-1])
  
  
  #### remove observations where outcome is NA ## should be done before
  # obs weights
  if(length(obs.weights)==1){
    obs.weights <- rep(1, length(y))
  } else{
    obs.weights <- obs.weights
  }
  
 
    x.all <- cbind(x.conf, get.env(env.label =env.if), x.snp)
  
  
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
  pred.inn <- predict(train.ranger,  data= data.frame(dat[ test.fin.ind2, ]) ) $predictions[ ,2]
  l <- list()
  l$train.ranger <- train.ranger
  l$importance.altmann <- importance.altmann
  l$pred.inn <- pred.inn
  l
}
obj.fin <- mclapply(1:5, function(x) do.for.imp.set(x), mc.cores=5) 

save(obj.fin, file=paste0("Results19/finalmodel/fin.models.", "dda" , "_", env.if ,"_", "ranger",".rda" ))
load(file=paste0("Results19/finalmodel/fin.models.", "dda" , "_", env.if ,"_", "ranger",".rda" ))
########################################## validation Innsbruck ############################################################################
####### validation
response.inn = dat.outcome$dd_asthma[test.fin.ind]



# unweighted
 print(roc.a.farm.inn <- lapply(1:5, function(z) loss.roc(y=response.inn, pred=obj.fin[[z]]$pred.inn, ci=TRUE)))

 
print(roc.b.i <- mclapply(1:5, function(z) ci.auc.bootstrap.w(roc=roc.a[[z]], boot.n=10000, parallel = FALSE, conf.level=.95, weights=gabriel.1707$weight_dna[test.fin.ind]) , mc.cores = 5 ))

# just for comparison
mult.preds <- lapply(1:5,  function(z) multiply.obs(y=response.inn, prediction=obj.fin[[z]]$pred.inn, weight.evals = gabriel.1707$weight_dna[test.fin.ind]))
print (roc.b.ii <- lapply(1:5,  function(z) loss.roc(y= mult.preds[[z]][ , "y"], pred=mult.preds[[z]][ , "predicted"], ci=TRUE)))

# class(roc.b.ii)
# prediction.object <- prediction(mult.preds[ , "predicted"], mult.preds[ , "y"])
# perf <- performance(prediction.object,"tpr","fpr")
# plot(perf)


### visualize
(inn.auc.mean <- mean(sapply(1:5, function(z) roc.b.i[[z]][2])))
# compare to replicated obs-AUC
mean(sapply(1:5, function(z) roc.b.ii[[z]]$auc)) # check

inn.auc.low <- mean(sapply(1:5, function(z) roc.b.i[[z]][1]))
inn.auc.up <- mean(sapply(1:5, function(z) roc.b.i[[z]][3]))


# pdf("../../../Results/2-validate-innsbruck/ROC_AUC_5imps_innsbruck.pdf")
plot(roc.b.ii[[1]], main = "Validation on GABRIELA (Innsbruck) - ROC curves for 5 imputations")
lapply(2:5, function(z) plot(roc.b.ii[[z]], add=T))
text(0.25,0.25, paste("AUC =", round(inn.auc.mean, 3), "(", round(inn.auc.low, 3) ,"-" , round(inn.auc.up, 3) , ")"))
# dev.off()


################### validate on farm only:
ind.farm.in.inn <- which(ind.part$all.inn %in% ind.part$farm.inn)


print(roc.a.farm.inn <- lapply(1:5, function(z) loss.roc(y= dat.outcome$dd_asthma[ind.part$farm.inn], pred=obj.fin[[z]]$pred.inn[ind.farm.in.inn], ci=TRUE)))
print(roc.b.i.farm.inn <- mclapply(1:5, function(z) ci.auc.bootstrap.w(roc=roc.a.farm.inn[[z]], boot.n=10000, parallel = FALSE, conf.level=.95, weights=gabriel.1707$weight_dna[ind.farm.in.inn]) , mc.cores = 5 ))
#################### validate on nonformonly


########################################## Validation PASTURE (not applicabel with snps) ####################################################################################
# 
# 
# # load imputed data sets:
# load(paste0("../../../5-validate-pasture/Data/5past.imps.rda"))
# load(paste0("../../../5-validate-pasture/Data/outcomes.rdata"))
# # load("../../../")
# # all(outcomes$idx == pasture$idx) #  check
# ####### data preparations for modeling
# 
# # Due to interpretation and literature we recode month of birth to season of birth
# for(i in 1:5){
#   eval(parse(text=paste0("season <- past.imp",i, "$child_b_month")))
#   levels(season) <- c("Winter", "Winter", rep("Spring", 3), rep("Summer", 3), rep("Autumn", 3), "Winter")
#   # summer as reference category for dummy coding:
#   season <- factor(season, levels = c( "Summer", "Spring","Autumn", "Winter"))
#   eval(parse(text=paste0("past.imp",i, "$child_b_month <- NULL" ))) # drop month of birth variable
#   eval(parse(text=paste0("past.imp",i, "$birth_b_season <- season" ))) # add season variable
# }
# 
# ## final data set preparation for RF (ranger)
# do.for.past.imp.set <- function(i){
#   
# 
#   eval(parse(text=paste0("x.env.imp <- past.imp",i)))
#   
#   
#   # NEWMARK
#   # x.env <- model.matrix(~.+ farm:Sex_female
#                         # ,x.env.imp)[,-1]
#   
#   # NEWMARK: for RF we can let factors be factors
#   x.env <- data.frame(x.env.imp, farm.Sex_female = as.factor(model.matrix(~.+ farm:Sex_female, x.env.imp)[ ,"farm2:Sex_female2"]))
#                         #                       
#   
#   
#   ######### change: we define several data sets replacing x.env
#   
#   ### confounders: always to include - center, age, month of birth
#   cnames.env <- colnames(x.env)
#   confs <- cnames.env[c(grep("center", cnames.env),  #grep("month", cnames.env), (let month of birth stay environmental)
#                         which(cnames.env=="age"), 
#                         which(cnames.env=="child_age"), which(cnames.env=="child_BMI"), which(cnames.env=="Sex_female"))]
#   
#   x.conf <- x.env[ , confs]
#   x.env <- x.env[ , -which(cnames.env %in% confs)]
#   
#   
#   ## 1. x.env.nofh (no family anamnese)
#   fam.an <- c("fhasthma", "fhhayfev", "fheczema", "FHx_Atopy" )
#   x.env.nofh <-  x.env[ , - which(colnames(x.env) %in% fam.an)]
#   
#   ## 2. x.fh (family anamnese)
#   x.fh <- x.env[ , fam.an]
#   
#   
#   
#   x.all <- cbind(x.conf, x.env.nofh, x.fh)
#   
#   
#   ################# set-up for algorithms
#   # form <- formula( ~. +0)
#   
#   y <- outcomes$asthma
#   
#   na.omit(data.frame(y=factor(y), x.all))
# }
# 
# past.imp.env <- lapply(1:5, do.for.past.imp.set)
# str(past.imp.env,1)
# 
# head(past.imp.env[[1]])
# 


################### do the prediction

# we assign one learnt ranger models on  one of the 5 imputations randomly to one of the 5 imputations of PASTURE (or in this case, 1 to 1, 2 to 2,...)
# pred.past <- lapply(1:5, function(z) predict(obj.fin[[z]]$train.ranger,  data= data.frame(past.imp.env[[z]]) )$predictions[ ,2])


################### do the validation

# 
# response.past = as.factor(na.omit(outcomes$asthma))
# 
# # life is easy here; no weighted ROC required
# print(roc.a.past <- lapply(1:5, function(z) loss.roc(y=response.past, pred=pred.past[[z]], ci=TRUE)))
# 
# ### visualize
# past.auc.mean <- mean(sapply(1:5, function(z) roc.a.past[[z]]$auc))
# past.auc.low <- mean(sapply(1:5, function(z) roc.a.past[[z]]$ci[1]))
# past.auc.up <- mean(sapply(1:5, function(z) roc.a.past[[z]]$ci[3]))
# 
# # pdf("../../../Results/5-validate-pasture/ROC_AUC_5imps_pasture.pdf")
# plot(roc.a.past[[1]], main = "Validation on PASTURE - ROC curves for 5 imputations")
# lapply(2:5, function(z) plot(roc.a.past[[z]], add=T))
# text(0.25,0.25, paste("AUC =", round(past.auc.mean, 3), "(", round(past.auc.low, 3) ,"-" , round(past.auc.up, 3) , ")"))
# # dev.off()

################################################# variable importance ###################################################################
vars.imp <- do.call(rbind, lapply(1:5, function(z) importance(obj.fin[[z]]$train.ranger)))

means.vars.imp <- colMeans(vars.imp)
means.vars.imp.sort <- sort(means.vars.imp, decreasing = F)
names.vars.sort <- names(sort(means.vars.imp, decreasing = F))

## plot best variables as boxplots
library(Hmisc)
sd.vars.imp.sort <- apply(vars.imp[ , names.vars.sort], 2, sd)

# pdf("../../../Results/2-validate-innsbruck/RF_final_model_fh_conf_env_all_variableimportance.pdf", height=9)
# pdf("3-onlyfarm/Results/ranger_snps_only_farm-only_final_importance.pdf", height=9)
errbar(x=names.vars.sort, y=means.vars.imp.sort, yplus = means.vars.imp.sort + sd.vars.imp.sort, yminus = means.vars.imp.sort- sd.vars.imp.sort)
abline(v=0)
# dev.off()



################### pvalues altmann

vars.imp.altmann <- do.call(rbind, lapply(1:5, function(z) obj.fin[[z]]$importance.altmann[ ,"importance" ]))

### the following code in comments gives the same as above...
 means.vars.imp.altmann <- colMeans(vars.imp.altmann)
# means.vars.imp.altmann.sort <- sort(means.vars.imp.altmann, decreasing = F)
 names.vars.sort.altmann <- names(sort(means.vars.imp.altmann, decreasing = F))
# 
# ## plot best variables as boxplots
# library(Hmisc)
# sd.vars.imp.altmann.sort <- apply(vars.imp.altmann[ , names.vars.sort.altmann], 2, sd)
# 
# # pdf("3-onlyfarm/Results/ranger_snps_only_farm-only_final_importance.pdf", height=9)
# errbar(x=names.vars.sort.altmann, y=means.vars.imp.altmann.sort, yplus = means.vars.imp.altmann.sort + sd.vars.imp.altmann.sort, yminus = means.vars.imp.altmann.sort- sd.vars.imp.altmann.sort)
# abline(v=0)
# # dev.off()

pvals.imp.altmann <- do.call(rbind, lapply(1:5, function(z) obj.fin[[z]]$importance.altmann[ ,"pvalue" ]))
pvals.imp.altmann

# which are always sign at alpha=0.05
always.sign.alt <- apply(pvals.imp.altmann, 2, function(x) all(x <0.05))

always.sign.alt[names.vars.sort.altmann]

# pdf("../../../Results/2-validate-innsbruck/RF_final_model_fh_conf_env_all_variableimportance_plus_sign.pdf", height=9)
par(mar = par('mar') + c(0,5,0,0))
barplot(means.vars.imp.altmann[names.vars.sort.altmann], las=2, horiz = T, col= always.sign.alt[names.vars.sort.altmann]+1,
        main = "Importance of final model (RF) on all children trained on fh + conf + env")
legend("bottomright", fill= c("red", "black"),legend= c("Sign. in all Imputations", "Insign. at least once (alpha=0.05)"))
# dev.off()

# pdf("../../../Results/2-validate-innsbruck/RF_final_model_fh_conf_env_all_variableimportance_only_sign.pdf", height=7)
par(mar = par('mar') + c(0,5,0,0))
barplot(means.vars.imp.altmann[names.vars.sort.altmann][always.sign.alt[names.vars.sort.altmann]] , las=2, horiz = T, col=2, 
        main = "Importance of sign. variables of final model (RF) on all children trained on fh + conf + env  (alpha=0.05)")
# dev.off()

# average pvalues
mean.pval.alt <- apply(pvals.imp.altmann, 2, function(x) mean(x))

# pdf("../../../Results/2-validate-innsbruck/RF_final_model_fh_conf_env_all_variableimportance_averaged_pvals.pdf", height=9)
par(mar = par('mar') + c(5,0,0,0))
barplot(mean.pval.alt[names.vars.sort.altmann], las=2, main = "Altmann-Importance-p-values averaged over imputations")
abline(h=.05, col=2)
# dev.off()