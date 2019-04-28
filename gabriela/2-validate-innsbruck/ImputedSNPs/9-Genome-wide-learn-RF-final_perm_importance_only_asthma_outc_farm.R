

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

# take only asthma! and no innsbruck
load("Data/ind.part.rda")
# no innsbruck
ind.noinn <- setdiff(1:1707, ind.part$all.inn)
train.fin.ind <- ind.noinn[gabriel.1707$dd_asthma[ind.noinn] == 1 ]
test.fin.ind <- ind.part$all.inn[gabriel.1707$dd_asthma[ind.part$all.inn]==1 ]


i=1

eval(parse(text=paste0("x.env.imp <- imp",i)))

# NEWMARK
# x.env.imp$farm <- NULL
# x.env.imp$center <- NULL

# NEWMARK
# x.env <- model.matrix(~. # + farm:Sex_female
#                       ,x.env.imp)[,-1]
# 


################# adjust snp.data, env.data and outcome (same number of observations etc)

### for variable importance we need candidate snps in there as well, for comparison (not all in anymore after ld pruning which was no 
# difference for prediction )

lit.snp.names <- colnames(dat.lit.snps.19)[-1]
gen.snp.names <- colnames(dt.imp.pruned) 
lit.not.in.gen.snp.names <- lit.snp.names[! lit.snp.names %in% gen.snp.names]

snp.data  <-  cbind(dt.imp.pruned, dat.lit.snps.19[ , lit.not.in.gen.snp.names])[ train.fin.ind,]

######### change: we define several data sets replacing x.env

### confounders: always to include - center, age, month of birth
# cnames.env <- colnames(x.env)
# confs <- cnames.env[c(grep("center", cnames.env),  #grep("month", cnames.env), (let month of birth stay environmental)
#                       which(cnames.env=="age"), 
#                       which(cnames.env=="child_age"), which(cnames.env=="child_BMI"), which(cnames.env=="Sex_female2"))]
# 
# x.conf <- x.env[ , confs]
# x.env <- x.env[ , -which(cnames.env %in% confs)]


## 1. x.env.nofh (no family anamnese)
# fam.an <- c("fhasthma2", "fhhayfev2", "fheczema2", "FHx_Atopy2" )
# x.env.nofh <-  x.env[ , - which(colnames(x.env) %in% fam.an)]

## 2. x.fh (family anamnese)
# x.fh <- x.env[ , fam.an]

## 3. x.interact (several interactions of snps with fam-hist and env)
# formula:

# NEWMARK
form.interact.char <- paste( "~(", paste(colnames(snp.data[,-1]), collapse="+"), ")", ":", 
                             "(", 
                             paste(c("Num_Sibs_12", 
                                     # "farm2", "farm2:Sex_female2",
                                     fam.an ), collapse="+"),
                             ")")

## old (age of onset wheeze-variable)
# x.wo <- as.matrix(x.wheeze_onset)[train.fin.ind,,drop=FALSE ]



# rm(x.env)

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

# env.if = "fh_conf_snps"
obs.weights=gabriel.1707$weight_dna[ train.fin.ind]
strata = gabriel.1707$strata_12[ train.fin.ind]
######################################################################################################
################################## dda: ##############################################################
######################################################################################################

# y <- get.y(y.label)

# y.which <- colnames(y)
# y <- y[ ,1]


#### note: we don't do this here since our final indices refere to the WHOLE data sets: we subset later s. below    
#     ind.y <- which(!is.na(y)) 
#     y <- y[ind.y]


# x.env <- get.env(env.if)
#### s. 3-wrp-learn-on-snps-env.R   
########## create x.snp

x.snp <- snp.data



  x.all = cbind(center=imp1$center[ train.fin.ind], Sex_female=imp1$Sex_female[ train.fin.ind], x.snp)

################# set-up for algorithms
dat<-data.frame(y=factor(gabriel.1707$farm)[ train.fin.ind], x.all)



################ compare this to ranger
ntrees=20000
set.seed(23423)




train.ranger <- ranger::holdoutRF(dependent.variable.name = "y", data = dat, write.forest = T, 
                       save.memory = TRUE, probability = TRUE, num.trees =ntrees, num.threads = 8, verbose = FALSE
) 
importance_janitza <- ranger::importance_pvalues(train.ranger, method = "janitza", dependent.variable.name = "y", data = dat) 

save(train.ranger, importance_janitza, file=paste0(
  "2-validate-innsbruck/ImputedSNPs/Results_genomewide/ranger_object_genomewide_final_asthma_outcomefarm.rda"))




########################################### prediction
# train ranger for prediction
library(ranger)
train.ranger2 <- ranger(dependent.variable.name = "y", data = dat, write.forest = T, 
                        save.memory = TRUE, probability = TRUE, num.trees =ntrees, num.threads = 10, verbose = FALSE
) 

#
x.snp.inn <-  cbind(dt.imp.pruned, dat.lit.snps.19[ , lit.not.in.gen.snp.names])[ test.fin.ind,]
x.all.inn = cbind(center=imp1$center[ test.fin.ind], Sex_female=imp1$Sex_female[ test.fin.ind], x.snp.inn)
dat.inn <- data.frame(y=factor(gabriel.1707$farm)[ test.fin.ind], x.all.inn)

 preds <- predict(train.ranger2, data=  dat.inn)$predictions[ ,2]
y.inn <- dat.inn$y
 
# validate:
(roc.obj <- roc(y.inn, preds, ci=TRUE, direction="<"))
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/5-prediction-and-validation.R")
# weighted:
ci.auc.bootstrap.w(roc.obj, weights = gabriel.1707$weight_dna[test.fin.ind], boot.n = 1000, parallel = FALSE, conf.level = .95)
# dim(dat.inn) # only 149 observations for validation --> big CIs


save(train.ranger2, y.inn, preds , file=paste0(
   "2-validate-innsbruck/ImputedSNPs/Results_genomewide/ranger_prediction_object_genomewide_final_asthma_outcomefarm.rda"))
 

####### compare significance with univariate models
# no weighting: 
dat.temp <- dat[ , c("y", "center", "Sex_female", "rs2792249")]
glm.temp <- glm(y ~., data = dat.temp, family="binomial")
summary(glm.temp)
# svyglm
library(survey)
dat.temp2 <- na.omit(data.frame(dat.temp, strata, obs.weights))
des <- svydesign(~1, strata = ~strata, weights = ~ obs.weights, data = dat.temp2)
summary(svyglm(y~ center + Sex_female+rs2792249, design=des ))

dat.temp <- dat[ , c("y", "center", "Sex_female", "rs17064697")]
glm.temp <- glm(y ~., data = dat.temp, family="binomial")
summary(glm.temp)


