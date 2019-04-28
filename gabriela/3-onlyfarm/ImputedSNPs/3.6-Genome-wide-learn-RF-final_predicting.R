

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
train.fin.ind <- setdiff(ind.part$farm, ind.part$farm.inn)




i=1

eval(parse(text=paste0("x.env.imp <- imp",i)))

# NEWMARK
x.env.imp$farm <- NULL
x.env.imp$center <- NULL

# NEWMARK
x.env <- model.matrix(~. # + farm:Sex_female
                      ,x.env.imp)[,-1]



################# adjust snp.data, env.data and outcome (same number of observations etc)

### for variable importance we need candidate snps in there as well, for comparison (not all in anymore after ld pruning which was no 
# difference for prediction )

lit.snp.names <- colnames(dat.lit.snps.19)[-1]
gen.snp.names <- colnames(dt.imp.pruned) 
lit.not.in.gen.snp.names <- lit.snp.names[! lit.snp.names %in% gen.snp.names]

snp.data  <-  cbind(dt.imp.pruned, dat.lit.snps.19[ , lit.not.in.gen.snp.names])[ train.fin.ind,]

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

x.snp <- snp.data


if(env.if=="fh_conf_snps"){
  x.all = cbind(get.env(env.if)[ train.fin.ind,], x.snp)
  pflist = list(c(1,1,1),  # first block: confounders and  weight always set to 1 
                c(2,1,1), c(1,1,2), # weight one block more
                c(2,1,2), # weight two blocks more
                c(3,1,1), c(1,1,3), # weight one block three times
                c(3,1,3) # weight two blocks three times
  )
}

################# set-up for algorithms
dat<-data.frame(y=factor(y)[ train.fin.ind], x.all)



################ compare this to ranger
ntrees=20000
set.seed(23423)




train.ranger.predict <- ranger::ranger(dependent.variable.name = "y", data = dat, write.forest = T, 
                       save.memory = TRUE, probability = TRUE, num.trees =20000, num.threads = NULL, verbose = FALSE
) 
# importance_janitza <- ranger::importance_pvalues(train.ranger, method = "janitza", dependent.variable.name = "y", data = dat) 


test.fin.ind <- ind.part$farm.inn
x.snp.test <- cbind(dt.imp.pruned, dat.lit.snps.19[ , lit.not.in.gen.snp.names])[test.fin.ind, ]
x.all.test = cbind(get.env(env.if)[ test.fin.ind,], x.snp.test)
dat.test <-data.frame(x.all)
pred.rf <- predict(train.ranger.predict,  data= dat.test ) $predictions[ ,2]

save(train.ranger.predict, pred.rf, file=paste0("3-onlyfarm/ImputedSNPs/Results_genomewide/finalmodel/ranger_predictobject_",env.if,"_imp",i,"_final.rda"))

rm(train.ranger.predict, pred.rf)


######################## same for snps
env.if = "snps_only"
dat.snp <-data.frame(y=factor(y)[ train.fin.ind], x.snp)

ntrees=20000
set.seed(23423)

train.ranger.predict <- ranger::ranger(dependent.variable.name = "y", data = dat.snp, write.forest = T, 
                               save.memory = TRUE, probability = TRUE, num.trees =20000, num.threads = NULL, verbose = FALSE
) 
# importance_janitza <- ranger::importance_pvalues(train.ranger, method = "janitza", dependent.variable.name = "y", data = dat) 



dat.snps.test <- data.frame(x.snp.test)
pred.rf <- predict(train.ranger.predict,  data= dat.snps.test ) $predictions[ ,2]


save(train.ranger.predict, pred.rf,  file=paste0("3-onlyfarm/ImputedSNPs/Results_genomewide/finalmodel/ranger_predictobject_",env.if,"_imp",i,"_final.rda"))





