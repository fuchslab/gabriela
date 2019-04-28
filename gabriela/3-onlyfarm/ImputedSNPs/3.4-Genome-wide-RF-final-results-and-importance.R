setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/3-wrp-learn-on-snp-env.R")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/0-functions.R")


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
test.fin.ind <- ind.part$farm.inn

i=1
##############################################################################################################################################
############################################################ "all predictors ('both')" #######################################################


library(ranger)

env.if="both"
var.imp <- list()
for(i in 1:5){
  load(paste0("3-onlyfarm/ImputedSNPs/Results_genomewide/finalmodel/ranger_object_",env.if,"_imp",i,".rda"))
  var.imp[[i]] <- train.ranger$variable.importance
}

mat.var.imp <- rbind(var.imp[[1]], var.imp[[2]], var.imp[[3]], var.imp[[4]], var.imp[[5]])


mean.var.imp <- colMeans(mat.var.imp)
# sd.var.imp <- apply(mat.var.imp, 2, sd)

n.top <- 50

mean.var.imp.top <- sort(mean.var.imp, decreasing = T)[n.top:1]
names.top <- names(sort(mean.var.imp, decreasing = T)[n.top:1])

## plot best variables as boxplots
library(Hmisc)
sd.var.imp.top <- apply(mat.var.imp[ , names.top], 2, sd)
pdf("3-onlyfarm/Results/ranger_both_farm-only_final_importance.pdf", height=9)
errbar(x=names.top, y=mean.var.imp.top, yplus = mean.var.imp.top + sd.var.imp.top, yminus = mean.var.imp.top- sd.var.imp.top)
dev.off()


##############################################################################################################################################
############################################################ same for "snps_only" ############################################################
env.if="snps_only"
var.imp <- list()
for(i in 1:5){
  load(paste0("3-onlyfarm/ImputedSNPs/Results_genomewide/finalmodel/ranger_object_",env.if,"_imp",i,".rda"))
  var.imp[[i]] <- train.ranger$variable.importance
}

mat.var.imp <- rbind(var.imp[[1]], var.imp[[2]], var.imp[[3]], var.imp[[4]], var.imp[[5]])


mean.var.imp <- colMeans(mat.var.imp)
# sd.var.imp <- apply(mat.var.imp, 2, sd)

n.top <- 50

mean.var.imp.top <- sort(mean.var.imp, decreasing = T)[n.top:1]
names.top <- names(sort(mean.var.imp, decreasing = T)[n.top:1])

## plot best variables as boxplots
library(Hmisc)
sd.var.imp.top <- apply(mat.var.imp[ , names.top], 2, sd)

pdf("3-onlyfarm/Results/ranger_snps_only_farm-only_final_importance.pdf", height=9)
errbar(x=names.top, y=mean.var.imp.top, yplus = mean.var.imp.top + sd.var.imp.top, yminus = mean.var.imp.top- sd.var.imp.top)
dev.off()

### check if best SNPs are correlated with other SNPs/ sex/ fh
dat.snps.top <- dt.imp.pruned[ , names.top, with=T]

dat.test <- dt.imp.pruned[test.fin.ind ,]
#### prediction on left out set
pred <- predict(train.ranger, dat.test)

