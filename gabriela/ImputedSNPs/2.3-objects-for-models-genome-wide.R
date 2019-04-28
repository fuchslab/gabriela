library(data.table)
library(parallel)

setwd("~/1-Projects/2-Asthma/ImputedSNPs/Data")

load("dt.imp.pruned.rda")
load("../../Data/gabriel.1707.Rdata")

##dim(dt.imp.pruned)
# [1]   1707 744908

# object.size(dt.imp.pruned)/10^9
# 10.303479008 bytes

library(glmnet)

#### try to make it computationally more efficient:
## to be checked if
# - standardizing up front makes glmnet faster
# - or if spare matrix but within lasso fit standardizing is faster.... 

### we use 60 GB RAM

X <- as.matrix(dt.imp.pruned)
save(X, file="mat.imp.pruned.rda")

# standardize before fitting
# dat.imp.pruned.scaled <- scale(dt.imp.pruned, center = T, scale = T)
# save(dat.imp.pruned.scaled, file="dat.imp.pruned.scaled.rda")

# create sparse matrix 
# mat.imp.pruned.sparse <- sparse.model.matrix(~., dt.imp.pruned)
# save(mat.imp.pruned.sparse, file="mat.imp.pruned.sparse.rda")


## Cria a matriz de entrada usando a primeira coluna
data <- dt.imp.pruned

X <- sparse.model.matrix(~data[,1]-1) # minus 1 for no intercept
for (i in 2:ncol(data)) {
    coluna <- sparse.model.matrix(~data[,i]-1)
    X <- cBind(X, coluna)
    if(i %% 10000==0) print(paste("SNP", i))
}
##### subsample into training data sets
# load("../../Data/fin.ind.rda")
# 
# y <- as.factor(gabriel.1707$dd_asthma)

