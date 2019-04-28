setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/Data/")
load("gabriel.1707.Rdata")
load("gabriel.1708.wheeze.lufu.Rdata")
load("gabriel.23vars.gen.Rdata")
load("gabriel.1708.sta.fmi.care.antib.Rdata")
load("../LiteratureSNPs/Data/dat.lit.snps.19.RData")


###### leave out innsbruck

## list index-partitioning
ind.part <- list()

# innsbruck
ind.part$all.inn <- which(gabriel.1707$center==1)

# farm
ind.part$farm <- which(gabriel.1707$farm==1)

# farm innsbruck
ind.part$farm.inn <- which(gabriel.1707$farm==1 & gabriel.1707$center==1)

# non-farm 
ind.part$nonfarm <- which(gabriel.1707$farm==0)

# non-farm innsbruck
ind.part$nonfarm.inn <- which(gabriel.1707$farm==0 & gabriel.1707$center==1)



### for inner CV

# 5-Fold-CV; in fin.train set: (exactly same code and seed as in 3-wrp-learn-on-snp-env.R)

partition.cv <- function(indices, k=5){
  n <- length(indices)
  set.seed(3534)
  ind <- sample(1:n)
  
  ind.test <- list()
  ind.train <- list()
  for(i in 1:k){
    if(i < k ) ind.test[[i]] <- ind[(1+(i-1)*floor(n/k)):(i*(floor(n/k))) ]
    else ind.test[[i]] <- ind[(1+(i-1)*floor(n/k)):n]
    ind.train[[i]] <- setdiff(ind, ind.test[[i]])
  }
  ind.test
}

ind.part$cvfolds.without.inn$all <- partition.cv(setdiff( 1:1707, ind.part$all.inn))
ind.part$cvfolds.without.inn$farm <-  partition.cv(setdiff(ind.part$farm, ind.part$farm.inn))
ind.part$cvfolds.without.inn$nonfarm <-  partition.cv(setdiff(ind.part$nonfarm, ind.part$nonfarm.inn))
ind.part$cvfolds.without.inn$entire_farm_pop <- partition.cv(ind.part$farm)

save(ind.part, file="ind.part.rda")
