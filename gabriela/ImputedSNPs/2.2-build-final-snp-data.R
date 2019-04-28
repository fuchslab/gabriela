###### 

# the code

## - combines the single chromosome data sets after LD pruning
## - removes the observation not existent for environment data set (from 1708 to 1707 observations)

library(data.table)
library(parallel)

setwd("~/1-Projects/2-Asthma/ImputedSNPs/Data/chromosomes")
snplist <- list()
save.dt.ldpruned <- function(i){
  load(paste0("ldpruned/thr95ldpruned_snplist_chr", i ,".rda"))
  load(paste0("dt.chr", i, ".clean.by.maf.rsq.rda"))
  dt.chri.ldpruned <- dt.chri.clean.by.maf.rsq[ , unlist(snpset), with=F]
  save(dt.chri.ldpruned, file=paste0("ldpruned/dt.thr95ldpruned_chr", i ,".rda"))
}
mclapply(1:22, save.dt.ldpruned, mc.cores = 5)




######## cbind them together
library(data.table)
setwd("~/1-Projects/2-Asthma/ImputedSNPs/Data/chromosomes")
load(paste0("ldpruned/dt.thr95ldpruned_chr", 1 ,".rda"))
dt.imp.pruned <- dt.chri.ldpruned
for(i in 2:22){
  load(paste0("ldpruned/dt.thr95ldpruned_chr", i ,".rda"))
  dt.imp.pruned <- cbind(dt.imp.pruned, dt.chri.ldpruned)
}

####### delete uncomplete observation
load("../../../Data/gabriel.23vars.gen.Rdata")
load("../../../Data/gabriel.1707.Rdata")

out.pos <- which(!gabriel.23vars.gen$serial_gen %in% gabriel.1707$serial_gen)

dt.imp.pruned <- dt.imp.pruned[-out.pos]

save(dt.imp.pruned, file="../dt.imp.pruned.rda")