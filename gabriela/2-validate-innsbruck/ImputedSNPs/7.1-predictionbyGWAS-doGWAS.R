### GWAS on CV folds

# we need info per imputation and per fold about
# - SNP
# - it's p-value
# - number of predictors (always the same)






### to be careful about here:

# each time another env block is used as env.only --> use another 'both' so as not to overwrite sth (or solve this in any way)
# for combination approach:
### just cbind the env in order: fam, env.new, interactions 
### bild list with the right.names and give it to function as argument

### lit snps
setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/3-wrp-learn-on-snp-env.R")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/0-functions.R")


library(parallel)



# load imputed data sets:
load(paste0("Data/5imps.rda"))

### data
load("LiteratureSNPs/Data/dat.lit.snps.19.RData")



# load("ImputedSNPs/Data/mat.imp.pruned.rda")

load("Data/gabriel.1707.Rdata")
load("Data/gabriel.23vars.gen.Rdata")
load("Data/dat.outcome.rda")
# load("../../Data/x.env.imp.rda")
load("Data/x.wheeze_onset.rda")

# new: leave out a final test set
# load("Data/fin.ind.rda")

#### MARKNEW
# take only farm and no innsbruck
load("Data/ind.part.rda")
train.fin.ind <- setdiff(1:1707, ind.part$all.inn)



# train.fin.ind <- fin.ind$train.fin.ind

# we use only imputation set 1: since only confounder values may (slightly) vary between the imputation data sets
i=1

  eval(parse(text=paste0("x.env.imp <- imp",i)))
  
  #### MARKNEW
  ## remove variables farm and center for the corresponding subgroupanalysis
  x.env.imp$center <- NULL
  x.env.imp$farm <- NULL
  
  # MARKNEW
  x.env <- model.matrix(~. #+ farm:Sex_female 
                        , x.env.imp)[ train.fin.ind,]
  
  # outcome:
  
  
  ################# adjust snp.data, env.data and outcome (same number of observations etc)
  # snp.data  <- X[ train.fin.ind,]
  # snp.data for interactions
  # snp.for.int <- dat.lit.snps.19[train.fin.ind, -1]
  
  
  ######### change: we define several data sets replacing x.env
  
  ### confounders: always to include - center, age, month of birth
  cnames.env <- colnames(x.env)
  confs <- cnames.env[c(grep("center", cnames.env),  #grep("month", cnames.env), (let month of birth stay environmental)
                        which(cnames.env=="age"), 
                        which(cnames.env=="child_age"), which(cnames.env=="child_BMI"), which(cnames.env=="Sex_female2"))]
  
  x.conf <- x.env[ , confs]
  pos.out <-  which(! gabriel.23vars.gen$serial_gen%in% gabriel.1707$serial_gen )
  
  
  dat.rest <- data.frame(data.frame(y=gabriel.1707$dd_asthma,
                    strata_dna = gabriel.23vars.gen$strata_dna[-pos.out],
                    weight_dna=gabriel.1707$weight_dna)[train.fin.ind, ], x.conf)
  
  
  
  # function for computing pvalues
  comp.snps.pvals <- function(k, snp.data, dat.rest){
    snp <- snp.data[ ,k, drop=FALSE, with=FALSE]
    snp.name <- colnames(snp)
    dat <- data.frame(dat.rest,snp)
    # for survey model
    des <- svydesign(ids=~1, strata = ~strata_dna, weights = ~weight_dna, data=dat, fcp=NULL)
    form <- as.formula(paste("y~", paste(c(confs, snp.name ), collapse = "+")))
    suppressWarnings(svyglm1 <- svyglm(form, design = des, family="binomial", data= dat))
    coef.mat <- summary(svyglm1)$coefficients
    snp.pval <- coef.mat[ snp.name , grep("Pr", colnames(coef.mat))]
    names(snp.pval) <- snp.name
    snp.pval
    # if(k %% 1000 == 0) print(paste("SNP", k, "done!" ))
  }
 
# do it:
require(survey)
library(data.table)
setwd("ImputedSNPs/Data/chromosomes")
l.snp.pvals <- list()
for(j in 1:5) l.snp.pvals[[j]] <- list()
for(i in 1:22){
    load(paste0("ldpruned/dt.thr95ldpruned_chr", i ,".rda"))
    snp.data <- dt.chri.ldpruned[ -pos.out, ][ train.fin.ind, ]
    
    for(j in 1:5){
      # NEWMARK: inner folds for noninnsbruck!!
      train.train.index <- setdiff( 1:length(unlist(ind.part$cvfolds.without.inn$all)), ind.part$cvfolds.without.inn$all[[j]])
      l.snp.pvals[[j]][[i]] <- unlist(mclapply(1:ncol(snp.data), function(x) comp.snps.pvals(x, snp.data = snp.data[ train.train.index, ], dat.rest =  dat.rest[ train.train.index, ] ),
                            mc.cores = 60 ))
      print(paste("CV-iteration", j, "done"))
      
    }
  save(l.snp.pvals, file=paste0("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/2-validate-innsbruck/ImputedSNPs/Data/l.snp.pvals_till_chr",i , ".rda"))
   print(paste("Chromosome", i, "done")) 
}

#### MARKNEW
## change path!

setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma")
save(l.snp.pvals, file = "2-validate-innsbruck/ImputedSNPs/Data/l.snp.pvals.rda")  
 
  
  # remember:
  # any(fin.ind$test.fin.ind %in% fin.ind$train.fin.ind[unlist(fin.ind$ind.within$ind.train)]) # FALSE (as it should)
  # any(fin.ind$train.fin.ind[fin.ind$ind.within$ind.test[[1]] ] %in% fin.ind$train.fin.ind[fin.ind$ind.within$ind.train[[1]] ]) # FALSE (as it should)
  # any(fin.ind$train.fin.ind[fin.ind$ind.within$ind.test[[2]] ] %in% fin.ind$train.fin.ind[fin.ind$ind.within$ind.train[[1]] ]) # TRUE (as it should)  
  # 
  
  
  


# 

