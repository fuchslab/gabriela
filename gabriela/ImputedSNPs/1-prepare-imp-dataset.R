
setwd("~/1-Projects/2-Asthma/ImputedSNPs/")
library(reshape2)
library(data.table)
library(parallel)
 
# use selection of SNPs from post-QC:
load("Data/mafandrsq.clean.table.rda")
snp.names <- row.names(mafandrsq.clean.table)
rm(mafandrsq.clean.table)


# load estimated dosages
# load("Data/dat.imp2.rda")

################################# clean by MAF and rsq #########################


# dat.imp.clean <-  subset(dat.imp2, select=snp.names)
# save(dat.imp.clean, file="dat.imp.clean.maf.rsq.rda")





# for(i in 1:22){
#     load(paste0("Data/chromosomes/dt.chr",i ,".rda" ))
#     count=1
#     snps.in.chr <- colnames(dt.chri)
# 
#     snp.names.clean.by.maf.rsq<- snp.names[snp.names %in% snps.in.chr]
#     dt.chri.clean.by.maf.rsq <- dt.chri[ , snp.names.clean.by.maf.rsq, with=FALSE]
#    
#     save(dt.chri.clean.by.maf.rsq, file=paste0("Data/chromosomes/dt.chr",i ,".clean.by.maf.rsq.rda"))
#     print(paste("Chromosome", i , "done"))
# }
# 

# done and saved


### Note - computational issues:
  # removing a spec



do.per.chr2 <- function(i){
  load(paste0("Data/chromosomes/dt.chr",i ,".clean.by.maf.rsq.rda" ))
  dt.chri.clean.by.maf.rsq[ , serial_gen:=NULL]
  # count=1
  snps.in.chr <- colnames(dt.chri.clean.by.maf.rsq)
  for(snp in snps.in.chr){
    if( !(snp %in% snp.names)){ # check if to delete because of MAF or call rate
      dt.chri.clean.by.maf.rsq[ ,eval(parse(text=snp)):=NULL]
      print(paste0("snp ", snp, "removed because of MAF or call rate"))
    } else{ # otherwise check if to delete because of correlation
      snp.pos <- which(colnames(dt.chri.clean.by.maf.rsq)==snp)
      j=snp.pos+1
    
      
      if(var(dt.chri.clean.by.maf.rsq[ , j, with=FALSE])==0){
        dt.chri.clean.by.maf.rsq <- dt.chri.clean.by.maf.rsq[ ,j:=NULL]
       # no index-increasing needed, since next snp is falling on current index
      }
      while(j <= ncol(dt.chri.clean.by.maf.rsq)){
        snp.next <- colnames(dt.chri.clean.by.maf.rsq)[j]
        if(abs(dt.chri.clean.by.maf.rsq[ , cor(eval(parse(text=snp)) , eval(parse(text=snp.next)) ) ]) > .95){
          dt.chri.clean.by.maf.rsq <- dt.chri.clean.by.maf.rsq[ , eval(parse(text=snp.next)):=NULL]
          j = ncol(dt.chri.clean.by.maf.rsq)+100
        } else j = j+1
        if(j %% 100 == 0)print(paste("comparing snp", j, "of", ncol(dt.chri.clean.by.maf.rsq),"snps" ))
      }
     
    }
    snps.in.chr <- colnames(dt.chri.clean.by.maf.rsq)
   print(paste("SNP", snp, "done. Dim of data is ", ncol(dt.chri.clean.by.maf.rsq),"(chrom.", i , ")"))
#     count <- count + 1
#     if(count %% 10 == 0) print(paste("iteration",count, "of chrom.", i , "done!"))
  }
  dt.chr.i.clean.by.cor.maf.rsq <- dt.chri.clean.by.maf.rsq
  save(dt.chr.i.clean.by.cor.maf.rsq , file=paste0("Data/chromosomes/dt.chr",i ,".clean.by.cor.maf.rsq.rda"))
}

parallel::mclapply(1:22, FUN = do.per.chr2, mc.cores = 3)


# ################################# clean by cors 
# ## to do...
# # cut-off: 0.97
# 
# # use code from:
# # http://www.r-bloggers.com/introduction-to-feature-selection-for-bioinformaticians-using-r-correlation-matrix-filters-pca-backward-selection/
# # and package 
# 
# dat.imp.clean.nocor
# 
# ################################# create sparse matrix
# 
# ## to do...
# 
# ## first: round values <0.05 to 0 (to make it even sparser)
# 
# # 
# Mat.imp.clean.nocor <- Matrix(..., sparse=TRUE)
#  
#  
# #  # 
# # # dat.imp2 <- dat.imp2[1:100, ]
# # # save(dat.imp2, file="dat.imp2.50obs.rda")
# # load("dat.imp2.50obs.rda")
# # 
# # load("gabriel.1707.Rdata")
# #   load("gabriel.23vars.gen.Rdata")
# # 
# # 
# 
# 
# gabriel.1707 <- gabriel.1707[idx.train, ]
# gabriel.23vars.gen <- gabriel.23vars.gen[idx.train,]
# 
# 
# strata <- gabriel.23vars.gen$strata_dna
# # perform tests (univariate survey models and save p-values)
