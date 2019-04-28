###### 
## do it like:
# http://corearray.sourceforge.net/tutorials/SNPRelate/#ld-based-snp-pruning

# # installed by:
# source("http://bioconductor.org/biocLite.R")
# biocLite("gdsfmt")
# biocLite("SNPRelate")
library("SNPRelate")
library(parallel)

# 
set.seed(3423)

# Load data
ldprune.per.chr <- function(i){
 
  ## SNPRelate uses integers and not the dosages... but this is ok just for the pruning (approved/suggested by mathias heinig)
  
  setwd("~/1-Projects/2-Asthma/ImputedSNPs/Data/chromosomes")
  
  load(paste0("dt.chr", i, ".clean.by.maf.rsq.rda"))
  
  load("../../../Data/gabriel.23vars.gen.Rdata")
  
  # use estimates and not dosages (just for the ld pruning)
  mat.snps.int <- round(dt.chri.clean.by.maf.rsq)
  
  # mat.snps <- t(as.matrix(dt.chri.clean.by.maf.rsq))
  mat.snps <- t(as.matrix(mat.snps.int))
  
  snp.id <- colnames(dt.chri.clean.by.maf.rsq)
  # Create a gds file
  snpgdsCreateGeno(paste0("test", i, ".gds"), genmat = mat.snps,
                   sample.id = gabriel.23vars.gen$serial_gen , snp.id = snp.id,
                   snp.chromosome = rep(i, length(snp.id)),
                   snpfirstdim=TRUE)
  
  
  # Open the GDS file
  genofile <- snpgdsOpen(paste0("test", i, ".gds"))
  # Close the GDS file
  
  snpset <- snpgdsLDpruning(genofile, ld.threshold=0.95, remove.monosnp = T)
  snpgdsClose(genofile)
  
  save(snpset, file=paste0("ldpruned/thr95ldpruned_snplist_chr", i ,".rda"))
  rm(list = ls())
}
### we had to re-run due to read stream errors: as follows
 mclapply( 1:22, ldprune.per.chr, mc.cores=11)
# mclapply( c(3,4,6,7,9:12,14:22 ), ldprune.per.chr, mc.cores=11)
# mclapply( c(3,4,7,9:12,14:18,20:22 ), ldprune.per.chr, mc.cores=11)
# mclapply( 4, ldprune.per.chr, mc.cores=1)
# mclapply( 7, ldprune.per.chr, mc.cores=1)
# mclapply( 22, ldprune.per.chr, mc.cores=1)
# mclapply( 9, ldprune.per.chr, mc.cores=1)
# mclapply( 10, ldprune.per.chr, mc.cores=1)
# mclapply( 11, ldprune.per.chr, mc.cores=1)
# mclapply( 12, ldprune.per.chr, mc.cores=1)
# mclapply( 14, ldprune.per.chr, mc.cores=1)
# mclapply( 15, ldprune.per.chr, mc.cores=1)
# mclapply( 16, ldprune.per.chr, mc.cores=1)
# mclapply( 17, ldprune.per.chr, mc.cores=1)
# mclapply( 20, ldprune.per.chr, mc.cores=1)
# mclapply( 21, ldprune.per.chr, mc.cores=1)
