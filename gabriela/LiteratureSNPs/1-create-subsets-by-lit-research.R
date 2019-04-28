setwd("~/Asthma-Auswertung/Data/")


# for all data sets, delete observation which is not given for all env. variables:
load("~/1-Projects/2-Asthma/Data/gabriel.1707.Rdata")
load("~/1-Projects/2-Asthma/Data/gabriel.23vars.gen.Rdata")

ser.1707 <- as.character(gabriel.1707$serial_gen)
ser.1708 <- as.character(gabriel.23vars.gen$serial_gen)

ser.out <- ser.1708[!(ser.1708 %in% ser.1707)]

# imputed data; only on server available!!
# load("Genetic imputed/dat.imp.RData")

# originally genotyped data, only those occuring in genotyped data filtered from 
# imputed data (dat.gen.imp)
# load("Genetic/dat.gen.or.red.RData")

# load classical variables (with outcome)


################################## Choices of Risk SNPs #############################################
# 1. SNPs found by literature research, wide and "tolerant"; any SNPs associated with asthma; 
#     tried to exclude only-adult-asthma publications 
# 2. SNPs by GWAS-catalog: category "childhood asthma"
# 3. SNPs by GWAS-catalog: where risk allele is known/given/displayed
# 4. Comparison of SNPs detected in Ege 2011

##################################################################################################################################
######################## 1. snps by update wide and "tolerant" literature search #################################################
library(xlsx)
risk.snp.info <- read.xlsx("~/Asthma-Auswertung/Data/Risk-SNPs/MyGWASS-asthma-child.xlsx", 
                           sheetIndex = 1, colIndex = c(3,21, 22))
# get list of risk snps uniquely
risk.snps <- unique( gsub(" ", replacement="", as.character(risk.snp.info$SNPs)) )

# risk snps in imputed data
pos.risk.snps.imp <- which(colnames(dat.imp) %in% risk.snps)

# get new data set containing only risk snps
dat.risk <- dat.imp[ ,pos.risk.snps.imp ]

dat.lit.snps.328 <- data.frame(serial_gen=dat.imp$serial_gen, dat.risk)[-which(dat.imp$serial_gen == ser.out), ]

save(dat.lit.snps.328, file="~/1-Projects//2-Asthma//LiteratureSNPs//Data/dat.lit.snps.328.RData")

##################### Create "correlation-pruned" version ########################################################################

##
x.cor.pruned.9 <- rem.cor(dat.lit.snps.328[,-1], keep.snps=NULL, t=.9)
x.cor.pruned.7 <- rem.cor(dat.lit.snps.328[,-1], keep.snps=NULL, t=.7)
dim(x.cor.pruned.9)
dim(x.cor.pruned.7) # doesn't change to much
dat.lit.snps.cor.pruned.9 <- data.frame(serial_gen=dat.lit.snps.328$serial_gen, x.cor.pruned.9)
save(dat.lit.snps.cor.pruned.9, file="~/1-Projects//2-Asthma//LiteratureSNPs//Data/dat.lit.snps.cor.pruned.9.RData")
##################################################################################################################################
######################## 2. SNPs by GWAS-catalog: category "childhood asthma" ####################################################
#... optionally
setwd("~/1-Projects/2-Asthma/LiteratureSNPs/")
load("Data/dat.lit.snps.328.RData")
library(xlsx)
dat.snps.31 <- read.xlsx("~/1-Projects/2-Asthma/LiteratureSNPs/Data/SNPsGWAScatChildhoodAsthma.xlsx", 
                         header=FALSE, sheetIndex=1, stringsAsFactors=FALSE)

snps.31 <- gsub("([0-9]+).*$", "\\1", dat.snps.31[,1])


length(snps.31)

which(!(snps.31 %in% colnames(dat.lit.snps.328))) # see below: they are indeed not in snps.all (all imputed snps)
# so only 19 snps here;

dat.lit.snps.19 <- dat.lit.snps.328[, c("serial_gen", snps.31[snps.31 %in% colnames(dat.lit.snps.328)])]
save(dat.lit.snps.19, file="~/1-Projects//2-Asthma//LiteratureSNPs/Data/dat.lit.snps.19.RData")

load("~/Asthma2/Data/dat.imp.info.rda")
snps.all <- dat.imp.info$SNP
length(snps.all)
cand.snp.pos.not.in.dat.imp <- which(!(snps.31 %in% snps.all))



##################################################################################################################################
################### 3. SNPs by GWAS-catalog: where risk allele is known/given/displayed ##########################################
# from GWAS catalog
# https://www.ebi.ac.uk/gwas/search?query=childhood%20asthma
# 2015-11-01 at 01:25
setwd("~/1-Projects/2-Asthma/LiteratureSNPs/")
load("Data/dat.lit.snps.328.RData")
# SNP associations to childhood asthma where risk allele known (displayed)
risk.snps <- t(matrix(c("rs7328278","C",
               "rs4658627", "A",
               "rs9815663","T",
               "rs7927044","A", 
               "rs10521233","G",
               "rs6967330","A" ,
               "rs2305480","G",
               "rs928413","G",
               "rs6871536","C",
               "rs3894194","A",
               "rs1295686","T",
               "rs7216389","T"  
               ), 2, 12))
risk.snps

load("../Data/gabriel.23vars.gen.Rdata")

all(risk.snps[,1] %in% colnames(dat.lit.snps.328)) # --> subselection of our 328 SNPs (as it should)

dat.lit.snps.12 <- dat.lit.snps.328[, c("serial_gen", risk.snps[ ,1])]
save(dat.lit.snps.12, file="~/1-Projects/2-Asthma/LiteratureSNPs/Data/dat.lit.snps.12.RData")

load("~/1-Projects/2-Asthma/LiteratureSNPs/Data//dat.lit.snps.12.RData")


## get info of all snps
load("~/Asthma2/Data/dat.imp.info.rda")
snps.all <- dat.imp.info$SNP
length(snps.all)
which(!(snps.31 %in% snps.all)) # ok. (check if we really find the 12 of the 31 childhood asthma SNPs in all the SNPs we have)

dat.imp.info[colnames(dat.lit.snps.12), ]

tabs <- apply(dat.lit.snps.12[ ,-1],2,function(x) table(round(x)))
tabs$rs7328278 <- c(0, tabs$rs7328278)
names(tabs$rs7328278) <- c(0,1,2)
tabs.mat <- sapply(tabs, "[")

cbind(dat.imp.info[colnames(dat.lit.snps.12)[-1], ], t(tabs.mat), risk.snps)

risk.snps
# by the 3 lines above, we'll manually determine which SNPs have to be convertet from 0-1-2 to 0-1-2
# 1 if to convert:
conv <- risk.snps[which(c(1,0,1, 0,1,0,1,1,0,0,1,1)==1),1]
dat.lit.snps.12.risk <- dat.lit.snps.12

### convert values in matrix:
dat.lit.snps.12.risk[  ,conv][round(dat.lit.snps.12.risk[ ,conv])==0] <- 3 
dat.lit.snps.12.risk[  ,conv][round(dat.lit.snps.12.risk[ ,conv])==2] <- 0 
dat.lit.snps.12.risk[  ,conv][round(dat.lit.snps.12.risk[ ,conv])==3] <- 2 

head(dat.lit.snps.12.risk, 20)

save(dat.lit.snps.12.risk, file="~/1-Projects/2-Asthma/LiteratureSNPs/Data//dat.lit.snps.12.risk.RData")



##### first analyses
dat.snp <- dat.lit.snps.12[,-1]

glm1 <- glm(gabriel.23vars.gen$ddasthma ~. , data=dat.snp, family="binomial")
summary(glm1)


dat <- cbind(subset(gabriel.23vars.gen, select=c(ddasthma, weight_dna, strata_dna)), dat.snp)
library(survey)
design <- svydesign(ids=~1, weights=~weight_dna, strata =~strata_dna, data=dat)
svyglm1 <- svyglm(ddasthma ~ rs7328278 + rs4658627+
                    rs9815663 + rs7927044 + rs10521233 + rs6967330 +
                    rs2305480 + rs928413 + rs6871536 + rs3894194 + rs1295686 + rs7216389,
                    design=design, family="quasibinomial", data=dat)
summary(svyglm1)


## training prediction error:
library(pROC)
roc(response = gabriel.23vars.gen$ddasthma, predictor = predict(glm1), direction="<")
roc(response = gabriel.23vars.gen$ddasthma, predictor = predict(svyglm1), direction="<")


# prediction on ind. test ....
# add environment...
# add interactions ...
# add moffat snp ...
# build summary statistics




##################################################################################################################################
################################# 4. Comparison of SNPs detected in Ege 2011 #####################################################

# Ege's significant SNPs from the GABRIEL asthma meta-analysis 
# with farming as a environmental exposure
snps.ege.meta <- c("rs3771166", "rs9273349", "rs1342326", "rs928413", "rs744910", "rs2305480",
                   "rs3894194")
# which are in the genotyped imp. snps
which(snps.ege.meta %in% colnames(dat.imp)) # all but one

### reconstruct Ege 2011:

# Ege's candidate SNPs
snps.ege.can <- c("rs4696480", "rs2569190", "rs2915863", "rs2075817", "rs10759932")
# which are in the genotyped imp. snps
which(snps.ege.can %in% colnames(dat.imp)) # two

# snps by newer meta-anylsis????



