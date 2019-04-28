# setwd("storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/Data")
setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/Data/")

load("gabriel.1707.Rdata")
load("gabriel.1708.wheeze.lufu.Rdata")
load("gabriel.23vars.gen.Rdata")
load("gabriel.1708.sta.fmi.care.antib.Rdata")
# we'll include the SNPs, too
load("../LiteratureSNPs/Data/dat.lit.snps.19.RData")

list.of.packages <- c("VIM", "Amelia", "mice", "survey", "MASS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# require(list.of.packages)
library("VIM")
library("Amelia")
library("mice")
library(survey)
library(MASS)
#######################################
################################# data preparation


dat.env.1707 <- merge(gabriel.1707, gabriel.1708.sta.fmi.care.antib, by="serial_gen")

##### variables more than 25% values missings:
given.vals.per.var <- apply(dat.env.1707, 2, function(x) length(na.omit(x)))
given.vals.per.var[ given.vals.per.var < 1707*.75] # non of the relevant variables


# check again
all(dat.env.1707$serial_gen==dat.lit.snps.19$serial_gen)


w

head(dat.env.1707)
attach(dat.env.1707)
dat.env.temp <- data.frame(  as.character(serial_gen), 
                              as.factor(strata_12), # strata accounting for ddasthma/notddasthma, farmer/farmexposed/notfarmer, ?
                             weight_dna, # weight for accounting the disproportionate stratified random sampling
                             as.factor(center), # study center
                             as.factor(farm),
                            # any_asthma, 
                             dd_asthma,
                            # as.factor(dd_hay_fever), --> outcome (develops even later than asthma)
                            # atopy035, 
                            # atopy07,
                            # atop_asthma,
                            # non_atop_asthma,
                           # as.factor(curr_wheeze), # wheezing during the last 12 months Gabriel g2m5 --> more an outcome
                            as.factor(dd_eczema), # Has a doctor ever diagnosed neurodermitis (atopic eczema, endogenous eczema, atopic dermatitis)? Gabriel g1r16
                          ### variables fieldwork:                          
                            as.factor(child_b_month),#  Month of birth
                            # child_b_year,
                            child_age,
                            child_height,
                            child_bodymass,
                            child_BMI,
                            fw_ph2_date,
                          ### NEWLY DEFINDED CONFOUNDER VARIABLES
                            as.factor(FHx_Atopy), # potential confounder parental atopy or sibling atopy
                            Num_Sibs_12, # potential confounder >=2 siblings
                            as.factor(School_Parents), # Potential confounder high parental education
                            as.factor(Sex_female), # 
                            as.factor(SmokingPregn),
                          ### NEWLY DEFINDED FARM VARIABLES
                            as.factor(AnyFarmMilk_12M), # Consumption of farm milk past 12 months
                            as.factor(AnyFarmMilk_0123), # Consumption of farm milk pregnancy to age 3yrs
                            as.factor(Contact_Cow_12M), # Contact with cows past 12 months
                            as.factor(Contact_Cow_0123), # Contact with cows pregnancy to age 3yrs
                            as.factor(Forage_Straw_12M), # Contact with straw past 12 months
                            as.factor(Forage_Straw_0123), # Contact with straw pregnancy to age 3yrs
                            as.factor(Combi_Hay_0123), # Contact with hay pregnancy to age 3yrs
                            as.factor(Combi_Hay_12M), # Contact with hay past 12 months
                            as.factor(Cow_Straw_0123), # Contact with cows and/or straw pregnancy to age 3yrs
                            as.factor(Cow_Straw_12M), # Contact with cows and/or straw past 12 months
                            as.factor(Barn_12M), # Stay in barn past 12 months
                            as.factor(Barn_0123), # Stay in barn pregnancy to age 3yrs
                            as.factor(StorageRoom_12M), # Stay in fodder storage room past 12 months
                            as.factor(StorageRoom_0123), # Stay in fodder storage room pregnancy to age 3yrs
                            as.factor(Stable_12M), # stay in the cattle stable independently of the activity of the adults last 12 months once a week
                          ### variables Living Environment
                            as.ordered(Traffic), # On week days, how often do trucks or busses drive by on the street where your child lives?
                            as.factor(Coalheat_12M), # In the last 12 months, did you use wood or coal to heat or cook that had a burning outlet in the living quarters?
                            as.factor(Mould_12M), # In the last 12 months, were there rooms in the house where your child spent time with visible mould, fungus or other humidity damages?
                            as.factor(Animals_Bedroom), # Where domestic animals ever allowed to stay in the room where your child sleeps?
                            as.factor(Dog_01), # Dog allowed to stay in the room where your child sleeps In the first year of life
                            as.factor(Cat_01), # Cat allowed to stay in the room where your child sleeps In the first year of life
                            as.factor(Other_Animal_01), # Other animal allowed to stay in the room where your child sleeps?  In the first year of life
                            as.factor(Dog_23), # Dog allowed to stay in the room where your child sleeps In the 2. and 3. year of life
                            as.factor(Cat_23), # Cat allowed to stay in the room where your child sleeps In the 2. and 3. year of life
                            as.factor(Other_Animal_23), # Other animal allowed to stay in the room where your child sleeps? In the 2. and 3. year of life
                            as.factor(Dog_45), # Dog allowed to stay in the room where your child sleeps In the 4. and 5. year of life
                            as.factor(Cat_45), # Cat allowed to stay in the room where your child sleeps In the 4. and 5. year of life
                            as.factor(Other_Animal_45),# Other animal allowed to stay in the room where your child sleeps? In the 4. and 5. year of life
                            as.factor(Dog_12M), # Dog allowed to stay in the room where your child sleeps  In the last 12 months
                            as.factor(Cat_12M), # Cat allowed to stay in the room where your child sleeps In the last 12 months 
                            as.factor(Other_Animal_12M), # Other animal allowed to stay in the room where your child sleeps? In the last 12 months
                          #  as.factor(Avoid_Any), # Did you change certain things in your living environment due to existing allergies or asthma in your family?
                          #  as.factor(Avoid_Dog),
                          #  as.factor(Avoid_Cat),
                          #  as.factor(Avoid_Other_Animal),
                          #  as.factor(Avoid_Carpet),
                          #  as.factor(Avoid_Eggs),
                          #  as.factor(Avoid_Milk),
                          #  as.factor(Avoid_Nuts),
                          #  as.factor(Avoid_Other_Food),
                          #  as.factor(Avoid_Other), # Did you change certain things in your living environment due to existing allergies or asthma in your family? Other changes
                            Children, # How many people (including children) live currently in your household? Children between 0-18 years
                            Adults, # How many people (including children) live currently in your household? Adults (18 years and older)
                            as.factor(Par_Smoke_Ever), # Do the parents smoke? father OR mother (1 OR 2)
                            as.factor(Par_Smoke_Curr),
                          as.factor(daycare),
                          as.factor(preg_antibiot),  
                          as.factor(stable_cattle_m12),  
                          as.factor(stable_cattle_1),
                           # as.factor(farmmilk_m12), same as AnyFarmMilk_12M (s.above)
                          as.factor(farmmilk_1),
                          stringsAsFactors = FALSE
)     
detach(dat.env.1707)                    

# adjust colnames
colnames(dat.env.temp) <- gsub("as.factor.", "", colnames(dat.env.temp))
colnames(dat.env.temp) <- gsub("as.ordered.", "", colnames(dat.env.temp))
colnames(dat.env.temp) <- gsub("as.character.", "", colnames(dat.env.temp))
colnames(dat.env.temp) <- gsub(".", "", colnames(dat.env.temp), fixed=TRUE)

dim(dat.env.temp)
dim(na.omit(dat.env.temp))



## include only real predictors (in abpsrache mit m.ege) 
colnames(dat.env.temp)

dat.env.fin <- data.frame(subset(dat.env.temp, select=c(serial_gen, strata_12, weight_dna, center, farm, dd_asthma, child_b_month, child_age,
                              child_BMI, FHx_Atopy, Num_Sibs_12, School_Parents, Sex_female, SmokingPregn, 
                              AnyFarmMilk_12M, AnyFarmMilk_0123, Contact_Cow_12M, Contact_Cow_0123, Forage_Straw_12M,
                              Forage_Straw_0123, Combi_Hay_0123, Combi_Hay_12M, Cow_Straw_0123, Cow_Straw_12M,
                              Barn_12M, Barn_0123, Stable_12M, Traffic, Mould_12M, Dog_01, Cat_01, 
                              Dog_23, Cat_23, Dog_45, Cat_45, Dog_12M, Cat_12M, Children, Adults, Par_Smoke_Ever,
                              Par_Smoke_Curr , 
                              # now new ones 
                              daycare,
                              preg_antibiot,
                              stable_cattle_m12,
                              stable_cattle_1,
                              # farmmilk_m12, 
                              farmmilk_1
                              )), dat.lit.snps.19[ ,-1])




# investigation of all variables
dim(dat.env.fin)
summary(dat.env.fin)


# add family anamnse
dat.fam <- subset(gabriel.23vars.gen, select=c(serial_gen, fhasthma, fhhayfev, fheczema))
dat.fam$fhasthma <- as.factor(dat.fam$fhasthma)
dat.fam$fhhayfev <- as.factor(dat.fam$fhhayfev)
dat.fam$fheczema <- as.factor(dat.fam$fheczema)

dat.env.plus.fam <- merge(x = dat.env.fin, y=dat.fam, by="serial_gen")
dim(dat.env.plus.fam)
head(dat.env.plus.fam)
##### merge final env.


save(dat.env.plus.fam, file="dat.env.plus.fam.rda")



#####################################################################################################################
####################################### data preparation ############################################################
#####################################################################################################################
covs.env <- subset(dat.env.plus.fam, select=-c(serial_gen, dd_asthma, strata_12, weight_dna))




###################### correct distribution through transformations #################################################
pdf(file="../Results/environment-preparation/missmap.env.pdf")
missmap(covs.env)
dev.off()
pdf(file="../Results/environment-preparation/summary.variables.pdf")
for(i in 1:ncol(covs.env)){
  barMiss(covs.env, pos = i) # (also) look at missing values
}
dev.off()

### no too extreme skewness in (less) continuous variables --> no transformations necessary
which(apply(covs.env, 2, is.numeric))

### which variables contain missings - table
num.miss <- apply(covs.env, 2, function(z) length(which(is.na(z))) )/nrow(covs.env)
vars.with.miss <- num.miss[ num.miss > 0]



source("../Results/Revision-plots-and-tables/0-help-objects/1-variable-coding-function.R")

names.vars.with.miss <- convert.variable.names(names(vars.with.miss))[,1]
names.vars.with.miss[1] <- "child birthmonth"
names(vars.with.miss) <- names.vars.with.miss
write.csv(vars.with.miss, file="../Results/environment-preparation/vars.with.miss.csv")

num.miss[ num.miss > 0.07]


##################### check for outliers in continous variables ######################################################
str(covs.env)
boxplot(covs.env$child_age) # makes sense --> no deletion
boxplot(covs.env$child_BMI) # makes sense --> no deletion
## all fine.

################################### imputation of missing values #####################################################
# Multivariate Imputation by Chained Equations (MICE)
set.seed(23432)
imputation <- mice(covs.env, m = 5, maxit = 10, print=FALSE)
snp.pos <- grep("rs", colnames(covs.env))
imp1 <- complete(imputation,1)[ , - snp.pos ]
imp2 <- complete(imputation,2)[ , - snp.pos ]
imp3 <- complete(imputation,3)[ , - snp.pos ]
imp4 <- complete(imputation,4)[ , - snp.pos ]
imp5 <- complete(imputation,5)[ , - snp.pos ]




# which percentage of values between imputation data sets differ 
m <- matrix(NA,5,5)
i = 1
j = 1
for(impi in list(imp1,imp2,imp3,imp4,imp5) ){
  for(impj in list(imp1,imp2,imp3,imp4,imp5)){
    prop.unequal.per.var <- rep(NA, ncol(impi))
    for(k in 1:ncol(impi)){
      prop.unequal.per.var[k]<- length(which(impi[ ,k]!=impj[ ,k]))/nrow(impi)
    }
    m[i,j] <- max(prop.unequal.per.var)
    j = j + 1
  }
  j = 1
  i = i +1
} 

colnames(m) <- paste0("imp", 1:5)
row.names(m) <- paste0("imp", 1:5)

 # maximum of values differing is 2% 
## => worst case: these 2% are in one column: 
## --> which is still no problem:
## --> we decide to take 1 set for our analysis

library(gridExtra)
pdf("../Results/environment-preparation/percentage.imputations.non-overlap.pdf", height=11, width=8.5)
grid.table(round(m,2)) 
dev.off()

x.env.imp <- imp1
save(x.env.imp, file="x.env.imp.rda")
save(imp1, imp2, imp3,imp4, imp5, file="5imps.rda")



# ############# additional variable age of onset wheeze:
# 
# plus.wheeze <- merge(dat.env.plus.fam, gabriel.1708.wheeze.lufu, by="serial_gen")
# dat.wheeze_onset <- subset(plus.wheeze, select=c(serial_gen, wheeze_onset))
# dat.wheeze_onset$wheeze_onset <- as.numeric(dat.wheeze_onset$wheeze_onset)
# pdf(file="../Results/environment-preparation/summary.wheeze_onset.pdf")
# barMiss(dat.wheeze_onset,pos=2)
# dev.off()
# 
# x.wheeze_onset <- dat.wheeze_onset[ ,2,drop=FALSE]
# save(x.wheeze_onset, file="x.wheeze_onset.rda")
# 
# ############## separate env and family anamnesis
# 
# 
# head(x.env.imp.no.fam)
# 
# # environment without family anamnesis
# x.env.imp.no.fam <- subset(x.env.imp, select=-c(fhasthma, fhhayfev,fheczema,FHx_Atopy))
# # model.matrix including interaction farm*sex
# x.env.may.new <- model.matrix(~. + farm:Sex_female, data = x.env.imp.no.fam)
# 
# # only family anamnesis
# x.fam <-  subset(x.env.imp, select=c(fhasthma, fhhayfev,fheczema,FHx_Atopy))
# 
# 
# ############# SNP-env. interaction matrix
# load("../LiteratureSNPs/Data/dat.lit.snps.19.RData")
# x.snps <- subset(dat.lit.snps.19, select = -serial_gen)
# x.snps.plus.env <- cbind(x.env.imp, x.snps)
# 
# # build formula
# form.snps <- paste(colnames(x.snps), collapse="+")
# form.rest <- "(Num_Sibs_12 + farm + farm:Sex_female + FHx_Atopy + fhasthma + fhhayfev + fheczema)"
# 
# # = 19 * 8 +1 (19 SNPs, 8 rest variables (two for farm:Sex_female), 1 intercept)
# 
# form.inter <- 
# x.snp.inter <- model.matrix(eval(parse(text=paste("~(", form.snps, "):", form.rest))) , data = x.snps.plus.env)[,-1] # no intercept
# 
# save(x.env.may.new, x.fam, x.snp.inter, x.snps, file="predictors_may_new_19SNPs.rda")
# 

######################### final-test set index and inner CV-indices:
set.seed(23234)
test.fin.ind <- sample(1:1707, 200)
train.fin.ind <- setdiff(1:1707, test.fin.ind)

# 5-Fold-CV; in fin.train set: (exactly same code and seed as in 3-wrp-learn-on-snp-env.R)
k <- 5
n <- length(train.fin.ind)
set.seed(3534)
ind <- sample(1:n)

ind.test <- list()
ind.train <- list()
for(i in 1:k){
  if(i < k ) ind.test[[i]] <- ind[(1+(i-1)*floor(n/k)):(i*(floor(n/k))) ]
  else ind.test[[i]] <- ind[(1+(i-1)*floor(n/k)):n]
  ind.train[[i]] <- setdiff(ind, ind.test[[i]])
}

fin.ind <- list(test.fin.ind=test.fin.ind, train.fin.ind=train.fin.ind, ind.within=list(ind.train=ind.train, ind.test=ind.test)) 
# note ind.within probably not so good to use, when data set is reduced because of missings in outcome....
save(fin.ind, file="fin.ind.rda" )


