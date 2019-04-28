### input: new data from pasture study
### output: data set(s) with variables from pasture study where modelfit from GABRIELA is applicable

setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma")
 
 # data set 1:
 load("Data/pasture_core_29vars_20170901.Rdata")
 colnames(pasture_core)
 
# data set 2:
load("Data/pasture_2017-09-22.Rdata")
colnames(pasture)

# some variables of core variables in extended selection, but not all: 
which(colnames(pasture_core) %in% colnames(pasture))


##### load original variables for comparsion/adjustment
load(paste0("Data/5imps.rda"))

colnames(imp1)

### create new data set:
table(pasture$mom_smoke_prg)
table(imp1$SmokingPregn) # check


# recode atopy-variable 
mat.fh.atopy <- cbind(pasture$sib_atopy_ev_prg , pasture$mom_atopy_ev_prg, pasture$dad_atopy_ev_m02)
pasture.FHx_Atopy <- apply(mat.fh.atopy, 1, function(x) as.numeric(any(x ==1, na.rm = F) ))

# recode School parents (education)
School_Parents <- pasture$mom_edu_prg 
School_Parents[School_Parents <= 2] <- 0
School_Parents[School_Parents > 2] <- 1
table(School_Parents)


## build cow straw variables
Contact_Cow_12M <- pasture$contactcow_rec_6y
Contact_Cow_0123 <- sign(pasture$contactcow_rec_2m + pasture$contactcow_rec_1y + pasture$contactcow_rec_2y + pasture$contactcow_rec_3y) 
Forage_Straw_12M <- as.factor(pasture$contactstraw_rec_4y)
Forage_Straw_0123 <- as.factor(sign( pasture$contactstraw_rec_1y + pasture$contactstraw_rec_2m+  pasture$contactstraw_rec_2y +  pasture$contactstraw_rec_3y))


Cow_Straw_12M  <- rep(NA, nrow(pasture)) 
Cow_Straw_12M[Contact_Cow_12M  == 0 & Forage_Straw_12M ==0] <- 0 # continue here  
Cow_Straw_12M[Contact_Cow_12M  == 0 & Forage_Straw_12M ==1] <- 1
Cow_Straw_12M[Contact_Cow_12M  == 1 & Forage_Straw_12M ==0] <- 2
Cow_Straw_12M[Contact_Cow_12M  == 1& Forage_Straw_12M ==1] <- 3



Cow_Straw_0123  <- rep(NA, nrow(pasture)) 
Cow_Straw_0123[Contact_Cow_0123 == 0 & Forage_Straw_0123 ==0] <- 0
Cow_Straw_0123[Contact_Cow_0123 == 0 & Forage_Straw_0123 ==1] <- 1
Cow_Straw_0123[Contact_Cow_0123 == 1 & Forage_Straw_0123 ==0] <- 2
Cow_Straw_0123[Contact_Cow_0123 == 1& Forage_Straw_0123 ==1] <- 3





## mark with "check" if variable/coding is matched 
pasture.new <- data.frame(farm=as.factor(pasture$farmer), #check
                          child_b_month=as.factor(pasture$birth_month), # check
                          child_age=pasture$age_y06/12, # check
                          child_BMI=pasture$bmi_y06, # check
                          FHx_Atopy=as.factor(pasture.FHx_Atopy), #  check
                          Num_Sibs_12=as.numeric(pasture$siblings_prg>=2), # check
                          School_Parents= as.factor(School_Parents), # check (fyi: distr. very different here)
                          Sex_female=as.factor(pasture$sex_female), # check
                          SmokingPregn=as.factor(pasture$mom_smoke_prg), # check
                          AnyFarmMilk_12M = as.factor(pasture$farm_milk_y06),  # check
                          AnyFarmMilk_0123=as.factor(sign(pasture$farm_milk_m12 + pasture$farm_milk_m24+pasture$farm_milk_m36)), # check
                          Contact_Cow_12M = as.factor(Contact_Cow_12M), # check
                          Contact_Cow_0123 = as.factor(Contact_Cow_0123), # check
                          Forage_Straw_12M= as.factor(Forage_Straw_12M), # no 6y --> 4y # check
                          Forage_Straw_0123=as.factor( Forage_Straw_0123), # check
                          Combi_Hay_0123   =as.factor( sign(pasture$regularhay_1y +  pasture$regularhay_18m+ pasture$regularhay_2y +  pasture$regularhay_3y)), # check
                          Combi_Hay_12M =as.factor( sign(pasture$regularhay_6y) ), # check
                          Cow_Straw_0123  =as.factor( Cow_Straw_0123  ), #check
                          Cow_Straw_12M =as.factor(Cow_Straw_12M), # check
                          Barn_12M   = as.factor( pasture$barn_6ym ),      # 12 wasn't given (mit Ege abgesprochen) # check
                          Barn_0123 =as.factor( sign(pasture$barn_18mm + pasture$barn_2ym + pasture$barn_3ym) ),  # check
                          Stable_12M  =as.factor(  pasture$stable_6y ), # check
                          Traffic =factor( rep(1, nrow(pasture)), levels = c(0,1,2,3)), # Ege: not in pasture/ we impute later
                          Mould_12M =as.factor( pasture$mould_y06),             # check                       
                          Dog_01  =as.factor( pasture$dog_house_m12),    # check
                          Cat_01   =as.factor( pasture$cat_house_m12), # check                        
                          Dog_23  =as.factor( sign(pasture$dog_house_m24 )), # missing: third year    # check      
                          Cat_23  =as.factor( sign(pasture$cat_house_m24 )),  # missing: third year        # check
                          Dog_45 =as.factor( sign(pasture$dog_house_ty_m48+pasture$dog_y05)),   # check
                          Cat_45 =as.factor( sign(pasture$cat_house_ty_m48 + pasture$cat_y05)), # check
                          Dog_12M =as.factor(       pasture$dog_house_y06 ), # check
                          Cat_12M =as.factor( pasture$cat_house_y06),        # check
                          Children = pasture$children_prg,           # check
                          Adults = pasture$adults_prg, # check
                          Par_Smoke_Ever =as.factor( sign(pasture$mom_smoke_ev_prg+ pasture$dad_smoke_ev_m02)),# check
                          Par_Smoke_Curr   =as.factor( sign(pasture$mom_smoke_m02+pasture$dad_smoke_m02 )),# check
                          daycare   =as.factor( sign(pasture$daycare_first + pasture$daycare_m12 + pasture$daycare_m24 + pasture$daycare_m36 + pasture$daycare_m48 + pasture$daycare_y05)),   # check     
                          preg_antibiot   =as.factor( pasture$mom_antibiot_prg ),# check
                          stable_cattle_m12 =as.factor( pasture$cattle_rec_6y), ## check
                          stable_cattle_1   =as.factor( sign(pasture$cattle_rec_2m + pasture$cattle_rec_1y)), # check
                          farmmilk_1  =as.factor( pasture$farm_milk_m12), # check
                          fhasthma =as.factor( sign(pasture$dad_asthma_ev_m02 + pasture$mom_asthma_ev_prg + pasture$sib_asthma_ev_prg)),          # check
                          fhhayfev  =as.factor( sign(pasture$dad_hayfever_ev_m02 + pasture$mom_hayfever_ev_prg + pasture$sib_hayfever_ev_prg)),       # check
                          fheczema =as.factor( sign(pasture$dad_eczema_ev_m02 + pasture$mom_eczema_ev_prg + pasture$sib_eczema_ev_prg))# check
                            )

### compare classes/str etc
summary(imp1[,-1])
summary(pasture.new)

save(pasture.new, file="5-validate-pasture/Data/pasture.new.rda")

################################### imputation of missing values #####################################################
# Multivariate Imputation by Chained Equations (MICE)
library(mice)
set.seed(23432)
imputation <- mice(pasture.new, m = 5, maxit = 10, print=FALSE)
past.imp1 <- complete(imputation,1)
past.imp2 <- complete(imputation,2)
past.imp3 <- complete(imputation,3)
past.imp4 <- complete(imputation,4)
past.imp5 <- complete(imputation,5)
save(past.imp1, past.imp2, past.imp3,past.imp4, past.imp5, file="5-validate-pasture/Data/5past.imps.rda")
summary(past.imp1)


