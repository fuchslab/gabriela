## ### lit snps
setwd("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/2-Asthma/Data/")
load("gabriel.1707.Rdata")
load("gabriel.1707.asthma-pheno.Rdata")
load("gabriel.1708.asthma.def.vars.Rdata")
load("gabriel.1708.wheeze.inhaler.Rdata")


### 1 variable for asthma phenotypes
asthma.pheno <- rep(NA, nrow(gabriel.1707.asthma))
asthma.pheno[gabriel.1707.asthma$transient_wheeze ==1 ] <- 1
asthma.pheno[gabriel.1707.asthma$persistent_wheeze ==1 ] <- 2
asthma.pheno[gabriel.1707.asthma$late_onset_wheeze ==1 ] <- 3
asthma.pheno[gabriel.1707.asthma$transient_wheeze  ==0 &gabriel.1707.asthma$persistent_wheeze==0 &
               gabriel.1707.asthma$late_onset_wheeze == 0] <- 0






## dd_asthma + hayfever
asthma_hay <- rowSums(subset(gabriel.1707, select=c(dd_hay_fever, dd_asthma)))
table(asthma_hay)
asthma_hay[asthma_hay==2] <- 1

## dd_asthma + eczema 
asthma_eczema <- rowSums(subset(gabriel.1707, select=c(dd_eczema, dd_asthma)))
table(asthma_eczema)
asthma_eczema[asthma_eczema==2] <- 1


## dd_asthma + hayfever + eczema
asthma_hay_eczema <- rowSums(subset(gabriel.1707, select=c(dd_hay_fever, dd_eczema, dd_asthma)))
table(asthma_hay_eczema)
asthma_hay_eczema[asthma_hay_eczema==2 | asthma_hay_eczema==3] <- 1

## dd_asthma vs no asthma + no hayfever + no eczema
asthma_pure <- gabriel.1707$dd_asthma
asthma_pure[gabriel.1707$dd_asthma==0 & (gabriel.1707$dd_hay_fever==1 | gabriel.1707$dd_eczema==1 ) ] <- NA
table(asthma_pure, useNA="always")

## non-atopic asthma vs none
asthma_nonatopic <- gabriel.1707$non_atop_asthma
asthma_nonatopic[gabriel.1707$dd_asthma == 1 & gabriel.1707$non_atop_asthma==0]<- NA
table(asthma_nonatopic, useNA="always")


################## new (relevant) definitions:
## we'll define:
# dda: doctor-diagnosed asthma
# dda_notany: doctor-diagnosed asthma vs not-anyasthma (=> set of setminus(anyasthma, dda) is excluded)
# ddahay_notdda: dda and hayfever vs not doctor-diagnosed asthma (=> setminus(hayfever and dda,dda ) excluded)
# ddaecz_notdda: dda and eczema vs not doctor-diagnosed asthma (=> setminus(eczema and dda,dda ) excluded)

dda_notany <- gabriel.1707$dd_asthma
dda_notany[ !((gabriel.1707$any_asthma == 1 & gabriel.1707$dd_asthma==1) | (gabriel.1707$any_asthma==0)) ] <- NA
# cbind(dda_notany, gabriel.1707$any_asthma, gabriel.1707$dd_asthma)

ddahay_notdda <- gabriel.1707$dd_asthma
ddahay_notdda[ (gabriel.1707$dd_asthma == 1 & gabriel.1707$dd_hay_fever ==0)  ] <- NA

ddaecz_notdda <- gabriel.1707$dd_asthma
ddaecz_notdda[ (gabriel.1707$dd_asthma == 1 & gabriel.1707$dd_eczema ==0)  ] <- NA



################# new outcome: dda==1 OR wheeze12==1 vs. anyasthma==0
# same as dda vs not_any; but take wheeze12(subset of anyasthma) to dda (so we loose less observations)



dat.with.wh12 <- merge(gabriel.1707[ , c("serial_gen", "dd_asthma","any_asthma")], gabriel.1708.wheeze_inhaler, by="serial_gen")[ ,1:4]
head(dat.with.wh12)

dat.with.wh12$dda.wheeze_notany <- NA

dat.with.wh12$dda.wheeze_notany[ dat.with.wh12$dd_asthma==1 | dat.with.wh12$ph1_wheeze_12m==1 ] <- 1 
dat.with.wh12$dda.wheeze_notany[ dat.with.wh12$any_asthma==0] <- 0 

#### asthma without bronchitis
### note: data set here is exceptionally not ordered by serial gen!!!!!!!!!!!!!!!
# load("gabriel.23vars.gen.Rdata")
# which(!gabriel.23vars.gen$serial_gen %in% gabriel.1707$serial_gen )
gabriel.1708.asthma.def.vars.ordered <-  gabriel.1708.asthma.def.vars[order(gabriel.1708.asthma.def.vars$serial_gen) ,] # reorder
rm(gabriel.1708.asthma.def.vars)
pos.out <- which(!gabriel.1708.asthma.def.vars.ordered$serial_gen %in% gabriel.1707$serial_gen )
asthma_wo_bronch <- gabriel.1708.asthma.def.vars.ordered$g1r12_01x[-pos.out]
asthma_wo_bronch[asthma_wo_bronch==2] <- 1
asthma_wo_bronch <- as.integer(asthma_wo_bronch)


######## genetic vs non genetic asthma

# genetic anyasthma
gen.any <- rep(NA, 1707)
gen.any[gabriel.1707$atopy07 == 1 & gabriel.1707$any_asthma==1] <- 1
gen.any[gabriel.1707$any_asthma ==0] <- 0

# nongenetic anyasthma
nongen.any <- rep(NA, 1707)
nongen.any[gabriel.1707$atopy07 == 0 & gabriel.1707$any_asthma==1] <- 1
nongen.any[gabriel.1707$any_asthma ==0] <- 0



# genetic anyasthma
gen.dda <- rep(NA, 1707)
gen.dda[gabriel.1707$atopy07 == 1 & gabriel.1707$dd_asthma==1] <- 1
gen.dda[gabriel.1707$dd_asthma ==0] <- 0

# nongenetic ddaasthma
nongen.dda <- rep(NA, 1707)
nongen.dda[gabriel.1707$atopy07 == 0 & gabriel.1707$dd_asthma==1] <- 1
nongen.dda[gabriel.1707$dd_asthma ==0] <- 0




# serial-number ok?
all(gabriel.1707$serial_gen == gabriel.1707.asthma$serial_gen)

dat.outcome <- data.frame(dd_asthma = gabriel.1707$dd_asthma, gabriel.1707.asthma[ ,-1], asthma.pheno, asthma_hay, asthma_eczema, 
                          asthma_hay_eczema, asthma_pure, asthma_nonatopic, dda_notany, ddahay_notdda, ddaecz_notdda,   dda.wh12_notany = dat.with.wh12$dda.wheeze_notany, asthma_wo_bronch, gen.any, nongen.any, gen.dda, nongen.dda)
save(dat.outcome, file="dat.outcome.rda")


## some checks:
attach(dat.outcome)
table(dd_asthma, dda.wh12_notany)[2,1] == 0
table(dd_asthma, asthma_wo_bronch)[1,2] ==0

