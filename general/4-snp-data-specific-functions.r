###### simple summary statistic risk scores/polygenic risk score (optional: weighted by inverse variance of each snp)

# polygenic risk score (additive model)
# polygenic risk score (additive model; using counts instead of est. dosages )
# polygenic risk score (recessive model; using counts instead of est. dosages)
# polygenic risk score (dominant model; using counts instead of est. dosages)
# variance
# sum over u- or t-statistics
  summary.rs <- function(data, type, w.inv.var=FALSE, gen.mod="add.estdos"){
    # data: SNPs as data frame
    # gen.mod = c("add.estdos", "add.counts", "rec", "dom"): genetic model (additive one can be with counts or estimated dosages (default))
    # type = c("pol.rs","var", "ustat", "tstat" ): way of building score
    # w.inv.var: with/without inverse variance weighting (variance of each SNP)
    # counts: use counts instead of estimated dosages
    
    if(gen.mod=="add.counts"){
      data <- round(data)
    } else if(gen.mod=="rec"){
      data <- round(data)
      data[data==1] <- 0
      data[data==2] <- 1
    } else if(gen.mod=="dom"){
      data <- round(data)
      data[data==2] <- 1
    }
    
    if(w.inv.var==TRUE){
      inv.var <- 1/apply( data, 2, var)
      inv.var[inv.var==Inf] <- 0 # i.e.: if var = 0, SNP is non-informative so we weight it by factor 0.
    } else inv.var <- 1
   
    ## apply methods of the score
    if(type=="pol.rs"){
      rs <- apply( data, 1, function(x) sum(inv.var*x)) 
    }else if(type=="var"){
      rs <- apply( data, 1, function(x) var(inv.var*x)) 
    }else if(type=="ustat"){
      rs <- apply( data, 1, function(x) var(inv.var*x)) 
    }else if(type=="tstat"){
      rs <- apply( data, 1, function(x) t.test(inv.var*x)$statistic) 
    }
    rs
  }


##### find risk allele for each snp
recode.risk.allele <-  function(x,y,output="new.snp.matrix"){
  # 1. round all snps to counts (reason: adding up a snp of e.g. estimated dosage 1.55 instead of 2 just because knowledge is vague, would bias such a risk score)
  x<-round(x)
  snps.to.recode <- character()
  for(snp in colnames(x)){
    # 2. fit univariate logistic regression models
    beta1 <- coefficients(glm(y ~ x[ ,snp ], family="binomial"))[2]
    # 3. check if beta_1 is greater or less than 0 and recode if lower 
    if(beta1 < 0){
      x[ , snp][ x[ , snp] == 2] <- 3
      x[ , snp][ x[ , snp] == 0] <- 2
      x[ , snp][ x[ , snp] == 3] <- 0
      snps.to.recode <- c(snps.to.recode, snp)
    }
  }
  if(output=="new.snp.matrix"){
    x
  } else if(output=="snps.to.recode"){
    snps.to.recode
  } 
}

##### remove correlated snps
rem.cor <- function(x, keep.snps=NULL, t=.9, p.limit=500){
  # 1. generate nx3 data.frame: all unique combination of snps + column for correlations
  # 2. extract snp pair from data set x and calculate correlation coefficient (stratified?)
  # 3. reduce matrix by removing rows with correlation less than threshold t
  # 4. compare list of pairs with selection of snps to be kept (keep.snps)
  # - if one snp of the pair in keep.snps --> remove the other one
  # - if none of the snps or both of the snps are in keep.snps --> remove randomly
  if(p.limit > ncol(x)){
    design <- t(combn(ncol(x) ,2))
    snp1 <-colnames(x)[design[, 1]]
    snp2 <-colnames(x)[design[, 2]]
    cors <- apply(cbind(snp1,snp2), 1, function(z) cor(x[ , z])[1,2] ) 
    snp1.cor <- snp1[cors >t]#
    snp2.cor <- snp2[cors > t]
    if(!is.null(keep.snps)){
      #if(...)
    } else{
      x.no.cors <- x[ , setdiff(colnames(x), snp1.cor)]
    }
  }else{
    stop("Number of variables too high; set p.limit higher if procedure should be performed anyway!")
  }
  x.no.cors
}
  

