
######## Wrapper functions for learning procedures specifically adjusted to snp and environment data
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/2-fcts-statlearning.R")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/4-snp-data-specific-functions.r")
source("~/../../../storage/cmbstore/norbert.krautenbacher/1-Projects/1-R-functions/5-prediction-and-validation.R")
require(parallel)

##### 
perform.asthma <- function(x.env=NULL, # environmental data: has to be a model matrix
                           x.conf, # data set of confounders or variables to adjust for
                           wheeze_onset=NULL, # 
                           snp.data, # snp data: has to be a model matrix
                           y, # outcome as data.frame or matrix with one column (in order to extract name)
                           env.if=c("snps_only", "env_only", "fh_only", "interactions_only", "both"), # env, snp, both or also wheeze onset age to include: 
                           blocks.ipf=NULL, # put colnumber of blocks in a list if env consists of more than one block (important for ipflasso)
                           snp.as=c("summary.", "estdos", "cat", "genomewide"), # how to include snps
                           method, # learning method
                           summary.type=NULL, # if snp.as=summary.: which kind of summary statistic to use
                           w.inv.var=NULL, #  if snp.as=summary.:should inverse variance of the snp be used as weight?
                           gen.mod=NULL, #  if snp.as=summary.: genetic model
                           ind.test=NULL, # if NULL: 5-Fold-CV, else: list (!!) of specification of the test sample; rest always train
                           k=5, # k of k-Fold CV; ignored if is.null(ind.test) ==FALSE
                           obs.weights=1, # observation weights 
                           cor.thr=NULL, # for glm_gwas: correlation threshold
                           interactions=FALSE,
                           needed.as.matrix=TRUE,
                           needed.as.df=TRUE,
                           conv.x.all.to.model.matrix=TRUE,
                           ntrees=500, 
                           n.cores=1, # number of cores used for CV/resampling procedure
                           snp.pvals=NULL, # snp.list containing pvals from gwas
                           index.cv.list.comparing=NULL, # optional but recommended for glm_gwas: compare indices used there with the ones used here (should be equal)
                           index.cv.list.comparing2=NULL, # see below (and here it is about the test index)
                           pflist.long=FALSE, # default pflist for ipflasso 
                           ncv.ipflasso=0,
                           to.be.validated=FALSE,
                           path.results="Results/" # path where to save the results
){ 
  # outcome name
  y.which <- colnames(y)
  y <- y[ ,1]
  

  
  # data 
  aucs.a <- aucs.b.i <- aucs.b.ii <- data.frame(method = character(), 
                                                snps = numeric(), snps.lo = numeric(), snps.up = numeric(),
                                                env = numeric(), env.lo = numeric(), env.up = numeric(),
                                                fh_only = numeric(), fh_only.lo = numeric(), fh_only.up = numeric(),
                                                interactions_only = numeric(), interactions_only.lo = numeric(), interactions_only.up = numeric(),
                                                both = numeric(), both.lo = numeric(), both.up = numeric()
                                                , stringsAsFactors=FALSE
  )
  # create file for saving aucs (if not already existent)
  file.path <- paste0(path.results, "aucs/aucs_" , y.which, "_", snp.as,  
                      summary.type, w.inv.var, gen.mod, ".rda")
  ifelse(!dir.exists(file.path(path.results, "aucs")), dir.create(file.path(path.results, "aucs")), FALSE)
  if(!file.exists(file.path)) save(aucs.a,aucs.b.i, aucs.b.ii, file=file.path)
#   
  
  ################# CV-set-up
  
  # remove NAs from outcome
  ind.y <- which(!is.na(y)) 
  y <- y[ind.y]
  


  
  ########## create x.snp
  if(snp.as == "estdos"){
    x.snp <-  model.matrix(~.+0, data.frame(snp.data)[,-1])
  } else if(snp.as == "genomewide"){
    x.snp <- snp.data
  } else if(snp.as == "cat"){
    x.snp.temp <-round(snp.data[,-1])
    for(i in 1:ncol(x.snp.temp)) x.snp.temp[ ,i] <- as.factor(x.snp.temp[ ,i])
    x.snp <- model.matrix(~.+0,x.snp.temp)
  } else if(class(summary.type)!="NULL"){                   # hier if summary.type=NULL?........
     x.snp <- as.matrix(summary.rs(data=snp.data[,-1], type=summary.type, w.inv.var=w.inv.var, gen.mod=gen.mod))
  }
  
  
  #### remove observations where outcome is NA ## should be done before
  # obs weights
  if(length(obs.weights)==1){
    obs.weights <- rep(1, length(y))
  } else{
    obs.weights <- obs.weights[ind.y]
  }
  
  
  # create x.all
  if(env.if == "both"){
    if(!(all(ind.y==1:1507))){
      x.all <- cbind(x.conf, x.env, x.snp)[ind.y, ]
      x.snp <- cbind(x.conf, x.snp)[ind.y, ] # since we need it for some pca-approaches
      x.env <- cbind(x.conf, x.env)[ind.y, ]
    } else{
      x.all <- cbind(x.conf, x.env, x.snp)
      x.snp <- cbind(x.conf, x.snp)
      x.env <- cbind(x.conf, x.env)
      
    }
  } 
  
  if(env.if=="conf_only"){
    x.all <- x.conf[ind.y, ]
    rm(x.snp)
  }
  
  if(env.if == "env_only" | env.if == "fh_only" |env.if ==  "interactions_only"){
    x.all <- cbind(x.conf, x.env)[ind.y, ]
    rm(x.snp)
  } 
  
  if(env.if == "snps_only"){
    if(!(all(ind.y==1:1507))){
      x.all <- cbind(x.conf, x.snp)[ind.y, ]
      
    } else{
      x.all <- cbind(x.conf, x.snp)
    }
    rm(x.env)
  } 
 
######### CV (or general arbitrary partitioning)  
  if(is.null(ind.test)){
    # 5-Fold-CV:
    k <- k
    n <- length(y)
    set.seed(3534)
    ind <- sample(1:n)
    
    ind.test <- list()
    ind.train <- list()
    for(i in 1:k){
      if(i < k ) ind.test[[i]] <- ind[(1+(i-1)*floor(n/k)):(i*(floor(n/k))) ]
      else ind.test[[i]] <- ind[(1+(i-1)*floor(n/k)):n]
      ind.train[[i]] <- setdiff(ind, ind.test[[i]])
    }
  } else{
    n <- length(y)
    ind <- 1:n
    k <- length(ind.test)
    
    ind.train <- list()
    for(i in 1:k){
      ind.train[[i]] <- setdiff(ind, ind.test[[i]])
    } 
  }
  
  # check indices
  if(!is.null(index.cv.list.comparing)){
    anyvec <- logical()
    for(i in 1:5) anyvec[i] <- all(index.cv.list.comparing$ind.within$ind.train[[i]] == ind.train[[i]])
    if(all(anyvec==TRUE)) print("Indices all fine") else stop("Indices of CV by this function unequal to the ones of fin.ind object")
  }
  # check indices more flexible
  # check indices
  if(!is.null(index.cv.list.comparing2)){
    anyvec <- logical()
    for(i in 1:5) anyvec[i] <- all(index.cv.list.comparing2[[i]] == ind.test[[i]])
    if(all(anyvec==TRUE)) print("Indices all fine") else stop("Indices of CV by this function unequal to the ones of fin.ind object")
  }
  
  ################# set-up for algorithms
  form <- formula(y ~.)
  if(needed.as.df==TRUE){
    dat <- data.frame(y=factor(y), x.all)
  }
  
  
  
  ######################################### perform ... ###########################################
# pred <- list() 
  # for(i in 1:k){
    perform.models <- function(i){
    # if x argument has to be set, make it formula specific
    if(needed.as.matrix==TRUE){
      if(conv.x.all.to.model.matrix==TRUE){
        x.pred <- as.matrix(model.matrix(form, data.frame(x.all))[ , -1]) # rem. intercept
      }else{
       x.pred = x.all 
      }
    }
    
    
    # 1. 
    ####  a) univariate snps selection by glm - final glm
    ####  b) univariate snps selection by svyglm - final svyglm
    
    if(method=="glm"){
      train.glm <- glm(formula=form,data=dat[ ind.train[[i]], ], family="binomial")
      pred <- predict(train.glm, newdata = dat[ ind.test[[i]], ], type="response")
    }
    
    if(method=="glm_weights"){
      
      form2 <- formula(y~.-w)
      train.glm <- glm(formula=form2,data=data.frame( dat, w=obs.weights)[ ind.train[[i]], ], family="binomial", weights=w)
      pred <- predict(train.glm, newdata = data.frame(dat, w=1)[ ind.test[[i]], ], type="response")
      # pred
      
      # form2 <- formula(y~.- w)
      # train.glm <- glm(formula=form2,data=data.frame(dat, w=obs.weights)[ ind.train[[i]], ], family="binomial", weights = w)
      # pred <- predict(train.glm, newdata = dat[ ind.test[[i]], ], type="response")
    }
    #   
    # 2.  a) LASSO    
      
      
      # test: lasso not standardized predictors
      if(method=="glmnet_weights_no_stand"){
        require(glmnet)
        train.glmnet <- cv.glmnet(y=y[ind.train[[i]]], x=x.pred[ ind.train[[i]], ], family="binomial", 
                                  type.measure="auc", alpha=1, weights=obs.weights[ ind.train[[i]]],
                                  standardize=FALSE)
        pred <- predict(train.glmnet, newx = x.pred[ ind.test[[i]], ], type="response")
      }
      
      
      
      
      
    if(method=="glmnet"){
      require(glmnet)
      train.glmnet <- cv.glmnet(y=y[ind.train[[i]]], x=x.pred[ ind.train[[i]], ], family="binomial", type.measure="auc", alpha=1)
      pred <- predict(train.glmnet, newx = x.pred[ ind.test[[i]], ], type="response")
    }
    #     b) LASSO-weights
    if(method=="glmnet_weights"){
      require(glmnet)
      train.glmnet <- cv.glmnet(y=y[ind.train[[i]]], x=x.pred[ ind.train[[i]], ], family="binomial", type.measure="auc", alpha=1, weights=obs.weights[ ind.train[[i]]])
      pred <- predict(train.glmnet, newx = x.pred[ ind.test[[i]], ], type="response")
    }
    
    ## ipf-lasso
    if(method=="ridge_weights"){
      require(glmnet)
      train.glmnet <- cv.glmnet(y=y[ind.train[[i]]], x=x.pred[ ind.train[[i]], ], family="binomial", type.measure="auc", alpha=0, weights=obs.weights[ ind.train[[i]]])
      pred <- predict(train.glmnet, newx = x.pred[ ind.test[[i]], ], type="response")
    }
    
    
    ### correct:
    # - confounders all right?? 
    # - weights --> use cvr2.ipflasso.ed
    
    
    if(method=="ipflasso" & env.if=="both"){
      require(ipflasso)
        if( !is.null(blocks.ipf) ){
          blocks = blocks.ipf
          pflist = list(c(1,1,1,1,1),  # first block: confounders and  weight always set to 1 
                        c(1,2,1,1,1), c(1,1,2,1,1), c(1,1,1,2,1), c(1,1,1,1,2), # weight one block more
                        c(1,2,2,1,1), c(1,2,1,2,1), c(1,2,1,1,2), c(1,1,2,2,1), c(1,1,2,1,2), c(1,1,1,2,2), # weight two blocks more
                        c(1,3,1,1,1), c(1,1,3,1,1), c(1,1,1,3,1), c(1,1,1,1,3), # weight one block three times
                        c(1,3,3,1,1), c(1,3,1,3,1), c(1,3,1,1,3), c(1,1,3,3,1), c(1,1,3,1,3), c(1,1,1,3,3) # weight two blocks three times
          )
        } else{
         stop("the env.blocks argument has to be not null if ipf-lasso should be performed")
        }
      
      train.glmnet <- cvr2.ipflasso(Y=y[ind.train[[i]]], X=x.pred[ ind.train[[i]], ], family="binomial", type.measure="auc", alpha=1, 
                                    blocks=blocks, pflist=pflist,
                                    nfolds=5,ncv=0)
      pred <- ipflasso.predict(train.glmnet, Xtest = x.pred[ ind.test[[i]], ])$probabilitiestest
    } 
    if(method=="ipflasso_weights" & env.if=="both"){
      require(ipflasso)
      if( !is.null(blocks.ipf) ){
        blocks = blocks.ipf
        if(pflist.long==FALSE){
          pflist = list(c(1,1,1,1,1),  # first block: confounders and  weight always set to 1 
                        c(1,2,1,1,1), c(1,1,2,1,1), c(1,1,1,2,1), c(1,1,1,1,2), # weight one block more
                        c(1,2,2,1,1), c(1,2,1,2,1), c(1,2,1,1,2), c(1,1,2,2,1), c(1,1,2,1,2), c(1,1,1,2,2), # weight two blocks more
                        c(1,3,1,1,1), c(1,1,3,1,1), c(1,1,1,3,1), c(1,1,1,1,3), # weight one block three times
                        c(1,3,3,1,1), c(1,3,1,3,1), c(1,3,1,1,3), c(1,1,3,3,1), c(1,1,3,1,3), c(1,1,1,3,3) # weight two blocks three times
          )
        } else if(pflist.long==TRUE){
          pflist = list(c(1,1,1,1,1),  # first block: confounders and  weight always set to 1 
                        c(1,2,1,1,1), c(1,1,2,1,1), c(1,1,1,2,1), c(1,1,1,1,2), # weight one block more
                        c(1,2,2,1,1), c(1,2,1,2,1), c(1,2,1,1,2), c(1,1,2,2,1), c(1,1,2,1,2), c(1,1,1,2,2), # weight two blocks more
                        c(1,3,1,1,1), c(1,1,3,1,1), c(1,1,1,3,1), c(1,1,1,1,3), # weight one block three times
                        c(1,3,3,1,1), c(1,3,1,3,1), c(1,3,1,1,3), c(1,1,3,3,1), c(1,1,3,1,3), c(1,1,1,3,3), # weight two blocks three times
                        c(1,4,1,1,1), c(1,1,4,1,1), c(1,1,1,4,1), c(1,1,1,1,4), # weight one block four times
                        c(1,4,4,1,1), c(1,4,1,4,1), c(1,4,1,1,4), c(1,1,4,4,1), c(1,1,4,1,4), c(1,1,1,4,4), # weight two blocks four times
                        c(1,5,1,1,1), c(1,1,5,1,1), c(1,1,1,5,1), c(1,1,1,1,5), # weight one block five times; another one twice
                        c(1,5,5,1,1), c(1,5,1,5,1), c(1,5,1,1,5), c(1,1,5,5,1), c(1,1,5,1,5), c(1,1,1,5,5), # weight two blocks five times
                        c(1,6,1,1,1), c(1,1,6,1,1), c(1,1,1,6,1), c(1,1,1,1,6), # weight one block five times; another one twice
                        c(1,6,6,1,1), c(1,6,1,6,1), c(1,6,1,1,6), c(1,1,6,6,1), c(1,1,6,1,6), c(1,1,1,6,6), # weight two blocks five times
                        c(1,7,1,1,1), c(1,1,7,1,1), c(1,1,1,7,1), c(1,1,1,1,7), # weight one block five times; another one twice
                        c(1,7,7,1,1), c(1,7,1,7,1), c(1,7,1,1,7), c(1,1,7,7,1), c(1,1,7,1,7), c(1,1,1,7,7) # weight two blocks five times
                                  )
        }
      } else{
        stop("the env.blocks argument has to be not null if ipf-lasso should be performed")
      }
   
      train.glmnet <- cvr2.ipflasso.ed(Y=y[ind.train[[i]]], X=x.pred[ ind.train[[i]], ], family="binomial", type.measure="auc", alpha=1, 
                                    blocks=blocks, pflist=pflist,
                                    nfolds=5,ncv=ncv.ipflasso, weights=obs.weights[ ind.train[[i]]])
      pred <- ipflasso.predict(train.glmnet, Xtest = x.pred[ ind.test[[i]], ])$probabilitiestest
    }
    
 
    
    # 2.  a) elnet    
    if(method=="elnet"){
      require(glmnet)
      train.glmnet <- cv.glmnet(y=y[ind.train[[i]]], x=x.pred[ ind.train[[i]], ], family="binomial", type.measure="auc", alpha=0.5)
      pred <- predict(train.glmnet, newx = x.pred[ ind.test[[i]], ], type="response")
    }
    #     b) LASSO-weights
    if(method=="elnet_weights"){
      require(glmnet)
      train.glmnet <- cv.glmnet(y=y[ind.train[[i]]], x=x.pred[ ind.train[[i]], ], family="binomial", type.measure="auc", alpha=0.5, weights=obs.weights[ ind.train[[i]]])
      pred <- predict(train.glmnet, newx = x.pred[ ind.test[[i]], ], type="response")
    }
    
    
    # 3. RF
    if(method=="ranger"){
      require(ranger)
      # train.ranger <- ranger(formula=form,data=dat[ ind.train[[i]], ], write.forest = TRUE, probability = TRUE)
  # for big data:
      train.ranger <- ranger(dependent.variable.name = "y", data = dat[ ind.train[[i]], ], write.forest = T, 
                             save.memory = TRUE, probability = TRUE, num.trees =ntrees, num.threads = 1, verbose = FALSE) 
      pred <- predict(train.ranger, data= dat [ind.test[[i]], ])$predictions[ ,2]
    }
    if(method=="parIPranger"){
      require(sambia)
      pred <- sambia::synthIPbag(data = data.frame(dat[ ind.train[[i]], ], obs.weights = obs.weights[ ind.train[[i]]]), # weights here required as stratum info
                                 weights = obs.weights[ ind.train[[i]]], type='parIP', 
                                 strata.variables = "obs.weights", # weights give info about stratum here
                                    learner='ranger',
                                 list.train.learner = list(formula=formula(as.factor(y)~. -obs.weights), num.threads = num.threads),
                                    list.predict.learner = list(data=dat [ind.test[[i]], ]),n.bs = ntrees)
    }
    
    # 4. a) PCA - followed by logit
    if(method=="pca_logit"){
      prc <- prcomp(x.all)
      prc.loads <- prc$x[ ,1:(which(summary(prc)$importance[3, ] > .8 )[1])]
      train.glm <- glm(formula=form, family="binomial", data=data.frame(y=factor(y), prc.loads)[ ind.train[[i]], ])
      pred <- predict(train.glm, newdata=data.frame(prc.loads)[ ind.test[[i]], ])
    }
    #    b) PCA-snps - followed by logit
    if(method=="pca_snps_logit" & (env.if=="both" | env.if=="all.wo") ){
      prc <- prcomp(x.snp)
      prc.loads <- prc$x[ ,1:(which(summary(prc)$importance[3, ] > .8 )[1])]
      
      if(env.if=="both"){
        data.pca <- data.frame(y=factor(y), x.env, prc.loads)
        new.data.pca <- data.frame(x.env, prc.loads)
      } else if(env.if=="all.wo"){
        data.pca <- data.frame(y=factor(y), x.env,x.wo, prc.loads)
        new.data.pca <- data.frame(x.env, x.wo,prc.loads)
      }
            
      
      train.glm <- glm(formula=form, family="binomial", data=data.pca[ ind.train[[i]], ])
      pred <- predict(train.glm, newdata=new.data.pca[ ind.test[[i]], ])
    }
    #    c) PCA - followed by RF
    if(method=="pca_rf"){
      prc <- prcomp(x.all)
      prc.loads <- prc$x[ ,1:(which(summary(prc)$importance[3, ] > .8 )[1])]
      require(ranger)
      train.ranger <- ranger(formula=form,data=data.frame(y=factor(y), prc.loads)[ ind.train[[i]], ],
                             write.forest = TRUE, probability = TRUE , num.threads = 8)
      pred <- predict(train.ranger, data= data.frame(prc.loads)[ ind.test[[i]], ])$predictions[ ,2]
    }
  
    # d) PCA-snps - followed by RF
    if(method=="pca_snps_rf" &  (env.if=="both" | env.if=="all.wo")){
      prc <- prcomp(x.snp)
      prc.loads <- prc$x[ ,1:(which(summary(prc)$importance[3, ] > .8 )[1])]
      
      if(env.if=="both"){
        data.pca <- data.frame(y=factor(y), x.env, prc.loads)
        new.data.pca <- data.frame(x.env, prc.loads)
       } 
      # else if(env.if=="all.wo"){
#         data.pca <- data.frame(y=factor(y), x.env,x.wo, prc.loads)
#         new.data.pca <- data.frame(x.env, x.wo,prc.loads)
#       }
      
      require(ranger)
      train.ranger <- ranger(formula=form,data= data.pca[ ind.train[[i]], ], write.forest = TRUE, probability = TRUE, num.threads = 8)
      pred <- predict(train.ranger, data= new.data.pca[ ind.test[[i]], ])$predictions[ ,2]
    }
      
      
   if(method=="glm_gwas"){ # careful: Wu 2013 uses j for subjects and j for snp
      top.snp.no <- 100
      
      top.snp.names <- names(sort(snp.pvals[[i]], decreasing = FALSE)[1:top.snp.no])
      
      if(!is.null(cor.thr)){
        #  to do: either prune wider set than top.snp.no via correlation or ld....
        # not used: pfeiffer paper uses 95% threshold (and data was pruned for this threshold from the beginning)
      }
      
      
      snp.coefs <- rep(NA, top.snp.no)
      for(j in 1:top.snp.no){
        ## only ajdusted for snps_only yet (--> change x.conf to general ...)
        dat.train.for.snp.j <- data.frame(y, x.pred[, c(top.snp.names[j], colnames(x.conf))], stringsAsFactors = FALSE)[ ind.train[[i]], ]
        # print( "line 326 worked")
        glmj <- glm(y ~ ., data=dat.train.for.snp.j, family = "binomial", weights = obs.weights[ind.train[[i]]])
        snp.coefs[j] <- coef(glmj)[2]
        names(snp.coefs)[j] <- names(coef(glmj)[2])
      }
      x.snp.top <- x.snp[ , names(snp.coefs)]
      genetic.score.train <-  apply(x.snp.top[ ind.train[[i]], ], 1,function(x) sum(snp.coefs*x))
      
      dat.for.gwas.pred.train <- data.frame(y=y[ind.train[[i]]], genetic.score = genetic.score.train, x.pred[ ,  which(!(colnames(x.pred) %in% colnames(x.snp))) ][ind.train[[i]], ])
      # print("line 335 worked")
      glm.fin <- glm(y ~ . - genetic.score + offset(I(1*genetic.score)), data = dat.for.gwas.pred.train, family="binomial", weights = obs.weights[ind.train[[i]]] ) 
      ## note: offset option used in order to keep the genetic score coefficient fixed
      # print("glm fit worked")
      genetic.score.test <-  apply(x.snp.top[ ind.test[[i]], ], 1,function(x) sum(snp.coefs*x))
      
#       if(env.if=="snps_only"){
#         pred <- genetic.score.test
#       } else if(env.if=="both"){
        dat.for.gwas.pred.test <- data.frame(genetic.score = genetic.score.test, x.pred[ ,  which(!(colnames(x.pred) %in% colnames(x.snp) )) ][ind.test[[i]], ])
        pred <- predict(glm.fin, newdata = dat.for.gwas.pred.test, type="response"  )
      # }
         print("all in gwas approach worked")
    } 
      
      if(method=="glm_mv_gwas"){ # multivariate glm after GWAS
        top.snp.no <- 100
        
        top.snp.names <- names(sort(snp.pvals[[i]], decreasing = FALSE)[1:top.snp.no])
        ## only ajdusted for snps_only yet (--> change x.conf to general ...)
        dat.for.top.snps <- data.frame(y, x.pred[, c(top.snp.names, colnames(x.conf))], stringsAsFactors = FALSE)
        glm.top.snps <- glm(y ~ ., data=dat.for.top.snps[ ind.train[[i]], ], family = "binomial", weights = obs.weights[ind.train[[i]]])
        
        pred <- predict(glm.top.snps, newdata =  dat.for.top.snps[ ind.test[[i]], ], type = "response")      
      } 
    
    if(method=="plsr"){
      require(pls)
      dat1 <- dat
      dat1$y <- as.numeric(dat$y) - 1
      train.plsr <- plsr(formula=form, data=dat1[ ind.train[[i]], ], validation="CV", ncomp=5)
      pred <- c(predict(train.plsr, newdata = dat[ ind.test[[i]], ], comps=5))
    }


    if(method=="ready.score"){
      pred <- x.pred[ ind.test[[i]], 1 ] # (CV has no effect here, however)
    }
    
    if(method=="ready.score.learn.risk.al"){
      x.snp <- round(x.snp)
      snp.list.train <- recode.risk.allele(x.snp[ ind.train[[i]], ], y=y[ ind.train[[i]] ], output="snps.to.recode")
      x.to.test <- x.snp # in order not do decode x.snp for the new iterations
      for(snp in snp.list.train){ # recode
        x.to.test[ ind.test[[i]], snp][ x.to.test[ ind.test[[i]], snp] == 2] <- 3
        x.to.test[ ind.test[[i]], snp][ x.to.test[ ind.test[[i]], snp] == 0] <- 2
        x.to.test[ ind.test[[i]], snp][ x.to.test[ ind.test[[i]], snp] == 3] <- 0
      }
      pred <- summary.rs(data=x.to.test[ ind.test[[i]], ], type=summary.type, w.inv.var=w.inv.var, gen.mod=gen.mod)
    }
    
    # 5.  a) pvclust - logit
    #     b) pvclust - RF
    print(paste("Iteration", i, "done!"))
    pred
  }
  
  pred <- mclapply(1:k, perform.models, mc.cores = n.cores)  
  # 
  if(length(pred) == 0){} else{
    ifelse(!dir.exists(file.path(path.results, "predictions")), dir.create(file.path(path.results, "predictions")), FALSE)
    save(pred, file=paste0(path.results, "predictions/pred.", method,"_", y.which,"_", snp.as , "_", env.if , 
                           summary.type, w.inv.var, gen.mod, ".rda"))
    
    if(to.be.validated==TRUE){
      # true outcome
      response <- y[unlist(ind.test)] 
      # response
      
      
      ### auc-file
      
      
      
      
      ###### evaluation: 
      #### a) AUC (biased but length is ok)
      
      load(file.path)
      
      print(roc.a <- loss.roc(y=response, pred=unlist(pred), ci=TRUE))
      #### b) weighted AUC
      # b.I) ipb (length is ok and unbiased)
      
      ##### note!! obs.weights is fine to be of size n since after CV we have exactly n observations to eval.
      print(roc.b.i <- ci.auc.bootstrap.w(roc=roc.a, boot.n=10000, parallel = FALSE, conf.level=.95, weights=obs.weights[unlist(ind.test)]))
      # b.II) multiplied (length is too small, but unbiased)
      mult.preds <- multiply.obs(y=response, prediction=unlist(pred), weight.evals = obs.weights[unlist(ind.test)])
      print (roc.b.ii <- loss.roc(y= mult.preds[ , "y"], pred=mult.preds[ , "predicted"], ci=TRUE))
      
      if(!(method %in% aucs.a$method)){
        rowno <- nrow(aucs.a) + 1
        aucs.a[rowno, ] <- aucs.b.i[rowno, ] <-aucs.b.ii[rowno, ] <- NA
        aucs.a$method[rowno] <- aucs.b.i$method[rowno] <- aucs.b.ii$method[rowno] <- method
      } else rowno <- which(aucs.a$method==method)
      if(env.if=="conf_only"){
        aucs.a$conf[rowno] <- roc.a$auc
        aucs.a$conf.lo[rowno] <- roc.a$ci[1]
        aucs.a$conf.up[rowno] <- roc.a$ci[3]
        
        aucs.b.i$conf[rowno] <- roc.b.i[2]
        aucs.b.i$conf.lo[rowno] <- roc.b.i[1]
        aucs.b.i$conf.up[rowno] <- roc.b.i[3]
        #   
        aucs.b.ii$conf[rowno] <- roc.b.ii$auc
        aucs.b.ii$conf.lo[rowno] <- roc.b.ii$ci[1]
        aucs.b.ii$conf.up[rowno] <- roc.b.ii$ci[3]
      } else if(env.if=="snps_only"){
        aucs.a$snps[rowno] <- roc.a$auc
        aucs.a$snps.lo[rowno] <- roc.a$ci[1]
        aucs.a$snps.up[rowno] <- roc.a$ci[3]
        
        aucs.b.i$snps[rowno] <- roc.b.i[2]
        aucs.b.i$snps.lo[rowno] <- roc.b.i[1]
        aucs.b.i$snps.up[rowno] <- roc.b.i[3]
        #   
        aucs.b.ii$snps[rowno] <- roc.b.ii$auc
        aucs.b.ii$snps.lo[rowno] <- roc.b.ii$ci[1]
        aucs.b.ii$snps.up[rowno] <- roc.b.ii$ci[3]
      } else if(env.if=="env_only"){
        aucs.a$env[rowno] <- roc.a$auc
        aucs.a$env.lo[rowno] <- roc.a$ci[1]
        aucs.a$env.up[rowno] <- roc.a$ci[3]
        
        aucs.b.i$env[rowno] <- roc.b.i[2]
        aucs.b.i$env.lo[rowno] <- roc.b.i[1]
        aucs.b.i$env.up[rowno] <- roc.b.i[3]
        
        aucs.b.ii$env[rowno] <- roc.b.ii$auc
        aucs.b.ii$env.lo[rowno] <- roc.b.ii$ci[1]
        aucs.b.ii$env.up[rowno] <- roc.b.ii$ci[3]
        
      } else   if(env.if=="fh_only"){
        aucs.a$fh_only[rowno] <- roc.a$auc
        aucs.a$fh_only.lo[rowno] <- roc.a$ci[1]
        aucs.a$fh_only.up[rowno] <- roc.a$ci[3]
        
        aucs.b.i$fh_only[rowno] <- roc.b.i[2]
        aucs.b.i$fh_only.lo[rowno] <- roc.b.i[1]
        aucs.b.i$fh_only.up[rowno] <- roc.b.i[3]
        
        aucs.b.ii$fh_only[rowno] <- roc.b.ii$auc
        aucs.b.ii$fh_only.lo[rowno] <- roc.b.ii$ci[1]
        aucs.b.ii$fh_only.up[rowno] <- roc.b.ii$ci[3]
      } else   if(env.if=="interactions_only"){
        aucs.a$interactions_only[rowno] <- roc.a$auc
        aucs.a$interactions_only.lo[rowno] <- roc.a$ci[1]
        aucs.a$interactions_only.up[rowno] <- roc.a$ci[3]
        
        aucs.b.i$interactions_only[rowno] <- roc.b.i[2]
        aucs.b.i$interactions_only.lo[rowno] <- roc.b.i[1]
        aucs.b.i$interactions_only.up[rowno] <- roc.b.i[3]
        
        aucs.b.ii$interactions_only[rowno] <- roc.b.ii$auc
        aucs.b.ii$interactions_only.lo[rowno] <- roc.b.ii$ci[1]
        aucs.b.ii$interactions_only.up[rowno] <- roc.b.ii$ci[3]
      }else   if(env.if=="both"){
        aucs.a$both[rowno] <- roc.a$auc
        aucs.a$both.lo[rowno] <- roc.a$ci[1]
        aucs.a$both.up[rowno] <- roc.a$ci[3]
        
        aucs.b.i$both[rowno] <- roc.b.i[2]
        aucs.b.i$both.lo[rowno] <- roc.b.i[1]
        aucs.b.i$both.up[rowno] <- roc.b.i[3]
        
        aucs.b.ii$both[rowno] <- roc.b.ii$auc
        aucs.b.ii$both.lo[rowno] <- roc.b.ii$ci[1]
        aucs.b.ii$both.up[rowno] <- roc.b.ii$ci[3]
      }   
      save(aucs.a, aucs.b.i, aucs.b.ii, file=file.path)
    }
  }
}

