##################################### functions ####################################################

#####################################
######## generally useful functions
#####################################
library(glmnet)
# library(CMA)
library(ROCR)
# library(plsgenomics)
library(ggplot2)
library(randomForest)
library(pROC)


####################################################################################################
######## prediction
####################################################################################################
# 
# # for evaluation weights (observation multiplication):
# multiply.obs <- function(y, prediction, weight.evals){
#   weights2<- (unique(weight.evals))
#   
#   pos <- numeric()
#   for(i in 1:length(weights2)){
#     pos <- c(pos, rep(which(weight.evals==weights2[i]), round(weights2[i])))  
#   }
#   y <- y[pos]
#   prediction <- prediction[pos]
#   cbind(y=y,predicted= prediction)
# }
# 
# 
# ### unbiased empirical loss according to Zadrozny 2004 (principle)
# unb.loss <- function(w, lossfct, pred, y, ci=FALSE){
#   loss.tot <- 0
#   w.u <- unique(w)
#   for(j in 1:length(w.u)){
#     prob.strat <- length(which(w == w.u[j])) * w.u[j]/sum(w)
#     loss.tot <- loss.tot + prob.strat * lossfct(y[w == w.u[j]], pred[w == w.u[j]])
# #     if(ci==TRUE){
# #       
# #     }
# #     #  sum(w * (x - weighted.mean(x, w))^2)/sum(w)
#   } 
#   loss.tot
# }
# 
# ## roc-auc as loss function
# loss.roc <- function(y, pred, ci=FALSE){
#   roc1 <- roc(response = y, predictor = pred, direction="<", ci=ci)
#   if(ci==TRUE){
#     roc1 
#   } else roc1$auc
# } 
# 
# 
# ## simple partitioner of an index vector into CV-folds
# cv.partition <- function(ind, k){
#   fold.size <- floor(length(ind)/k)
#   folds <- list()
#   for(i in 1:(k-1)){
#     folds[[i]] <- ind[(1+(i-1)*fold.size):(fold.size*i) ]
#   }
#   folds[[k]] <- ind[(1 + (k-1)*fold.size): length(ind)]
#   folds
# }

############################ glmnet:


### my function for performing lasso
glmnet.my.cv <- function(y, x, nfolds, niters=1, weights=1, plot=TRUE, alpha=1, standardize=TRUE, measure.cv="auc",
                         weight.evals=1, just.auc=TRUE){
  
  learningsets <- list()
  for(i in 1:niters){
    learningsets[[i]] <- GenerateLearningsets(y=y, method="CV", strat=TRUE, fold=nfolds)
  }
  
  # tuningstep <- CMA::tune(X=x, y=as.factor(y), learningsets=learningssets, classifier=LassoCMA,
  #                         family=binomial, weights=weights, 
  #                         fold=3) 
  #### note: apparently it was not taken care of the weights argument in CMA
  
  if(length(weights)==1){
    weights <- rep(1, length(y))
  } 
  if(length(weight.evals)==1){
    weight.evals <- rep(1, length(y))
    weight.evals.use=FALSE
  } else {
    weight.evals.use=TRUE
  }
  
  ## lasso-glmnet
  prediction_glmnet <- matrix(NA, ncol=niters, nrow=length(y))
  
  classif <- list()
  
  for(i in 1:niters){
    for(j in 1:nfolds){
      idx_train <- learningsets[[i]]@learnmatrix[j, ]
      idx_test <- setdiff(1:length(y), idx_train)
      glmnet <- cv.glmnet(x=x[idx_train, ], y=(y[idx_train]), 
                          weights=weights[idx_train], alpha=alpha,
                          family='binomial', nfolds=3, standardize=standardize, type.measure=measure.cv)
      prediction_glmnet[idx_test, i] = predict(glmnet, x[idx_test, ])
      classif[[ nfolds*(i-1) + j]] <- cbind(y=y[idx_test], predicted = prediction_glmnet[idx_test, i], w.eval = weight.evals[idx_test])
      print(paste0("Iteration ", i, ", fold ", j, " done."))
    }
  }

  # should evaluations be weighted
  if(weight.evals.use==TRUE){
    weights2<- (unique(weight.evals))
    # for just.auc=TRUE:
    pos <- numeric()
    for(i in 1:length(weights2)){
      pos <- c(pos, rep(which(weight.evals==weights2[i]), round(weights2[i])))  
    }
    y <- y[pos]
    prediction_glmnet <- as.matrix(prediction_glmnet[ pos,])
    # otherwise:
    classif <- lapply(classif, function(x) multiply.obs(x[ , "y"], x[ , "predicted" ], x[ , "w.eval"]))
  }  else {
    classif <- lapply(classif, function(x) x[ ,1:2])
  }
  
  
  
  if(just.auc==TRUE){
    if(niters!=1) warning("iterations (niters) are ignored - set just.auc=FALSE to include them")
    pred <- ROCR::prediction(prediction_glmnet[ ,1], (y))
    perf <- performance(pred, x.measure = "fpr", measure = "tpr") 
    if(plot==TRUE) plot(perf)
    performance(pred,'auc')
  } else{
    classif
  }
}

######## function for unlisting fold-predictions
unlist.folds <- function(classif){
  cbind(
    unlist(lapply(classif, function(x) x[ ,1])),
    unlist(lapply(classif, function(x) x[ ,2]))
    )
  
}

######## function for AUC, CI and plot
my.roc <- function(prediction, cma=FALSE, plot=FALSE){
  if(cma==TRUE){ # if more than one fold: apply join() on cma-object 
    prediction.matrix <- cbind(prediction@y, prediction@prob[ ,2])
  } else{
    prediction.matrix <- prediction
  }
  plot1 <- plot
  pROC::roc(response=prediction.matrix[ ,1], predictor=prediction.matrix[ ,2], auc=TRUE , ci=TRUE, plot=plot1, direction="<")
}

# 

# ### my function for performing lasso
# glmnet.my.cv2 <- function(y, x, nfolds, niters=1, weights=1, plot=TRUE, alpha=1, standardize=TRUE, measure.cv="auc"){
#   
#   learningsets <- list()
#   for(i in 1:niters){
#     learningsets[[i]] <- GenerateLearningsets(y=y, method="CV", strat=TRUE, fold=nfolds)
#   }
#   
#   # tuningstep <- CMA::tune(X=x, y=as.factor(y), learningsets=learningssets, classifier=LassoCMA,
#   #                         family=binomial, weights=weights, 
#   #                         fold=3) 
#   #### note: apparently it was not taken care of the weights argument in CMA
#   
#   if(length(weights)==1){
#     weights <- rep(1, length(y))
#   } 
#   
#   ## lasso-glmnet
#   prediction_glmnet <- list()
#   for(i in 1:niters){
#     for(j in 1:nfolds){
#       
#       idx_train <- learningsets[[i]]@learnmatrix[j, ]
#       idx_test <- setdiff(1:length(y), idx_train)
#       glmnet <- cv.glmnet(x=x[idx_train, ], y=(y[idx_train]), 
#                           weights=weights[idx_train], alpha=alpha,
#                           family='binomial', nfolds=3, standardize=standardize, type.measure=measure.cv)
#       prediction_glmnet[idx_test, i] = predict(glmnet, x[idx_test, ])
#       print(paste0("Iteration ", i, ", fold ", j, " done."))
#     }
#   }
#   pred <- ROCR::prediction(prediction_glmnet[ ,1], (y))
#   perf <- performance(pred, x.measure = "fpr", measure = "tpr") 
#   if(plot==TRUE) plot(perf)
#   performance(pred,'auc')
# }
# 
# 


############################ PLS followed by LDA ######################
pls.lda.my.cv <- function(y, x, nfolds, niters=1, weights=1, plot=TRUE){
  # Generate Learning Sets
 learningsets <- GenerateLearningsets(y=y, method = "CV", fold = nfolds, 
                                      niter = niters, strat = TRUE)
  
  ############## Partial Least Squares with linear discriminant analysis #############################
  
  # perform it
  class.pls_lda <- classification(X=x, y=y, learningsets=learningsets,
                                  classifier=pls_ldaCMA)
  res <- join(class.pls_lda)

 pred <- ROCR::prediction(res@prob[ ,2], res@y)
 perf <- performance(pred, x.measure = "fpr", measure = "tpr")
 if(plot==TRUE) plot(perf)
 performance(pred,'auc')
}

############################ RandomForst ##############################

rf.my.cv <- function(y, x, nfolds, niters=1, weights=1, plot=TRUE){
  
  learningsets <- list()
  for(i in 1:niters){
    learningsets[[i]] <- GenerateLearningsets(y=y, method="CV", strat=TRUE, fold=nfolds)
  }
  
  # tuningstep <- CMA::tune(X=x, y=as.factor(y), learningsets=learningssets, classifier=LassoCMA,
  #                         family=binomial, weights=weights, 
  #                         fold=3) 
  #### note: apparently it was not taken care of the weights argument in CMA
  
  if(length(weights)==1){
    weights <- rep(1, length(y))
  } 
  
  ## rf
  prediction_rf <- matrix(NA, ncol=niters, nrow=length(y))
  for(i in 1:niters){
    for(j in 1:nfolds){
      idx_train <- learningsets[[i]]@learnmatrix[j, ]
      idx_test <- setdiff(1:length(y), idx_train)
    rf <- randomForest(x=x[idx_train, ], y=(y[idx_train]))
      prediction_rf[idx_test, i] = predict(rf, x[idx_test, ])
      print(paste0("Iteration ", i, ", fold ", j, " done."))
    }
  }
  pred <- ROCR::prediction(prediction_rf[ ,1], (y))
  perf <- performance(pred, x.measure = "fpr", measure = "tpr") 
  if(plot==TRUE) plot(perf)
  performance(pred,'auc')
}

############################ PLS followed by RFs ######################

############################ PCA followed by naiveBayes using first PC only:
# 
# for(i in 1:k){
#   idx_train = CVinds!=i
#   idx_test = CVinds==i
#   PCA_C1_tr = prcomp(countsC1_cc[idx_train,])
#   PCA_C1_ts = predict(PCA_C1_tr,countsC1_cc[idx_test,])
#   gnb <- naiveBayes(x=as.matrix(PCA_C1_tr$x[,1]),y=as.factor(lables[idx_train])) 
#   prediction[idx_test,] = predict(gnb, as.matrix(PCA_C1_ts[,1]), type='raw')
# }
# 
# pred <- prediction(prediction[ ,1], as.factor((lables==1)))
# perf <- performance(pred, x.measure = "fpr", measure = "tpr") # false postive rate, true positive rate
# performance(pred,'auc') # interpretation? -->AUC 
# plot(perf)
# 
# 
# pca.nbayes.my.cv <- function(y, x, nfolds, niters=1, weights=1, plot=TRUE, alpha=1){
#   
#   learningsets <- list()
#   for(i in 1:niters){
#     learningsets[[i]] <- GenerateLearningsets(y=y, method="CV", strat=TRUE, fold=nfolds)
#   }
#   
#   if(length(weights)==1){
#     weights <- rep(1, length(y))
#   } 
#   
#   ## lasso-glmnet
#   prediction_glmnet <- matrix(NA, ncol=niters, nrow=length(y))
#   for(i in 1:niters){
#     for(j in 1:nfolds){
#       idx_train <- learningsets[[i]]@learnmatrix[j, ]
#       idx_test <- setdiff(1:length(y), idx_train)
#       
#       PCA_tr = prcomp(countsC1_cc[idx_train,])
#       PCA_ts = predict(PCA_C1_tr,countsC1_cc[idx_test,])
#       gnb <- naiveBayes(x=as.matrix(PCA_tr$x[,1]),y=as.factor(y[idx_train])) 
#       
#       prediction[idx_test,] = predict(gnb, as.matrix(PCA_ts[,1]), type='raw')
#       
#       prediction_glmnet[idx_test, i] = predict(glmnet, x[idx_test, ])
#       print(paste0("Iteration ", i, ", fold ", j, " done."))
#     }
#   }
#   pred <- ROCR::prediction(prediction_glmnet[ ,1], (y))
#   perf <- performance(pred, x.measure = "fpr", measure = "tpr") 
#   if(plot==TRUE) plot(perf)
#   performance(pred,'auc')
# }
# 



####################################################################################################
######## visualization
####################################################################################################




########################## density plot for genetic risk
densityplot <- function(title= "Genetik risk for childhood asthma",
                        df, # data.frame containing:  categorial outcome in first column
                          #                           counts in second column
                        labels1 = c("no", "yes"), # vector of labels, should have same length as number categories
                        name.labels = FALSE
                        ){
      if(name.labels == FALSE){
        name.labels <- colnames(df)[1]
      } 
      df[ ,1] <- factor(df[ ,1], labels=labels1)
      names(df) <- c("y","x")
      ggplot(df, aes(x, fill = y ) ) +
      geom_density(alpha = 0.2) + ylab("Density") + xlab("Genetik risk of asthma in counts") +
      ggtitle(title) + 
      scale_fill_discrete(name=name.labels) 
}

## examples
# (df <- data.frame(asthma = c(0,0,0,0,1,1,1,1), counts =c(34,35,26,23,78,98,56,56)))
# densityplot(df=df)
# 
# df2 <- data.frame(asthma = c(0,0,0,0,1,1,1,1,2,2,2,2), counts =c(34,35,26,23,78,98,56,56,102, 204,302,203))
# densityplot(df=df2, labels1=c("no", "def1", "def2"), name.labels="asthma-types")

