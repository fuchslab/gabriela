library(pROC)
library(ROCR)
library(parallel)

####################################################################################################
######## prediction and validation
####################################################################################################


# bootstrapped CIs
# loss.boot <- function(y, prediction, weights){
#   require(sampling)
# }





############# for evaluation weights (observation replication):
multiply.obs <- function(y, prediction, weight.evals){
  weights2<- (unique(weight.evals))
  if(length(weights2)!=1){ # only do procedure only if there is more than one weight
    pos <- numeric()
    for(i in 1:length(weights2)){
      pos <- c(pos, rep(which(weight.evals==weights2[i]), round(weights2[i])))  
    }
    y <- y[pos]
    prediction <- prediction[pos]
  } 
  cbind(y=y,predicted= prediction)
}

multiply.obs2 <- function(roc, weights, get.auc=TRUE){ # roc-object (ROCR); observation weights
  predictor <- roc$predictor
  response <- roc$response
  y.repl <- rep(response, weights)
  pred.repl <- rep(predictor, weights)
  if(get.auc==TRUE) roc(response=y.repl, predictor= pred.repl, direction="<")$auc
  else cbind(y=y.repl, predicted = pred.repl)
} 


############## unbiased empirical loss according to Zadrozny 2004 (principle)
unb.loss <- function(w, lossfct, pred, y, ci=FALSE){
  loss.tot <- 0
  w.u <- unique(w)
  for(j in 1:length(w.u)){
    prob.strat <- length(which(w == w.u[j])) * w.u[j]/sum(w)
    loss.tot <- loss.tot + prob.strat * lossfct(y[w == w.u[j]], pred[w == w.u[j]])
    #     if(ci==TRUE){
    #       
    #     }
    #     #  sum(w * (x - weighted.mean(x, w))^2)/sum(w)
  } 
  loss.tot
}

################# roc-auc as loss function
loss.roc <- function(y, pred, ci=FALSE){
  roc1 <- roc(response = y, predictor = pred, direction="<", ci=ci)
  if(ci==TRUE){
    roc1 
  } else roc1$auc
} 


##################################################################################################################################
################################ weighted bootstrap ROC-AUC with CIs by modified functions from pROC-package #####################
##################################################################################################################################

#### implementation of weighted bootstrap-CI (also point estimate for AUC calculated by bootstrap (as the 50%-quantile))


# example
# set.seed(235)
# y <- c(rep(0,50), rep(1,50))
# yhat <- c(rnorm(50), rnorm(50, mean=1) )
# a<- roc(response = y,predictor = yhat, direction="<")
# a
# 
# my.weights <- c(100, rep(1,99))

# 
# require(pROC)
# methods(ci.auc)
# getAnywhere(ci.auc.auc)
# getAnywhere(ci.auc.roc)
# 
#   getAnywhere(ci.auc.delong)
#   getAnywhere(ci.auc.bootstrap)
#    getAnywhere(auc.roc)
# 
#     getAnywhere(nonstratified.ci.auc) # there it is: the sample function which we have to ## modify ###    
nonstratified.ci.auc.w <- function (n, roc, weights, n.cases.min) # n: not appearing (probably=boot.n) # roc: roc-object (just to get predictor and response back)
{
  tmp.idx <- sample(1:length(roc$predictor), replace = TRUE, prob=weights) # here we have to set the weights; NOTE: not to mix up with using 1/w:
  ######## probabilities are needed but in terms of calculating the weights relatively! (1/weights would be inverse-inverse-probability weights -> wrong)
 
  # if high skewness should be prevented make a minimum number of cases obligatory: 
  if(!(is.null(n.cases.min)) ){
    while( min(table(roc$response[tmp.idx])) / sum(table(roc$response[tmp.idx])) < n.cases.min){
      tmp.idx <- sample(1:length(roc$predictor), replace = TRUE, prob=weights)
    } 
 } 
  predictor <- roc$predictor[tmp.idx]
  response <- roc$response[tmp.idx]
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(roc$levels[1])]]
  cases <- splitted[[as.character(roc$levels[2])]]
  thresholds <- pROC:::roc.utils.thresholds(c(controls, cases))
  perfs <- roc$fun.sesp(thresholds = thresholds, controls = controls, 
                        cases = cases, direction = roc$direction)
  roc$sensitivities <- perfs$se
  roc$specificities <- perfs$sp
  pROC:::auc.roc(roc, partial.auc = attr(roc$auc, "partial.auc"), 
                 partial.auc.focus = attr(roc$auc, "partial.auc.focus"), 
                 partial.auc.correct = attr(roc$auc, "partial.auc.correct"), 
                 allow.invalid.partial.auc.correct = TRUE)
}



# nonstratified.ci.auc.w(n=2000, roc=a, weights=weights) # works.
# getAnywhere(roc.utils.thresholds)

###### --> one level higher --> modify this function correspondingly --> get the CI   
# getAnywhere(ci.auc.bootstrap)
ci.auc.bootstrap.w <- function (roc, conf.level, boot.n, # boot.stratified, 
                                progress= "none", 
                                parallel, weights, return.aucs=FALSE, n.cases.min=NULL, ...) 
{
  if (class(progress) != "list") 
    progress <- pROC:::roc.utils.get.progress.bar(progress, title = "AUC confidence interval", 
                                                  label = "Bootstrap in progress...", ...)
  #   if (boot.stratified) {
  #     aucs <- unlist(plyr:::llply(1:boot.n, .fun = pROC:::stratified.ci.auc, 
  #                          roc = roc, .progress = progress, .parallel = parallel))
  #   }
  #   else {
  aucs <- unlist(plyr:::llply(1:boot.n, .fun = nonstratified.ci.auc.w, 
                              roc = roc, .progress = progress, .parallel = parallel, weights=weights, n.cases.min=n.cases.min ))
  #   }
  if (sum(is.na(aucs)) > 0) {
    warning("NA value(s) produced during bootstrap were ignored.")
    aucs <- aucs[!is.na(aucs)]
  }
  if(return.aucs==TRUE){
    return(aucs)
  } else{
    return(quantile(aucs, c(0 + (1 - conf.level)/2, 0.5, 1 - 
                              (1 - conf.level)/2)))
  }
 
}




##### version 2 ("by hand")
my.bs.2samples.roc <- function(yhat1,yhat2=NULL,y,B, prob=NULL, alpha.2samples=.05, cores=1){
  if(is.null(yhat2)){
    ind <- lapply(1:B, function(z) sample(1:length(y), size=length(y), replace=T, prob=prob))
    aucs <- mclapply(ind, function(z){
      roc(y[z], yhat1[z], direction = "<")$auc
    }, mc.cores = cores)
    quantile(unlist(lapply(aucs, function(x) x)), probs=c(.025, .5, .975))
  } else{
    ind <- lapply(1:B, function(z) sample(1:length(y), size=length(y), replace=T, prob=prob))
    aucs <- mclapply(ind, function(z){
      auc1 <- roc(y[z], yhat1[z], direction = "<")$auc
      auc2 <- roc(y[z], yhat2[z], direction = "<")$auc 
      # print(auc1, auc2)
      c(auc1, auc2)
    }, mc.cores = cores )
    diffs <- unlist(lapply(aucs, function(x) x[2] - x[1] ))
    results <- list()
    
    results$aucs1 <- quantile(unlist(lapply(aucs, function(x) x[1])), probs=c(.025, .5, .975))
    results$aucs2 <- quantile(unlist(lapply(aucs, function(x) x[2])), probs=c(.025, .5, .975))
    results$diffaucs <- quantile(diffs, probs = c(alpha.2samples/2, 0.5, 1-alpha.2samples/2))
    results 
  }
}



#### Differently implmented function giving back CI and all boostrapped AUCs (e.g. useful for Bayesfactor calculation)
boot_AUC = function(event, marker, nboot = 200, siglev = 0.95, weights, cores=1){
  tt1 = res = list()
  auci = numeric(nboot)
  iisamp = matrix(NA, ncol = nboot, nrow = length(event))
  if(cores!=1){
    stop("Careful: Don't use 'parallel' here - gets fucked up when trying to have the same assignments by same seed!!!")
    auc.boot.fct <- function(b){
      ii = sample(1:length(event),size = length(event),replace = T, prob = weights)
      #auci[i] = pROC::roc(event[ii], marker[ii], direction = "<")$auc
      rr = ROCR::prediction(marker[ii], event[ii])
      unlist(ROCR::performance(rr, measure = "auc")@y.values)
    }
    res$AUC = unlist(mclapply(1:nboot, auc.boot.fct, mc.cores=cores ))
  } else if(cores==1){
    for (i in 1:nboot){
      ii = sample(1:length(event),size = length(event),replace = T)
      #auci[i] = pROC::roc(event[ii], marker[ii], direction = "<")$auc
      rr = ROCR::prediction(marker[ii], event[ii])
      auci[i] = unlist(ROCR::performance(rr, measure = "auc")@y.values)
    }
    res$AUC = auci
  }
  res$ciAUC = c(sort(res$AUC)[floor(nboot*((1-siglev)/2))],sort(res$AUC)[floor(nboot*(1-(1-siglev)/2))])
  return(res)
}
### application:
# set.seed(123)
# ROboot = boot_AUC(teddy$mult_second_AB_6y, teddy$rsRO, nboot=nb)
# set.seed(123)
# WIboot = boot_AUC(teddy$mult_second_AB_6y, teddy$rs40, nboot=nb)
# 
# # Bayes factors
# (sum(WIboot$AUC>ROboot$AUC)/nb)/(sum(WIboot$AUC<=ROboot$AUC)/nb)






### test
# ci.auc.bootstrap.w(roc=a,  # boot.stratified = FALSE, 
#                    boot.n=2000, parallel = FALSE, conf.level=.95, weights=my.weights)



# roc(y,yhat, boot.stratified = FALSE, method="bootstrap", , ci.auc=TRUE)

## siehe auch Trello zu erklärungen








################## simple partitioner of an index vector into CV-folds
cv.partition <- function(ind, k){
  fold.size <- floor(length(ind)/k)
  folds <- list()
  for(i in 1:(k-1)){
    folds[[i]] <- ind[(1+(i-1)*fold.size):(fold.size*i) ]
  }
  folds[[k]] <- ind[(1 + (k-1)*fold.size): length(ind)]
  folds
}

################# Precision-recall AUC (function from dream challenge)
score_q2<-function(pred,y)
{
  ## remove subjects whose discontinuation status is NA
  pred <- pred[!is.na(y)]
  y <- y[!is.na(y)]
  
  prf=ROCR::performance(ROCR::prediction(pred,y),"prec","rec")
  x=prf@x.values[[1]]
  y=prf@y.values[[1]]
  auc=sum((x[-1]-x[-length(x)])*y[-1])
  auc
}
