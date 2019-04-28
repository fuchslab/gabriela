##### ipf-lasso with weights


library(ipflasso)
##### extending functions there

# fix(cvr.ipflasso)
cvr.ipflasso.ed <- function (X, Y, family, type.measure, standardize = TRUE, alpha = 1, 
                             blocks, pf, nfolds, ncv, ...) 
{
  M <- length(blocks)
  if (M != length(pf)) {
    stop("blocks and pf must have the same length.")
  }
  ulblocks <- as.numeric(unlist(blocks))
  if (!setequal(ulblocks, 1:ncol(X))) {
    stop("Each predictor should be included in exactly one block.")
  }
  if (family == "gaussian") {
    if (type.measure != "mse") 
      warning("type.measure is set to mse.")
    type.measure <- "mse"
  }
  if (family == "cox") {
    if (type.measure != "deviance") 
      warning("type.measure is set to partial likelihood.")
  }
  if (family == "binomial" & !is.element(type.measure, c("auc", 
                                                         "class"))) {
    warning("type.measure is set to class")
    type.measure <- "class"
  }
  if (any(pf <= 0)) {
    stop("pf should have positive entries.")
  }
  pfvector <- numeric(ncol(X))
  for (m in 1:M) {
    pfvector[blocks[[m]]] <- pf[m]
  }
  a <- cvr.glmnet(X = X, Y = Y, family = family, standardize = standardize, 
                  alpha = alpha, nfolds = nfolds, ncv = ncv, type.measure = type.measure, 
                  penalty.factor = pfvector,...)
  coeff <- a$coeff
  if (family != "cox") {
    rownames(coeff)[1] <- "intercept"
  }
  if (type.measure == "auc") {
    ind.bestlambda <- which.max(a$cvm)
  }
  if (type.measure != "auc") {
    ind.bestlambda <- which.min(a$cvm)
  }
  nzero <- apply(coeff[-1, ], FUN = function(x) return(sum(x != 
                                                             0)), MARGIN = 2)
  return(list(coeff = coeff, ind.bestlambda = ind.bestlambda, 
              lambda = a$lambda, cvm = a$cvm, nzero = nzero, family = family))
}

# fix(cvr2.ipflasso)
cvr2.ipflasso.ed <- function (X, Y, family, type.measure, standardize = TRUE, alpha = 1, 
                              blocks, pflist, nfolds, ncv, nzeromax = +Inf, plot = FALSE,...) 
{
  M <- length(blocks)
  nw <- length(pflist)
  if (!setequal(M, sapply(pflist, FUN = length))) {
    stop("The length of the entries of argument pflist must equal the number of blocks.")
  }
  ulblocks <- as.numeric(unlist(blocks))
  if (!setequal(ulblocks, 1:ncol(X))) {
    stop("Each predictor should be included in exactly one block.")
  }
  if (family == "gaussian") {
    if (type.measure != "mse") 
      warning("type.measure is set to mse.")
    type.measure <- "mse"
  }
  if (family == "cox") {
    if (type.measure != "deviance") 
      warning("type.measure is set to partial likelihood.")
  }
  if (family == "binomial" & !is.element(type.measure, c("auc", 
                                                         "class"))) {
    warning("type.measure is set to class")
    type.measure <- "class"
  }
  a <- list()
  cvmin <- +Inf
  for (j in 1:nw) {
    a[[j]] <- cvr.ipflasso.ed(Y = Y, X = X, family = family, 
                               type.measure = type.measure, standardize = standardize, 
                               alpha = alpha, blocks = blocks, pf = pflist[[j]], 
                               nfolds = nfolds, ncv = ncv, ...)
    allowedindices <- which(as.numeric(a[[j]]$nzero) <= nzeromax)
    if (type.measure == "auc") {
      ajcvm <- -a[[j]]$cvm
    }
    if (type.measure != "auc") {
      ajcvm <- a[[j]]$cvm
    }
    mincvmj <- min(ajcvm[allowedindices])
    if (mincvmj < cvmin) {
      ind.bestpf <- j
      ind.bestlambda <- which.min(ajcvm[allowedindices])[1]
      bestlambda <- a[[j]]$lambda[ind.bestlambda]
      cvmin <- mincvmj
    }
  }
  if (plot == TRUE) {
    par(mfrow = c(1, 2))
    plot(a[[1]]$nzero, a[[1]]$cvm, type = "l", ylim = c(0, 
                                                        1), ylab = type.measure, xlab = "total number of included variables")
    if (nw > 1) {
      for (j in 2:nw) {
        points(a[[j]]$nzero, a[[j]]$cvm, type = "l", 
               col = j)
      }
      abline(v = nzeromax, lty = 2)
    }
    pfnames <- c()
    for (j in 1:nw) {
      pfnames <- c(pfnames, paste(pflist[[j]], collapse = "-"))
    }
    legend(col = 1:nw, lty = 1, legend = pfnames, y = 0.95, 
           x = max(a[[1]]$nzero)/3)
    plot(a[[ind.bestpf]]$nzero, apply(a[[ind.bestpf]]$coeff[blocks[[1]] + 
                                                              1, ], MARGIN = 2, FUN = function(x) sum(x != 0)), 
         ylim = c(0, max(a[[ind.bestpf]]$nzero)), col = ind.bestpf, 
         pch = 2, xlab = "total number of included variables", 
         ylab = "number of variables")
    for (b in 2:length(blocks)) {
      points(a[[ind.bestpf]]$nzero, apply(a[[ind.bestpf]]$coeff[blocks[[b]] + 
                                                                  1, ], MARGIN = 2, FUN = function(x) sum(x != 
                                                                                                            0)), pch = b + 1, col = ind.bestpf)
    }
    abline(a = 0, b = 1)
    legend(pch = 1 + (1:length(blocks)), legend = paste("block", 
                                                        1:length(blocks)), x = 0, y = max(a[[ind.bestpf]]$nzero), 
           col = ind.bestpf)
  }
  coeff <- a[[ind.bestpf]]$coeff
  return(list(coeff = coeff, ind.bestlambda = ind.bestlambda, 
              bestlambda = bestlambda, ind.bestpf = ind.bestpf, a = a, 
              family = family))
}




# 
# 
# ## example for testing:
# n <- 100
# y <- rnorm(n)
# weights <- sample(1:3, replace=T, size=n)
# 
# x <- cbind(rnorm(100), rnorm(100), rexp(100), rnorm(100))
# 
# 
# set.seed(3412)
# cvr.ipf <- cvr.glmnet(X=x, Y=y, family = "gaussian",  type.measure = "mse", nfolds = 4, ncv=0) # still works in or. version
# cvr.ipf$coeff[ , "s5"]
# 
# set.seed(3412)
# cvr.ipf <- cvr.glmnet(X=x, Y=y, family = "gaussian",  type.measure = "mse", nfolds = 4, ncv=0,
#            weights=weights) # still works in or. version
# cvr.ipf$coeff[ , "s5"] # different, as it should
# 
# # 
# 
# set.seed(3434)
# cvr.ipf <- cvr.ipflasso.ed(X=x, Y=y, family = "gaussian", blocks=list(1:2,3:4) , type.measure = "mse", nfolds = 4, ncv=0, pf = c(1,1))
# cvr.ipf$coeff[ , "s5"]
# set.seed(3434)
# cvr.ipf <- cvr.ipflasso.ed(X=x, Y=y, family = "gaussian", blocks=list(1:2,3:4) , type.measure = "mse", nfolds = 4, ncv=0, pf = c(1,1), weights = weights)
# cvr.ipf$coeff[ , "s5"]
# 
# 
# set.seed(2343)
# cv.ipf <- cvr2.ipflasso.ed(X=x, Y=y, family = "gaussian", blocks=list(1:2,3:4), pflist = list(c(1,1), c(2,1) ), type.measure = "mse", nfolds = 4, ncv=0)
# cv.ipf$a[[2]]$coeff[ , "s6"]
# set.seed(2343)
# cv.ipf <- cvr2.ipflasso.ed(X=x, Y=y, family = "gaussian", blocks=list(1:2,3:4), pflist = list(c(1,1), c(2,1) ), type.measure = "mse", nfolds = 4, ncv=0, weights=weights)
# cv.ipf$a[[2]]$coeff[ , "s6"]

