library(tidyverse)
source("algo.R")

#' Generate the biological variation in a scRNA-seq data model.
gen_Lambda <- function(n=200,p=100, k=1, Fs, intercepts) {
if (k==1) {
  L <- rnorm(n)
  L <- L - mean(L)
}

LFT <- L%*%t(Fs)
logLambda <- t(apply(LFT, 1, function(u) u+intercepts))
Lambda <- exp(logLambda)

return(Lambda)
return(list(dat=X,mean=EX, coeffs=Fs, intercepts=intercepts))
}

#' Generate the observed data in a scRNA-seq data model.
gen_pois_data <- function(gammas, Lambda) {
  EX <- apply(Lambda, 2, function(u) u*gammas)
  X <- apply(EX,2,function(u) rpois(length(u),u))
  return(X)
}

one_trial <- function(
    n,
    p,
    filename, 
    foldername,
    k=1,
    K=2,
    propImp=0.1, 
    sig_strength=5, 
    propLowMedHigh = c(1/2,1/2), 
    eps=c(0.5),
    L=50,
    verbose=FALSE) {
  intercepts <- sample(c(log(3), log(25)), prob=propLowMedHigh, replace=TRUE, size=p)
  num_non_null <- floor(propImp * p)
  numNull <- p - num_non_null
  
  Fs <- c(rep(sig_strength, num_non_null), rep(0, numNull))
  
  Lambda <- gen_Lambda(n, p, k, Fs, intercepts)
  gammas <- rgamma(n, 10, scale=1/10)
  
  X <- gen_pois_data(gammas, Lambda)
  
  if (file.exists(filename)) {
    result <- read.csv(filename)
  } else {
    result <- NULL
  }
  
  for (ep in eps) {
    res <- cbind(
        1:p, 
        Fs, 
        intercepts, 
        datathin.multisplit(
            X,
            gammas,
            Lambda,
            K=K,   
            eps=ep,
            L=L,
            verbose=verbose,
            foldername=foldername      
        ), 
        ep, 
        "known", 
        propImp, 
        n,
        p, 
        propLowMedHigh[1]
    )

    res <- res %>%
      as_tibble() %>%
      setNames(c("j", "trueCoeff", "intercept", "pval",
        "estPopuPara", "cor", "eps", "type", "propImp",
        "n", "p", "prop1"))

    result <- rbind(result, res)
  }

  write.csv(result, filename, row.names=FALSE)
}

