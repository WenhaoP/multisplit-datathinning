## -----------------------------------------
## Credit: https://github.com/richardkwo/MultiSplit, https://github.com/anna-neufeld/countsplit_paper
## -----------------------------------------

source("multisplit.R")

countsplit <- function(X, Lambda, ep, gammas=rep(1,NROW(X))) {
    # apply the 2-fold data thinning to the subsample
    
    Xtrain <- apply(X,2,function(u) rbinom(n=length(u), size=u, p=ep))
    Xtest <- X-Xtrain

    # apply the latent variable estimation and differential expression analysis to train and test set
    # Note this is for continuous trajectories only
    hXtrain <- log(diag(1/gammas)%*%(Xtrain+1))
    hXtraincenter <- apply(hXtrain,2,function(u) u-mean(u))
    pseudotime <- svd(hXtraincenter)$u[,1]

    est.coeff <- apply(Xtest, 2, function(u) summary(glm(u~pseudotime, offset=log(gammas), family="poisson"))$coefficients[2,1])
    true_coeffs <- suppressWarnings(apply(Lambda, 2, function(u) summary(glm(u~pseudotime, family="poisson"))$coefficients[2,1]))

    cbind(pvals_pseudotime,true_coeffs)
}


#' Combination of data thinning and rank-transformed subsampling algorithms
#' 
#' @param data A data frame or matrix with iid rows.
#' #' @param gammas Size factors in a scRNA-seq data model.
#' @param Lambda Biological variation in a scRNA-seq data model.
#' @param K Number of folds for the data thinning procedure. The default value is 2.
#' @param eps Data thinning tuning parameter. The default value is 1/K.
#' @param S Aggregation function(s). It can be either a single aggregation function
#'    (e.g., `S=mean`) or a list of aggregation functions (e.g., `S=c(mean, max)`).
#'    For the latter case, a p-value that automatically adapts to the most 
#'    powerful aggregation function in the list will be returned. The default value is mean.
#' @param reject.larger The side of the test. Use `FALSE` if `test.single` returns
#'    a p-value. The default value is TRUE.
#' @param m The size of a subsample. The default value is floor(n / ln(n)) where n is 
#' the sample size.
#' @param L Number of repetitions of applying the test statistic function test.single 
#' to a subsample. If K = 2, the default value is 50. If K > 2, set L to be K.
#' @param J Number of collections of subsamples.
#' @param verbose If `TRUE`, will print details. 
datathin.multisplit <- function(
    data,
    gammas,
    Lambda,
    K=2,
    eps=1/K,
    S=function(x) mean(x, na.rm=TRUE),
    reject.larger=TRUE,
    m=NULL,
    J=5,
    L=50,
    verbose=FALSE
) {
    n <- nrow(data) # sample size
    if (is.null(m)) { # set subsample size
        m <- floor(n / log(n))
    } else {
        m <- as.integer(m)
    }
    
    B <- J * floor(n / m) # number of subsamples
    if (K != 2) { L = K } # set the number of test repetitions

    popu_coeffs <- c()

    
    tuple.mat <- get.m.out.n.tuples(m, n, B) # generate the indices of elements in each sample
    
    if (verbose) {
        message(sprintf("Subsampling with %d tuples (m = %d) and %d splits ...", B, m, n.splits))
    }

    pval.jobs <- plyr::llply(1:B, function(b) {
        future::future({
            data.sub <- data[tuple.mat[b,], ]  # of size m
            print(data.sub)
            print(data)
            if (K == 2){
                # iterate over splits
                replicate(L, {

                    # apply the 2-fold data thinning to the subsample
                    result <- countsplit(data.sub, Lambda, eps, gammas=rep(1,NROW(data.sub)))

                    # store the population coefficient
                    popu_coeffs = c(popu_coeffs, result[2])

                    result[1]
                })
            } else {
                # apply the K-fold data thinning to the subsample
                # TODO

                # iterate over splits
                replicate(L, {
                    # aggregate thinned sets into train and test set 
                    # TODO

                    # apply the latent variable estimation and differential expression analysis to train and test set
                    # TODO

                    # store the population coefficient
                    # TODO
                })
            }
        }, seed=TRUE)

    }, .progress = ifelse(verbose, "text", "none"))

    # retrieve value
    T.mat <- plyr::laply(pval.jobs, future::value)
    T.mat.transformed <- (rank(T.mat, ties.method = "random") - 1/2) / length(T.mat)  # rank transform
    T.mat.transformed <- matrix(T.mat.transformed, nrow=B)

    # apply inverse CDF
    T.mat.transformed <- stats:qnorm(T.mat.transformed)
    if (verbose) {
        message(sprintf("Max change from rank transform = %g", max(abs(T.mat.transformed - T.mat))))
    }

    # observed T
    T.obs.vec <- replicate(L, countsplit(data, Lambda, eps, c=1, gammas=rep(1,NROW(data)))[1])

    # aggregation 
    if (!is.list(S)) {
        # single aggregation function ------
        if (verbose) {
        message(sprintf("Single aggregation function (reject for %s values)", ifelse(reject.larger, "larger", "smaller")))
        }
        T.obs <- S(T.obs.vec)
        T.sub <- apply(T.mat.transformed, 1, S)
        if (reject.larger) {
        p.value <- mean(T.sub > T.obs)
        } else {
        p.value <- mean(T.sub < T.obs)
        }
    } else {
        # multiple aggregation functions -----
        if (verbose) {
        message(sprintf("Adapting to the best among %d aggregation functions (reject for %s values)", length(S), ifelse(reject.larger, "larger", "smaller")))
        }
        T.obs <- sapply(S, function(.s) .s(T.obs.vec))
        T.sub <- sapply(S, function(.s) apply(T.mat.transformed, 1, .s))
        p.value <- get.smart.agg.pval(T.obs, T.sub, reject.larger=reject.larger)
    }

    popu_coeffs.mean = mean(popu_coeffs) # Might need to change the dimension the average is w.r.t. for higher-dimensioned beta_1

    cbind(p.value, popu_coeffs.mean, NULL) # placeholder for correlation
}