## -----------------------------------------
## Credit: https://github.com/richardkwo/MultiSplit, https://github.com/anna-neufeld/countsplit_paper
## -----------------------------------------

source("multisplit.R")

datathin <- function(X, Lambda, ep, gammas=rep(1,NROW(X)), K=2) {

    if (K == 2) {
        # apply the k-fold data thinning to the subsample
        Xtrain <- apply(X, 2, function(u) rbinom(n=length(u), size=u, p=ep))
        Xtest <- X - Xtrain

        # apply the latent variable estimation and differential expression analysis to train and test set
        # Note this is for continuous trajectories only
        hXtrain <- log(diag(1/gammas) %*% (Xtrain+1))
        hXtraincenter <- apply(hXtrain, 2, function(u) u-mean(u))
        pseudotime <- svd(hXtraincenter)$u[,1]
        coeff.est <- apply(
            Xtest, 2, function(u) summary(glm(u~pseudotime, offset=log(gammas), family="poisson"))$coefficients[2,3]
        )
        popu.para.est <- suppressWarnings(apply(Lambda, 2, function(u) summary(glm(u~pseudotime, family="poisson"))$coefficients[2,1]))

        trueprincomp <- svd(apply(log(Lambda), 2, function(u) u-mean(u)))$u[,1]
        true.cor <- cor(pseudotime, trueprincomp)

        abs(cbind(coeff.est, popu.para.est, true.cor))
    } else {
        # apply the k-fold data thinning to the subsample
        X.folds <- apply(X, 2, function(u) rmultinom(n=length(u), size=u, prob=rep(1/K, K)))
    }

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
#' @param reject.twoside Use 'FALSE' if the test is single-sided. The default value is
#' TRUE.
#' @param m The size of a subsample. The default value is floor(n / ln(n)) where n is 
#' the sample size.
#' @param L Number of repetitions of applying the test statistic function test.single  # nolint
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
    reject.twoside=TRUE,
    m=NULL,
    J=5,
    L=50,
    verbose=FALSE
) {
    n <- nrow(data) # sample size
    p <- ncol(data)
    if (is.null(m)) { # set subsample size
        m <- floor(n / log(n))
    } else {
        m <- as.integer(m)
    }
    
    B <- J * floor(n / m) # number of subsamples
    if (K != 2) { L = K } # set the number of test repetitions

    tuple.mat <- get.m.out.n.tuples(m, n, B) # generate the indices of elements in each sample
    
    if (verbose) {
        message(sprintf("Subsampling with %d tuples (m = %d) and %d splits ...", B, m, L))
    }

    pval.jobs <- plyr::llply(1:B, function(b) {
        future::future({
            data.sub <- data[tuple.mat[b,], ]  # of size m
            Lambda.sub <- Lambda[tuple.mat[b,], ]   
            gammas.sub <- gammas[tuple.mat[b,]] 

            if (k == 2) {
                replicate(L, {datathin(data.sub, Lambda.sub, eps, gammas.sub, K)})  
            } else {
                
            }
               
        }, seed=TRUE)
    }, .progress=ifelse(verbose, "text", "none"))

    T.mat <- plyr::laply(pval.jobs, future::value)

    # observed T
    T.obs.vec <- replicate(L, datathin(data, Lambda, eps, gammas)[, 1])
    p.values <- numeric(p)

    # p.values <- apply(
    #     T.mat[1:B, 1:p, 1, 1:L],
    #     2,
    #     function(array) {
    #         T.mat.transformed <- (rank(array, ties.method = "random") - 1/2) / length(array)  # rank transform
    #         T.mat.transformed <- matrix(T.mat.transformed, nrow=B)

    #          # apply inverse CDF
    #         T.mat.transformed <- stats::qnorm(T.mat.transformed)
    #         if (verbose) {
    #             message(sprintf("Max change from rank transform = %g", max(abs(T.mat.transformed - array))))
    #         }
    #         # aggregation 
    #         if (!is.list(S)) {
    #             # single aggregation function ------
    #             if (verbose) {
    #                 message(sprintf("Single aggregation function (reject for %s values)", ifelse(reject.larger, "larger", "smaller")))
    #             }
    #             T.obs <- S(T.obs.vec[j,])
    #             T.sub <- apply(T.mat.transformed, 1, S)
    #             if (reject.larger) {
    #                 if (reject.twoside) {
    #                     p.value <- mean(abs(T.sub) > abs(T.obs))
    #                 } else {
    #                     p.value <- mean(T.sub > T.obs)
    #                 }
    #             } else {
    #                 if (reject.twoside) {
    #                     p.value <- mean(abs(T.sub) < abs(T.obs))
    #                 } else {
    #                     p.value <- mean(T.sub < T.obs)
    #                 }
    #             }
    #         } else {
    #             # multiple aggregation functions -----
    #             if (verbose) {
    #             message(sprintf("Adapting to the best among %d aggregation functions (reject for %s values)", length(S), ifelse(reject.larger, "larger", "smaller")))
    #             }
    #             T.obs <- sapply(S, function(.s) .s(T.obs.vec))
    #             T.sub <- sapply(S, function(.s) apply(T.mat.transformed, 1, .s))
    #             p.value <- get.smart.agg.pval(T.obs, T.sub, reject.larger=reject.larger)
    #         }
    #         p.value
    #     }
    # )

    abs.coeff.ests <- T.mat[1:B, 1:p, 1, 1:L]
    for (j in seq_len(p)) {
        # retrieve value
        T.mat.j <- abs.coeff.ests[1:B, j, 1:L]

        T.mat.j.transformed <- (rank(T.mat.j, ties.method = "random") - 1/2) / length(T.mat.j)  # rank transform
        T.mat.j.transformed <- matrix(T.mat.j.transformed, nrow=B)

        # apply inverse CDF
        T.mat.j.transformed <- fdrtool::qhalfnorm(T.mat.j.transformed)
        if (verbose) {
            message(sprintf("Max change from rank transform = %g", max(abs(T.mat.j.transformed - T.mat.j))))
        }

        # aggregation 
        if (!is.list(S)) {
            # single aggregation function ------
            if (verbose) {
                message(sprintf("Single aggregation function (reject for %s values)", ifelse(reject.larger, "larger", "smaller")))
            }
            T.obs <- S(T.obs.vec[j,])
            T.sub <- apply(T.mat.j.transformed, 1, S)
            if (reject.larger) {
                if (reject.twoside) {
                    p.value <- mean(abs(T.sub) > abs(T.obs))
                } else {
                    p.value <- mean(T.sub > T.obs)
                }
            } else {
                if (reject.twoside) {
                    p.value <- mean(abs(T.sub) < abs(T.obs))
                } else {
                    p.value <- mean(T.sub < T.obs)
                }
            }
        } else {
            # multiple aggregation functions -----
            if (verbose) {
            message(sprintf("Adapting to the best among %d aggregation functions (reject for %s values)", length(S), ifelse(reject.larger, "larger", "smaller")))
            }
            T.obs <- sapply(S, function(.s) .s(T.obs.vec))
            T.sub <- sapply(S, function(.s) apply(T.mat.j.transformed, 1, .s))
            p.value <- get.smart.agg.pval(T.obs, T.sub, reject.larger=reject.larger)
        }

        p.values[j] <- p.value

    }

    abs.popu.para.ests.mean <- apply(T.mat[1:B, 1:p, 2, 1:L], 2, mean)
    abs.true.corrs.mean <- apply(T.mat[1:B, 1:p, 3, 1:L], 2, mean)

    cbind(p.values, abs.popu.para.ests.mean, abs.true.corrs.mean) # placeholder for correlation
}