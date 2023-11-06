# Credit: https://github.com/richardkwo/MultiSplit

#' Generate m-out-of-n subsamples
#' @param m size of subsample
#' @param n size of full sample
#' @param B number of subsamples 
#' @export
#' @return A (B x m) matrix where each row is a subsample.
#' @examples
#' get.m.out.n.tuples(10, 100, 500)
get.m.out.n.tuples <- function(m, n, B) {
  n.groups <- floor(n / m)
  n.perms <- ceiling(B / n.groups)
  tuple.mat <- replicate(n.perms, {
  idx <- sample(n, n.groups * m)
      matrix(idx, ncol=m)}, simplify = FALSE)
  tuple.mat <- do.call(rbind, tuple.mat)
  tuple.mat <- tuple.mat[1:B, ]
  return(tuple.mat)
}

# smart aggregation -----
get.smart.agg.pval <- function(obs.vec, ref.mat, reject.larger=TRUE) {
  # obs.vec: observed agg. stat of length W
  # ref.mat: subsampling agg. stat of size B x W
  stopifnot(length(obs.vec) == ncol(ref.mat))
  W <- length(obs.vec)
  if (!reject.larger) {
    obs.vec <- -obs.vec
    ref.mat <- -ref.mat
  }
  # columnwise calibration
  ecdf.funs <- lapply(1:W, function(w) stats::ecdf(ref.mat[,w]))
  ref.mat.calibrated <- sapply(1:W, function(w) ecdf.funs[[w]](ref.mat[,w]))
  R.ref <- apply(ref.mat.calibrated, 1, max)
  R.obs <- max(sapply(1:W, function(w) ecdf.funs[[w]](obs.vec[w])))
  return(mean(R.ref > R.obs))
}