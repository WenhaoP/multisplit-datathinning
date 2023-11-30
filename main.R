## -----------------------------------------
## Original credit: https://github.com/anna-neufeld/countsplit_paper
## -----------------------------------------


## -----------------------------------------
## load user-defined functions, packages
## -----------------------------------------

## tidy stuff
library("tidyr")
library("fdrtool")
## Needed functions.
source("run_sims.R")

## -----------------------------------------
## set up a grid of parameters to cycle over
## -----------------------------------------

propLowMedHighs <- 1:3
regCoeffs <- c(
    log(1.2), log(1.35), log(1.5), log(1.6), log(1.7), log(1.85), 
    log(2), log(2.5), log(3), log(4), log(5), log(7), log(10), 
    log(15), log(20))
propImps <- c(0.1)
ns <- c(200)
ps <- c(10)

## number of monte-carlo iterations per job
nreps_per_combo <- 10

## simulation name
simname <- "test"

## epsilon candidates
# eps=c(0.1,0.25,0.5,0.75,0.9)
eps=c(0.5)

## set up grid of parameters
param_grid <- expand.grid(
    regCoeff=regCoeffs,
    propImp=propImps,
    propLowMedHigh = propLowMedHighs,
    n=ns,
    p=ps
)

probMatrix <- rbind(
    c(0.5, 0.5),
    c(0,1),
    c(1,0)
)

for (i in seq_len(nrow(param_grid))){
# for (i in seq_len(1)){
    set.seed(i)
    cat("Working on", i, "th parameter combination\n")
    current_dynamic_args <- param_grid[i,]

    filename <- paste("res/", simname, i, ".csv", sep="")

    replicate(nreps_per_combo, 
        one_trial(
        n=current_dynamic_args$n,
        p=current_dynamic_args$p,
        filename,
        k=1,
        K=2,
        propImp=current_dynamic_args$propImp, 
        eps=eps, 
        L=15,
        sig_strength=current_dynamic_args$regCoeff,
        propLowMedHigh = probMatrix[current_dynamic_args$propLowMedHigh,],
        verbose=TRUE
    ))
}
