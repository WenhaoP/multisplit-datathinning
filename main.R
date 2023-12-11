## -----------------------------------------
## Original credit: https://github.com/anna-neufeld/countsplit_paper
## -----------------------------------------


## -----------------------------------------
## load user-defined functions, packages
## -----------------------------------------

## tidy stuff
library("tidyr")
library("fdrtool")
library("reticulate")
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
ps <- c(100)

## number of monte-carlo iterations per job
nreps_per_combo <- 10 # number of monte-carlo iterations per job
K <- 2 # number of data thinning folds
m <- NULL # subsample size 
J <- 5 # number of collections
L <- 50 # number of splits

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

dir.name <- paste("res/", simname)

if (!file.exists(dir.name)) {dir.create(dir.name)}

for (i in seq_len(nrow(param_grid))){
# for (i in seq_len(1)){
    set.seed(i)
    cat("Working on", i, "th parameter combination\n")
    current_dynamic_args <- param_grid[i,]

    filename <- paste(dir.name, "/", i, ".csv", sep="")

    for (j in seq_len(nreps_per_combo)) {

        foldername <- paste(dir.name, "/jobid_", i, sep="")
        if (!dir.exists(foldername)) {
            dir.create(foldername)
        }

        foldername <- paste(foldername, "/iter_", j, sep="")
        if (!dir.exists(foldername)) {
            dir.create(foldername)
        }

        one_trial(
            n=current_dynamic_args$n,
            p=current_dynamic_args$p,
            filename,
            foldername,
            k=1,
            K=K,
            propImp=current_dynamic_args$propImp, 
            eps=eps, 
            m=m,
            J=J,
            L=L,
            reject.twoside=TRUE,
            sig_strength=current_dynamic_args$regCoeff,
            propLowMedHigh = probMatrix[current_dynamic_args$propLowMedHigh,],
            verbose=TRUE
        )
    }
}
