## -----------------------------------------
## Original credit: https://github.com/anna-neufeld/countsplit_paper
## -----------------------------------------

## -----------------------------------------
## install packages
## -----------------------------------------

## define packages to install
packages <- c("argparse", "tidyr", "findpython", "mclust", "fdrtool")

## install all packages that are not already installed
user_name = "wenhaop"
lib_dir = paste("/home/users/", user_name, "/R_lib", sep="")
install.packages(
    setdiff(packages, rownames(installed.packages(lib_dir))),
    repos = "http://cran.us.r-project.org",
    lib=lib_dir
) # install packages only once to avoid errors when parallel computing

## -----------------------------------------
## load user-defined functions, packages
## -----------------------------------------

lapply(packages, library, character.only=TRUE, lib="/home/users/wenhaop/R_lib")

## Needed functions.
source("run_sims.R")

## -----------------------------------------
## load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--simname", default = "robust_ses",
                    help = "name of simulation")
parser$add_argument("--nreps", type = "double", default = 200,
                    help = "number of replicates for each set of params")
parser$add_argument("--K", type = "integer", default = 2,
                    help = "number of data thinning folds")
parser$add_argument("--m", type = "integer", default = NULL,
                    help = "subsample size")
parser$add_argument("--J", type = "integer", default = 5,
                    help = "number of collections")
parser$add_argument("--L", type = "integer", default = 50,
                    help = "number of splits")
args <- parser$parse_args()

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

## extract the hyperparameters
nreps_per_combo <- args$nreps # number of monte-carlo iterations per job
K <- args$K # number of data thinning folds
m <- args$m # subsample size 
J <- args$J # number of collections
L <- args$L # number of splits

## set up grid of parameters
param_grid <- expand.grid(
    regCoeff=regCoeffs,
    propImp=propImps,
    propLowMedHigh = propLowMedHighs,
    n=ns,
    p=ps
)

probMatrix <- rbind(
    c(0.5,0.5),
    c(0,1),
    c(1,0)
)

## -----------------------------------------
## get current dynamic arguments
## -----------------------------------------
## get job id from scheduler

jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## current args
current_dynamic_args <- param_grid[jobid, ]

## -----------------------------------------
## run the simulation nreps_per_job times
## -----------------------------------------
current_seed <- jobid
set.seed(current_seed)


filename <- paste("res/", 
    args$simname, "_", 
    jobid, ".csv", sep="")

## epsilon candidates
# eps=c(0.5)
eps=c(0.1,0.25,0.5,0.75,0.9)
 
system.time(replicate(nreps_per_combo, 
    one_trial(
        n=current_dynamic_args$n,
        p=current_dynamic_args$p,
        filename,
        k=1,
        K=K,
        propImp=current_dynamic_args$propImp, 
        eps=eps, 
        m=m,
        J=J,
        L=L,
        reject.twoside=FALSE,
        sig_strength=current_dynamic_args$regCoeff,
        propLowMedHigh = probMatrix[current_dynamic_args$propLowMedHigh,],
        verbose=TRUE
    )))
