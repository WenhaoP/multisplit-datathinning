## -----------------------------------------
## Original credit: https://github.com/anna-neufeld/countsplit_paper
## -----------------------------------------

## -----------------------------------------
## install packages
## -----------------------------------------

## define packages to install
packages <- c("argparse", "tidyr", "findpython", "mclust")

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

## get command line arguments
library("argparse", lib="/home/users/wenhaop/R_lib")
## tidy stuff
library("tidyr", lib="/home/users/wenhaop/R_lib")
## find an acceptable python binary
library("findpython", lib="/home/users/wenhaop/R_lib")
## Gaussian Mixture Modelling
library("mclust", lib="/home/users/wenhaop/R_lib")
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

## number of monte-carlo iterations per job
nreps_per_combo <- args$nreps

## set up grid of parameters
param_grid <- expand.grid(regCoeff=regCoeffs,
                          propImp=propImps,
                          propLowMedHigh = propLowMedHighs,
                          n=ns,
                          p=ps)

probMatrix <- rbind(c(0.5,0.5),
                    c(0,1),
                    c(1,0))

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

filename <- paste("res/", args$simname, jobid, ".csv", sep="")

## epsilon candidates
# eps=c(0.5)
eps=c(0.1,0.25,0.5,0.75,0.9)
 
system.time(replicate(nreps_per_combo, 
                      one_trial(n=current_dynamic_args$n,
                                            p=current_dynamic_args$p,
                                            filename,
                                            k=1,
                                            K=2,
                                            propImp=current_dynamic_args$propImp, 
                                            eps=eps, 
                                            # L=50,
                                            sig_strength=current_dynamic_args$regCoeff,
                                 propLowMedHigh = probMatrix[current_dynamic_args$propLowMedHigh,])))
