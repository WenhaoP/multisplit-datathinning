#!/bin/bash
#$ -cwd

rm -rf ./res/*.csv

Rscript main_cluster.R --simname n_200_p_100 --nreps 100
