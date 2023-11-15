#!/bin/bash
#$ -cwd

rm -rf ./res/*.csv

Rscript main_cluster.R --simname test --nreps 100
