#!/bin/bash
#$ -cwd

nreps=50
K=2
m=NULL
J=5
L=50

# NOTE: If m is NULL, do not include "--m $m" in the line
# Rscript main_cluster.R --simname n_200_p_100_nreps_${nreps}_K_${K}_m_${m}_J_${J}_L_${L} --nreps $nreps --K $K --m $m --J $J --L $L
Rscript main_cluster.R --simname n_200_p_100_nreps_${nreps}_K_${K}_m_${m}_J_${J}_L_${L}_dummy --nreps $nreps --K $K --J $J --L $L