#!/bin/bash

njobs=45

sbatch -p witten-12c128g -t 72:00:00 --array=1-$njobs -e ./out/s-%A_%a.out -o ./out/s-%A_%a.out ./call_sim.sh 


