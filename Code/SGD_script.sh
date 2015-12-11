#!/bin/bash

#run this with qsub script.sh

#use current directory
#$ -cwd

#name job
#$ -N SGD 

#allow use of environment variables, like $HOME
#$ -V

#range of indices. 0 can't be an index
#$ -t 1-100

#max number of concurrent jobs (be a good citizen). Current Stanford limit is 3000
#$ -tc 1000

#combine output and error into output file, to reduce clutter
#$ -j y

INDEX=$((SGE_TASK_ID-1))
X_INDEX = $((INDEX % 10))
N_INDEX = $((INDEX / 10))

X_ARRAY=(0.01 0.02 0.05 0.1 0.2 0.5 1.0 2.0 5.0 10.0)
N_ARRAY=(10 20 30 50 75 100 125 150 175 200)
python SGD_neural.py -e 200 -i 10 -x ${X_ARRAY[$X_INDEX]} -n ${N_ARRAY[$N_INDEX]}
