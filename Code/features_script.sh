#!/bin/bash

#run this with qsub script.sh

#use current directory
#$ -cwd

#name job
#$ -N anti-kt

#allow use of environment variables, like $HOME
#$ -V

#range of indices. 0 can't be an index
#$ -t 1-1000

#max number of concurrent jobs (be a good citizen). Current Stanford limit is 3000
#$ -tc 1000

#combine output and error into output file, to reduce clutter
#$ -j y

INDEX=$((SGE_TASK_ID-1))

python particle_features.py -e $INDEX
