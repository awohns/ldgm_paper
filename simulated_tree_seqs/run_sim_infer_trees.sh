#!/bin/bash

# UGER script to allocate resources to run sim_infer_trees.py
# Part of analysis shown in figure 3d

#############################
### Default UGER Requests ###
#############################

# This section specifies uger requests.
# This is good for jobs you need to run multiple times so you don't forget what it needs.

# Name for job
#$ -N simulations_inf

# Memory request for 2G
#$ -l h_vmem=2G

# Cores
#$ -pe smp 10 -binding linear:10

# Runtime request.  Usually 30 minutes is plenty for me and helps me get backfilled into reserved slots.
#$ -l h_rt=10:00:00

# I don't like the top level of my homedir filling up.
#$ -o out/
#$ -e err/

#task array
#$ -t 1-10

######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse

# Use your dotkit
use Anaconda
source activate ldgm_paper

echo "starting job" $SGE_TASK_ID
cd ldgm_paper/simulated_tree_seqs/

python sim_infer_trees.py --seed $SGE_TASK_ID
