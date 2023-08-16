#!/bin/bash

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
#$ -o /broad/oconnor/trees/final_ldgm_paper/ldgm_paper/simulated_tree_seqs/out/
#$ -e /broad/oconnor/trees/final_ldgm_paper/ldgm_paper/simulated_tree_seqs/err/

#task array
#$ -t 1-10

######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse

# Use your dotkit
use Anaconda
source activate ldgm_paper_test

echo "starting job" $SGE_TASK_ID
cd /broad/oconnor/trees/final_ldgm_paper/ldgm_paper/simulated_tree_seqs/


python sim_infer_trees.py --seed $SGE_TASK_ID
