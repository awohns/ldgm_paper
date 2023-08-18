#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# This section specifies uger requests.
# This is good for jobs you need to run multiple times so you don't forget what it needs.

# Name for job
#$ -N LDGM_chr22

# Memory request for 2G
#$ -l h_vmem=1G#25G

# Cores
#$ -pe smp 5 -binding linear:5

# Runtime request.  Usually 30 minutes is plenty for me and helps me get backfilled into reserved slots.
#$ -l h_rt=00:50:00

# I don't like the top level of my homedir filling up.
#$ -o /broad/oconnor/trees/final_ldgm_paper/ldgm_paper/inferred_ldgms/out/
#$ -e /broad/oconnor/trees/final_ldgm_paper/ldgm_paper/inferred_ldgms/err/

#task array
#$ -t 1-23

######################
### Dotkit section ###
######################

chr="21"
treepath="/broad/oconnor/trees/final_ldgm_paper/ldgm_paper/tree_seqs/1kg_chr"$chr".trees"
bedpath="EUR_LD_blocks.chr"$chr".bed"
outpath="."

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse

# Use your dotkit
use Anaconda
source activate ldgm_paper_test

echo "starting job" $SGE_TASK_ID
cd /broad/oconnor/trees/final_ldgm_paper/ldgm_paper/inferred_ldgms
python make_ldgm.py $chr ALL $bedpath $SGE_TASK_ID $treepath $outpath --path-threshold 8 --prune-snps 0.01 --recombination-threshold 0.01 --progress --num-processes 5

