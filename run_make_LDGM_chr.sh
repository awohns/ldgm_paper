#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# This section specifies uger requests.
# This is good for jobs you need to run multiple times so you don't forget what it needs.

# Name for job
#$ -N LDGM_ALL_removedtrios_22

# Memory request for 1G
#$ -l h_vmem=10G

# Cores
#$ -pe smp 1
#$ -binding linear:1

# Runtime request.  Usually 30 minutes is plenty for me and helps me get backfilled into reserved slots.
#$ -l h_rt=24:00:00

# I don't like the top level of my homedir filling up.
#$ -o /broad/oconnor/trees/nygc/out_ldgm/
#$ -e /broad/oconnor/trees/nygc/err_ldgm/

#task array
#$ -t 1-22

######################
### Dotkit section ###
######################

chr="22"
treepath="/broad/oconnor/trees/nygc/nygc_removedtrios_chr"$chr".trees"
bedpath="/broad/oconnor/trees/ld_intervals_GRCh38/fourier_ls-chr"$chr".GRCh38.bed "
outpath="/broad/oconnor/trees/chr22exploratory/LDGMs/ALL_notrios_softmin/"

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse

# Use your dotkit
use Anaconda
source activate /home/unix/awohns/ldgm

echo "starting job" $SGE_TASK_ID
python /broad/oconnor/trees/nygc/make_ldgm.py $chr ALL $bedpath $SGE_TASK_ID $treepath $outpath --path-threshold 8 --prune-snps 0.01 --recombination-threshold 0.01 --softmin

