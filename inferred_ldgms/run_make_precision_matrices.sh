#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# This section specifies uger requests.
# This is good for jobs you need to run multiple times so you don't forget what it needs.

# Name for job
#$ -N make_precision_matrices

# Memory request for 2G
#$ -l h_vmem=12G

# Cores
#$ -pe smp 1
#$ -binding linear:1

# Runtime request.  Usually 30 minutes is plenty for me and helps me get backfilled into reserved slots.
#$ -l h_rt=150:00:00

# I don't like the top level of my homedir filling up.
#$ -o /broad/oconnor/trees/final_ldgm_paper/ldgm_paper/inferred_ldgms/out/
#$ -e /broad/oconnor/trees/final_ldgm_paper/ldgm_paper/inferred_ldgms/err/

#task array
#$ -t 1-4

######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse

# Use your dotkit
use Matlab

cd /../ldgm/MATLAB/precision
matlab -r "estimatePrecision('/broad/oconnor/trees/final_ldgm_paper/ldgm_paper/inferred_ldgms/',\
    'genos_file_index',$SGE_TASK_ID,\
    'output_dir','precisionMatrices/',\
    'population_data_file','/broad/oconnor/trees/final_ldgm_paper/ldgm_paper/data/1kg_nygc_trios_removed_All_pops_geno_ids_pops.csv',\
    'population_name','AFR',\
    'path_distance_threshold',4,\
    'l1_penalty',0.10,\
    'minimum_maf',0.01,\
    'unpenalized_iterations', 25,\
    'penalized_iterations', 5,\
    'banded_control_ldgm',0);"

