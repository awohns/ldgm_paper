#!/bin/bash

# UGER script to infer LDGM precision matrices from LDGMs created
# from msprime simulations and inferred tree sequences from these simulations.
# Part of the analysis shown in Figure 3d

#############################
### Default UGER Requests ###
#############################

# This section specifies uger requests.
# This is good for jobs you need to run multiple times so you don't forget what it needs.

# Name for job
#$ -N sim

# Memory request for 2G
#$ -l h_vmem=4G

# Cores
#$ -pe smp 1
#$ -binding linear:1

# Runtime request.  Usually 30 minutes is plenty for me and helps me get backfilled into reserved slots.
#$ -l h_rt=82:00:00

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
use Matlab

cd /broad/oconnor/trees/final_ldgm_paper/ldgm_paper/ldgm/MATLAB/precision
# Run LDGM precision matrix inference for LDGMs created from the simulated tree seqs
matlab -r "estimatePrecision('ldgm_paper/simulated_tree_seqs/',\
    'data_type','genotypes',\
    'data_file_index',SGE_TASK_ID,\
    'data_pattern','sim_*_MAF_0.01_RF_0.01_T_8',\
    'data_file_extension','.genos',\
    'output_dir','ldgm_paper/simulated_tree_seqs/',\
    'AF_column_name','0',\
    'path_distance_threshold',4,\
    'l1_penalty',0.10,\
    'minimum_maf',0.01,\
    'unpenalized_iterations', 25,\
    'penalized_iterations', 5,\
    'banded_control_ldgm',0);"

# Run LDGM precision matrix inference for LDGMs created from the inferred tree seqs
matlab -r "estimatePrecision('ldgm_paper/simulated_tree_seqs/',\
    'data_type','genotypes',\
    'data_file_index',SGE_TASK_ID,\
    'data_pattern','inferred_sim_*_MAF_0.01_RF_0.01_T_8',\
    'data_file_extension','.genos',\
    'output_dir','ldgm_paper/simulated_tree_seqs/',\
    'AF_column_name','0',\
    'path_distance_threshold',4,\
    'l1_penalty',0.10,\
    'minimum_maf',0.01,\
    'unpenalized_iterations', 25,\
    'penalized_iterations', 5,\
    'banded_control_ldgm',0);"

