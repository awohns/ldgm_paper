# Code for "Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies"

This repository hosts the source code and resources to reproduce the analysis and figures in the paper "Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies". The paper can be cited as follows:

> Pouria Salehi Nowbandegani, Anthony Wilder Wohns, Jenna L. Ballard, Eric S. Lander, Alex Bloemendal, Benjamin M. Neale, and Luke J. O'Connor (2022) _Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies_. Nat Genet. DOI: 10.1038/s41588-023-01487-8


### Respository Contents:
* a pipeline to generate inferred tree sequences from high coverage 1000 Genomes data
* scripts to generate LDGMs from the 1000 Genomes tree sequences and subsequently infer LDGM precision matrices from LDGMs
* code to run analyses in the paper


### Getting Started
First clone this repo:

```
git clone https://github.com/awohns/ldgm_paper.git
```

### Downloading LDGMs inferred from the 1000 Genomes Project

The `data` directory can be populated with data from this [Zenodo repository](https://zenodo.org/record/8157131) containing the precision matrices, LDGMs, tree sequences, correlation matrices, and genotypes used and produced by this study. 

Most users will only want to download precision matrices relevant to the ancestry group(s) they are working with. These can be downloaded as follows:

```
make continent 
```

Where AFR, AMR, EAS, EUR, or SAS can be specified as "continent".

To download all of this data (Warning, this is 43 GB) run:

```
cd data
make download_all
```

### Required Software

#### Installing required python modules

```
python -m pip install -r requirements.txt
```

#### LDGM MATLAB Software

Please clone the [`ldgm`](https://github.com/awohns/ldgm) Github repository in the home directory of this repo:

```
git clone https://github.com/awohns/ldgm.git
```

Many of the analyses require the MATLAB scripts contained in this directory (the `ldgm` python package is installed via the requirements.txt file).

#### Required software for preparing real data

Additionally, we require [BCFtools, SAMtools, and HTSlib](http://www.htslib.org/download/).
to prepare the variant data files from the 1000 Genomes Project.


### Inferring Tree Sequences from Real Data

To infer tree sequences from 1000 Genomes data, navigate to the `tree_seqs` directory and run the following for each chromosome (replacing "chr1" with the chromosome you wish to infer):

```
make 1kg_chr1.trees
```

This command runs the following process:
1. Download NYGC high-coverage 1000 Genomes VCF data
2. Create tsinfer SampleData files from this data (excluding trios)
3. Infer a tree sequence from this data with all biallelic SNPs, mapping indels to the tree sequence with parsimony


### Inferring LDGMs and LDGM precision matrices

1. After inferring tree sequences, navigate to the `inferred_ldgms` directory to generate LDGMs and LDGM precision matrices.

2. Run the following to generate bed files of the LD blocks (this will create one bed file per chromosome).

```
make EUR_LD_blocks.chr.bed
```

3. Next, `run_make_LDGM_chr.sh` contains code to create LDGMs from the inferred tree sequences. This bash script is written for the Broad Institute's internal UGER research computing cluster and is simply used to allocate resources to run `make_LDGMs.py` for every LD block. The script is designed to be run for each chromosome: currently chromsome 21 is specified, but any autosome can be run by modifying the value of `CHR=` and modifying the task list to reflect the number of blocks in that chromsome. You can modify the script to run in your preferred environment.

4. After creating LDGMs with this code, create LDGM precision matrices by running the code in 'run_make_precision_matrices.sh`. You will need to change the commented line of code to point towards your local version of the `ldgm` code repository. This script is also written for the Broad Institute's internal UGER research computing cluster, but can be modified for your preferred environment: this file is simply allocating resources to run the [estimatePrecision function from `ldgm`](https://github.com/awohns/ldgm/blob/42719a699097a7bb08a66b548a96243d1499d129/MATLAB/precision/estimatePrecision.m).


### Running analyses

The `MATLAB` directory to run analyses in the paper.

The analysis shown in Figure 3d can be run from the `simulated_tree_seqs` directory: navigate to the `simulated_tree_seqs` directory and run the code in `run_sim_infer_trees.sh` to simulate and infer tree sequences and then `run_precision_sim.sh` to infer LDGMs and LDGM precision matrices as well as calculate error metrics. As before, these bash scripts allocate resources in the Broad research computing cluster to run the computationally demanding Python and MATLAB software.

### Support
Please feel free to create an issue or contact the paper's corresponding authors with any questions.
