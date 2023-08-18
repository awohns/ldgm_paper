# Code for "Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies"

This repo contains code to reproduce analysis and figures in "Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies". The preprint can be cited as follows:

> Pouria Salehi Nowbandegani, Anthony Wilder Wohns, Jenna L. Ballard, Eric S. Lander, Alex Bloemendal, Benjamin M. Neale, and Luke J. O'Connor (2022) _Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies_. Nat Genet. DOI: 10.1038/s41588-023-01487-8


The repo includes:
* a pipeline to generate inferred tree sequences from high coverage 1000 Genomes data
* scripts to generate LDGMs and infer LDGM precision matrices from 1000 Genomes tree sequences
* code to run analyses in the paper

### Getting Started
First clone the repo:

```
git clone https://github.com/awohns/ldgm_paper.git
```

### Downloading LDGMs inferred from the 1000 Genomes Project

The `data` directory can be populated with data from this [Zenodo repository](https://zenodo.org/record/8157131) containing the precision matrices, LDGMs, tree sequences, correlation matrices, and genotypes used and produced by this study. To download all of this data (Warning! This is 43 GB) run:

```
cd data
make download_all
```

Most users will only want to download precision matrices relevant to the ancestry group(s) they are working with. These can be downloaded as follows:

```
make continent 
```

Where AFR, AMR, EAS, EUR, or SAS can be specified as "continent".

### Required Software

#### Installing required python modules

```
python -m pip install -r requirements.txt
```


#### Required software for preparing real data

Please clone the `ldgm` Github repository: many of the following analyses require the MATLAB scripts contained in this directory (the python package is installed in the previous command).
We require [BCFtools, SAMtools, and HTSlib](http://www.htslib.org/download/), as well as
[Picard](https://broadinstitute.github.io/picard/),
to prepare the variant data files from
the 1000 Genomes Project.


### Inferring Tree Sequences from Real Data

To infer tree sequences from 1000 genomes data, navigate to the `tree_seqs` directory and run the following for each chromosome:

```
make nygc_notrios_chr1.trees
```

This command runs the following process:
1. Download NYGC high-coverage 1000 Genomes VCF data
2. Create tsinfer SampleData files from this data (excluding trios)
3. Infer a tree sequence from this data with all biallelic SNPs, mapping indels to the tree sequence with parsimony

Additionally, run the following (for any chromosome) to generate bed files of the LD blocks:

```
make EUR_LD_blocks.chr1.bed
```


### Inferring LDGMs and LDGM precision matrices

After inferring tree sequences, navigate to the `inferred_ldgms` directory to generate LDGMs and LDGM precision matrices.
`run_make_LDGM_chr.sh` contains code to create LDGMs from the inferred tree sequences. This bash script is written for the Broad Institute's internal UGER research computing system and is simply used to allocate resources to run `make_LDGMs.py` for every LD block. The script is designed to be run for each chromosome: currently chromsome 21 is specified, but any autosome can be run by modifying the value of `CHR=` and modifying the task list to reflect the number of blocks in that chromsome. You can modify the script to run in your preferred environment. Further instructions can be found in the `README` in this directory.

After creating LDGMs with this code, create LDGM precision matrices by running the code in 'run_make_precision_matrices.sh`. You will need to change the commented line of code to point towards your local version of the `ldgm` code repository. This script is also written for the Broad Institute's internal UGER research computing system but can be modified for your preferred environment.


### Running analyses

Refer to the code in the `MATLAB` directory to run analyses in the paper.
For the analysis shown in Figure 3d; this can be run from the `simulated_tree_seqs` directory: navigate to the `simulated_tree_seqs` directory and run the code in `run_sim_infer_trees.sh` to simulate and infer tree sequences and then `run_precision_sim.sh` to infer LDGMs and LDGM precision matrices as well as calculate error metrics.

Please feel free to create an issue or contact the corresponding authors with any questions!
