# Code for "Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies"

This repo contains code to reproduce analysis and figures in "Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies". The preprint can be cited as follows:

> Pouria Salehi Nowbandegani, Anthony Wilder Wohns, Jenna L. Ballard, Eric S. Lander, Alex Bloemendal, Benjamin M. Neale, and Luke J. O’Connor (2022) _Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies_. Biorxiv doi: https://doi.org/10.1101/2022.09.06.506858


The repo includes:
* a pipeline to generate an inferred tree sequence of high coverage 1000 Genomes Data
* scripts to generate LDGMs and infer LDGM precision matrices from the 1000 genomes tree sequence
* code to run analyses in the paper
* code to produce non-schematic figures in the paper

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

We require [BCFtools, SAMtools, and HTSlib](http://www.htslib.org/download/), as well as
[Picard](https://broadinstitute.github.io/picard/),
to prepare the variant data files from
the 1000 Genomes Project.


### Inferring Tree Sequences from Real Data

To infer tree sequences from 1000 genomes data, run the following for each chromosome:

```
cd tgp
make nygc_notrios_chr1.trees
```


### Inferring LDGMs and LDGM precision matrices

`make_ldgm.py` contains code to create LDGMs from the inferred tree sequences.


### Running validation analyses in MATLAB

Refer to the code in the `MATLAB` directory to run analyses (save for the analysis shown in Figure 3d) and make all plots in the paper.


