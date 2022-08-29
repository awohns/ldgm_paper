# Code for "Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies"

This repo contains code used in "Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies" (see forthcoming preprint).

This includes:
* a pipeline to generate an inferred tree sequence of high coverage 1000 Genomes Data
* scripts to generate LDGMs and infer LDGM precision matrices from the 1000 genomes tree sequence
* an analysis comparing the accuracy of LDGMs from simulated vs. inferred (using [tsinfer](https://tsinfer.readthedocs.io/)) tree sequences
* code to produce all non-schematic figures in the paper


### Getting Started
You must first clone this repo:

```
git clone https://github.com/awohns/ldgm_paper.git
```

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

`make_ldgm.py` contains code to create LDGMs from the inferred tree sequence


### Running validation analyses in MATLAB

...

### Plotting figures

...


