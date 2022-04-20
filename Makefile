NUM_THREADS ?= 0

help:
	@echo Makefile to create dated tree sequences used in paper

all: 


# Save all intermediate files
.SECONDARY:

# Allow filtering in prerequisites
.SECONDEXPANSION:

####################################################
# Standard pipeline for samples file to .dated.trees
####################################################


#############################################
# Ancestral states from Ensembl
#############################################

# HGDP is in GRCh38, and tgp has a GRCh38 liftover available. Others we can lift over. 
# So we download the ancestral states for GRCh38. 

# Recorded in the sample file provenance.
REFERENCE_NAME=GRCh38

ANCESTRAL_STATES_PREFIX=homo_sapiens_ancestor_GRCh38
ANCESTRAL_STATES_TARBALL=${ANCESTRAL_STATES_PREFIX}.tar.gz
ANCESTRAL_STATES_URL=ftp://ftp.ensembl.org/pub/release-100/fasta/ancestral_alleles/${ANCESTRAL_STATES_TARBALL}

${ANCESTRAL_STATES_TARBALL}:
		curl ${ANCESTRAL_STATES_URL} -o ${ANCESTRAL_STATES_TARBALL}

${ANCESTRAL_STATES_PREFIX}/README: ${ANCESTRAL_STATES_TARBALL}
		rm -fR ${ANCESTRAL_STATES_PREFIX}
		tar -xvzf ${ANCESTRAL_STATES_TARBALL}
		# Update access times or we'll keep rebuilding this rule. Have to make sure 
		# that the README we touch is older than the actual fa files.
		touch $@
		touch ${ANCESTRAL_STATES_PREFIX}/*.fa

chr%_ancestral_states.fa: ${ANCESTRAL_STATES_PREFIX}/README
		ln -sf ${ANCESTRAL_STATES_PREFIX}/homo_sapiens_ancestor_$*.fa $@

chr%_ancestral_states.fa.fai: chr%_ancestral_states.fa
		samtools faidx $^

####################
# NYGC Data
####################

nygc_samples.ped:
	curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped -o $@

nygc_phased_genotypes_%.vcf.gz:
	curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_$*.filtered.shapeit2-duohmm-phased.vcf.gz -o $@

nygc_phased_genotypes_%.vcf.gz.tbi: nygc_phased_genotypes_%.vcf.gz
	/home/unix/awohns/tools/bin/tabix -f -p vcf $^

nygc_%.samples: nygc_phased_genotypes_%.vcf.gz.tbi %_ancestral_states.fa.fai nygc_samples.ped
	python3 convert.py 1kg -p \
			nygc_phased_genotypes_$*.vcf.gz \
			$*_ancestral_states.fa \
			-m nygc_samples.ped \
			--ancestral-states-url=${ANCESTRAL_STATES_URL} \
			--reference-name=${REFERENCE_NAME} \
			--num-threads=${NUM_THREADS} \
			$@ > $@.report

trios_info.txt:
	curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index -o $@

nygc_notrios_%.samples: nygc_%.samples
	python3 tsutil.py remove_trios $^ trios_info.txt $@

nygc_removedtrios_%.trees: nygc_%.trees
	python3 tsutil.py remove_trios_ts $^ trios_info.txt $@

####################################################
# Standard pipeline for samples file to .dated.trees
####################################################

%.ancestors: %.samples
	python3 -m tsinfer ga -vp -t ${NUM_THREADS} $^

%.ancestors.trees: %.ancestors
	python3 -m tsinfer ma -vp -t ${NUM_THREADS} $*.samples

%.nosimplify.trees: %.ancestors.trees
	python3 -m tsinfer ms -vp -t ${NUM_THREADS} $*.samples -O $@ --no-simplify

%.nosimplify.nopc.trees: %.ancestors.trees
	python3 -m tsinfer ms -vp -t ${NUM_THREADS} $*.samples -O $@ --no-simplify --no-path-compression

%.trees: %.nosimplify.trees
	python3 tsutil.py simplify $^ $@

%.trees.gz: %.trees
	gzip -c $^ > $@

%.trees.tsz: %.trees
	tszip -k $^ 

%.trees.bcf: %.trees
	msp vcf -P 2 $^ | bcftools view - -O b -o $@

%.snipped.trees: %.trees ${CENTROMERES_CSV}
	python3 tsutil.py snip-centromere $< $@ $* ${CENTROMERES_CSV}


###################################
# List of SNPs in the tree sequence
###################################

# The NYGC vcf doesn't have rsids, so 
bed_chr_%.bed.gz:
	curl ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/bed_chr_$*.bed.gz -o bed_chr_$*.bed.gz

chr%_snplist.txt: bed_chr_%.bed.gz
	python find_snplist.py $*

