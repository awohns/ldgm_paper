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
# 1kg NYGC Data
####################

1kg_samples.ped:
	curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped -o $@

1kg_phased_genotypes_%.vcf.gz:
	curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_$*.filtered.shapeit2-duohmm-phased.vcf.gz -o $@

1kg_phased_genotypes_%.vcf.gz.tbi: 1kg_phased_genotypes_%.vcf.gz
	/home/unix/awohns/tools/bin/tabix -f -p vcf $^

1kg_inctrios_%.samples: 1kg_phased_genotypes_%.vcf.gz.tbi %_ancestral_states.fa.fai 1kg_samples.ped
	python3 convert.py 1kg -p \
			1kg_phased_genotypes_$*.vcf.gz \
			$*_ancestral_states.fa \
			-m 1kg_samples.ped \
			--ancestral-states-url=${ANCESTRAL_STATES_URL} \
			--reference-name=${REFERENCE_NAME} \
			--num-threads=${NUM_THREADS} \
			$@ > $@.report

trios_info.txt:
	curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index -o $@

1kg_%.samples: 1kg_inctrios_%.samples trios_info.txt
	python3 tsutil.py remove_trios $^ trios_info.txt $@


####################################################
# Standard pipeline for samples file to .trees
####################################################

1kg_%.trees: 1kg_%.samples
	python3 run_tsinfer_excluded.py $^ ${NUM_THREADS}


