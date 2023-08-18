NUM_THREADS ?= 0

help:
	@echo Makefile to create dated tree sequences used in paper

all: 


# Save all intermediate files
.SECONDARY:

# Allow filtering in prerequisites
.SECONDEXPANSION:

######################################################
# Create LD blocks so all variants fall within a block
######################################################

EUR_LD_blocks.bed:
	curl https://raw.githubusercontent.com/jmacdon/LDblocks_GRCh38/70ad21e53e85543b6fa85cbb517e6ee5ecf8d359/data/EUR_LD_blocks.bed -o $@
	tr -s ' ' '\t' < EUR_LD_blocks.bed > EUR_LD_blocks.tab.bed
	mv EUR_LD_blocks.tab.bed EUR_LD_blocks.bed

EUR_LD_blocks.chr.bed: EUR_LD_blocks.bed
	awk 'BEGIN{FS="\t";OFS="\t"}NR>1{print > "EUR_LD_blocks."$$1".bed"}' $^

###################################
# List of SNPs in the tree sequence
###################################

# The NYGC vcf doesn't have rsids, so 
bed_chr_%.bed.gz:
	curl ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/bed_chr_$*.bed.gz -o bed_chr_$*.bed.gz

chr%_snplist.txt: bed_chr_%.bed.gz
	python find_snplist.py $*

