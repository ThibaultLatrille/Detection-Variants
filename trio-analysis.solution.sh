#!/bin/bash

WORK_DIR=~/variant_calling
# If on Biosphere (IFB), use '/mnt/data/variant_calling' instead of '~/variant_calling' as a working directory
# WORK_DIR=/mnt/data/variant_calling

########################################################################################################################
# Requirements:
#   Java (version 8)
#   GATK (version 3.3)
########################################################################################################################

java -version
java -jar ${GATK} --help
java -jar ${PICARD}

############################
### Joint variant calling ##
############################

REF_GENOME=Homo_sapiens.Chr20.fa
PEDIGREE=20130606_g1k.ped

# Perform joint variant calling
# Command: gatk GenotypeGVCFs
# Input : list of genomic variant calling files (.g.vcf) + reference genome (.fa)
# Output: Variant calling file (.vcf)

java -jar ${GATK} -T GenotypeGVCFs \
	-R ${REF_GENOME} \
	--variant HG02024.g.vcf \
	--variant HG02025.g.vcf \
	--variant HG02026.g.vcf \
	-o trio.vcf


##########################
### Recalibrate variants #
##########################

# Not possible here because exomes
# => Skip this part


########################################
### Post-processing: Analysis of trios #
########################################

# Modify PEDIGREE to keep only first 6 columns and no header
cut -f 1-6 ${PEDIGREE} |sed '1,1d' > ${PEDIGREE}.txt

# Phase trios
# Command: gatk PhaseByTransmission
# Input: Variant calling file (.vcf) + reference genome (.fa) + PEDIGREE file (.ped)
# Output: Phased variant calling file (.vcf)
java -jar ${GATK} -T PhaseByTransmission \
	-R ${REF_GENOME} \
	--variant trio.vcf \
	-ped ${PEDIGREE}.txt \
	-o trio.phased.vcf

# Evaluate variants by computing control metrics
# Command: gatk VariantEval
# Input: List of variant calling files to evaluate (.vcf) + reference genome (.fa)
# Output: Variant evaluation file (.txt)
java -jar ${GATK} -T VariantEval \
	-R ${REF_GENOME} \
	--eval:set1 trio.vcf \
	--eval:set2 trio.phased.vcf \
	-o trio.phased.VE.txt

# Tabulate the number of sites which overlap and share alleles
# Command: gatk GenotypeConcordance
# Input: 2 variant callings files (.vcf) + reference genome (.fa)
# Output: Report file (.txt)
java -jar ${GATK} -T GenotypeConcordance \
		   -R ${REF_GENOME} \
		   -eval trio.phased.vcf \
		   -comp trio.vcf \
		   -o trio.GC.txt
