#!/bin/bash
cd ~/TP-mardi

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

ref_genome=Homo_sapiens.Chr20.fa
pedigree=./resources/20130606_g1k.ped

# Perform joint variant calling
# Command: gatk GenotypeGVCFs
# Input : list of genomic variant calling files (.g.vcf) + reference genome (.fa)
# Output: Variant calling file (.vcf)
java -jar ${GATK} -T GenotypeGVCFs \
	-R ${ref_genome} \
	--variant father.g.vcf \
	--variant mother.g.vcf \
	--variant daughter.g.vcf \
	-o trio.vcf


##########################
### Recalibrate variants #
##########################

# Not possible here because exomes
# => Skip this part


########################################
### Post-processing: Analysis of trios #
########################################

# Modify pedigree to keep only first 6 columns and no header
cut -f 1-6 ${pedigree} |sed '1,1d' > ${pedigree}.txt

# Phase trios
# Command: gatk PhaseByTransmission
# Input: Variant calling file (.vcf) + reference genome (.fa) + pedigree file (.ped)
# Output: Phased variant calling file (.vcf)
java -jar ${GATK} -T PhaseByTransmission \
	-R ${ref_genome} \
	--variant trio.vcf \
	-ped ${pedigree}.txt \
	-o trio.phased.vcf

# Evaluate variants by computing control metrics
# Command: gatk VariantEval
# Input: List of variant calling files to evaluate (.vcf) + reference genome (.fa)
# Output: Variant evaluation file (.txt)
java -jar ${GATK} -T VariantEval \
	-R ${ref_genome} \
	--eval:set1 trio.vcf \
	--eval:set2 trio.phased.vcf \
	-o trio.phased.VE.txt

# Tabulate the number of sites which overlap and share alleles
# Command: gatk GenotypeConcordance
# Input: 2 variant callings files (.vcf) + reference genome (.fa)
# Output: Report file (.txt)
java -jar ${GATK} -T GenotypeConcordance \
		   -R ${ref_genome} \
		   -eval trio.phased.vcf \
		   -comp trio.vcf \
		   -o trio.GC.txt



