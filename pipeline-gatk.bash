#!/bin/bash
cd ~/Shared/Detection-Variants/data
source ~/.profile


########################################################################################################################
# Requirements:
#   Java (version 8)
#   GATK (version 3.3)
########################################################################################################################

java -version
ls $GATK
ls $PICARDTOOL


###################################################
### Overall workflow of variant-calling with GATK #
###################################################

# GATK overall workflow:
# 1. Prepare GATK input data (pre-processing)
# 2. Local realignment around indels
# 3. Base quality recalibration
# 4. Call variants
# 5. Recalibrate variants
# 6. Post-processing


#####################
# Extract resources #
#####################

# Make directory for resources
mkdir ./resources/
cd ./resources/

# Download resources
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed variant calling files (.vcf.gz)

# Known indels
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

# Known SNPs
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz.tbi

# Pedigree
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped


# Choose variable names
cd ../
known_indels=./resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
known_SNPs=./resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz
pedigree=./resources/20130606_g1k.ped
ref_genome=Homo_sapiens.Chr20.fa
file_name=father


#######################################
### Prepare reference genome for GATK #
#######################################

# Index reference
# Command: samtools faidx
# Input: reference genome (.fa)
# Output: indexed genome (.fa.fai) (no need to mention it)
samtools faidx ${ref_genome}

# Create a sequence dictionnary for the reference genome
# Command: CreateSequenceDictionary (PICARDtools)
# Input: reference genome (.fa)
# Sequence dictionnary (.dict)
java -jar $PICARDTOOL CreateSequenceDictionary \
	R=${ref_genome} \
	O=${ref_genome/.fa/.dict}


#############################
### Prepare GATK input data #
#############################

# Fix mate pair information (prevention in case of errors) + sort (in case not done)
# Command: FixMateInformation (PICARDtools)
# Input: alignment (.bam)
# Ouput: alignment with fixed mate information (.bam)
java -jar $PICARDTOOL FixMateInformation \
	I=${file_name}.bam \
	O=${file_name}.fixed_mate.bam \
	SO=coordinate \
	VALIDATION_STRINGENCY=LENIENT

# Mark Duplicate reads
# Command: MarkDuplicates (PICARDtools)
# Input: alignment (.bam)
# Ouput: alignment with duplicates marked (.bam)
java -jar $PICARDTOOL MarkDuplicates \
	REMOVE_DUPLICATES=FALSE \
	I=${file_name}.fixed_mate.bam \
	O=${file_name}.marked_dups.bam \
	M=${file_name}.dup_metrics.txt

# Make sure that read group already added. Else, add it now
# Command: gatk AddOrReplaceReadGroups
# Input: alignment (.bam) and read group
# Ouput: alignment (.bam)

# Make sure that index already obtained. Else, do it now
# Command: samtools index / BuilBamIndex (PICARDtools)
# Input: alignment (.bam)
java -jar $PICARDTOOL BuildBamIndex \
	INPUT=${file_name}.marked_dups.bam


#####################################
### Local realignemnt around indels #
#####################################

# Find regions that need to be realigned
# Command: gatk RealignerTargetCreator
# Input: preprocessed alignment (.bam) + compressed known indels (.vcf.gz) + reference genome (.fa)
# Output: list of intervals (.list / .txt)
#gatk RealignerTargetCreator \
java -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator \
	-R ${ref_genome} \
	-known ${known_indels} \
	-I ${file_name}.marked_dups.bam \
	-o ${file_name}.target_intervals.list 

# Perform local realignement
# Command: gatk IndelRealigner
# Input: preprocessed alignment (.bam) + compressed known indels (.vcf.gz) + list of intervals (.list / .txt) + reference genome (.fa)
# Output: realigned alignment (.bam)
#gatk IndelRealigner \
java -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner \
	-R ${ref_genome} \
	-known ${known_indels} \
	-I ${file_name}.marked_dups.bam \
	-targetIntervals ${file_name}.target_intervals.list \
	-o ${file_name}.realigned_reads.bam



################################
### Base quality recalibration #
################################

# Generate a base recalibration table to analyse patterns of covariation in the dataset
# Command: gatk BaseRecalibrator
# Input: realigned alignment (.bam) + compressed known indels (.vcf.gz) + compressed known snps (.vcf.gz) + reference genome (.fa)
# Output: base recalibration table (.table / .txt)
java -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator \
	-R ${ref_genome} \
	-knownSites ${known_SNPs} \
	-knownSites ${known_indels} \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov ContextCovariate \
	-I ${file_name}.realigned_reads.bam \
	-o ${file_name}.recal_data.table

# (Optional)
# Analyse covariation remaining after recalibration (2nd pass)
# Command: gatk BaseRecalibrator
# Input: realigned alignment (.bam) + compressed known indels (.vcf.gz) + compressed known snps (.vcf.gz) + reference genome (.fa) + base quality recalbration table (.table / .txt)
# Output: post recalibration table (.table / .txt)
java -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator \
	-R ${ref_genome} \
	-knownSites ${known_SNPs} \
	-knownSites ${known_indels} \
	-I ${file_name}.realigned_reads.bam \
	-BQSR ${file_name}.recal_data.table \
	-o ${file_name}.post_recal_data.table 

# (Optional)
# Generate before/after plots to assess the recalibration effect
# Command: gatk AnalyzeCovariates
# Input: base quality recalbration table (.table / .txt) + post recalibration table (.table / .txt) + reference genome (.fa)
# Output: plots (.pdf)
java -jar $GATK/GenomeAnalysisTK.jar -T AnalyzeCovariates \
	-R ${ref_genome} \
	-before ${file_name}.recal_data.table \
	-after ${file_name}.post_recal_data.table \
	-plots ${file_name}.recalibration_plots.pdf

# Perform base quality recalibration
# Command: gatk PrintReads + BQSR option
# Input: realigned alignment (.bam) + reference genome (.fa) + base quality recalbration table (.table / .txt)
# Output: base quality recalibrated alignement (.bam)
java -jar $GATK/GenomeAnalysisTK.jar -T PrintReads \
	-R ${ref_genome} \
	-I ${file_name}.realigned_reads.bam \
	-BQSR ${file_name}.recal_data.table \
	-o ${file_name}.recal_reads.bam



###################
### Call variants #
###################

# Perform variant calling
# Command: gatk HaplotypeCaller
# Input: base quality recalibrated alignement (.bam) + reference genome (.fa)
# Output: Genomic variant calling file (.g.vcf)
java -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller \
	-R ${ref_genome} \
	-I ${file_name}.recal_reads.bam \
	-ERC GVCF \
	--genotyping_mode DISCOVERY \
	--variant_index_type LINEAR \
	-variant_index_parameter 128000 \
	-o ${file_name}.g.vcf 


# Redo the whole procedure (starting from pre-processing of data) for daughter and mother files
# ...


# Perform joint variant calling
# Command: gatk GenotypeGVCFs
# Input : list of genomic variant calling files (.g.vcf) + reference genome (.fa)
# Output: Variant calling file (.vcf)
java -jar $GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs \
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
cut -f 1-6 $pedigree |sed '1,1d' > $pedigree.txt

# Phase trios
# Command: gatk PhaseByTransmission
# Input: Variant calling file (.vcf) + reference genome (.fa) + pedigree file (.ped)
# Output: Phased variant calling file (.vcf)
java -jar $GATK/GenomeAnalysisTK.jar -T PhaseByTransmission \
	-R ${ref_genome} \
	--variant trio.vcf \
	-ped ${pedigree}.txt \
	-o trio.phased.vcf

# Evaluate variants by computing control metrics
# Command: gatk VariantEval
# Input: List of variant calling files to evaluate (.vcf) + reference genome (.fa)
# Output: Variant evaluation file (.txt)
java -jar $GATK/GenomeAnalysisTK.jar -T VariantEval \
	-R ${ref_genome} \
	--eval:set1 trio.vcf \
	--eval:set2 trio.phased.vcf \
	-o trio.phased.VE.txt

# Tabulate the number of sites which overlap and share alleles
# Command: gatk GenotypeConcordance
# Input: 2 variant callings files (.vcf) + reference genome (.fa)
# Output: Report file (.txt)
java -jar $GATK/GenomeAnalysisTK.jar -T GenotypeConcordance \
		   -R ${ref_genome} \
		   -eval trio.phased.vcf \
		   -comp trio.vcf \
		   -o trio.GC.txt



