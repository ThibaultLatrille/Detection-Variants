#!/bin/bash
cd ~/Shared/Detection-Variants/data
source ~/.profile

########################################################################################################################
# Requirements:
#	Java (version 8)
#	FastQC (version 0.11.7)
#	BWA-MEM (version 0.7.17-r1194-dirty)
#	SAMtools (version 1.9)
#	IGV (version 2.4.14)
#	GATK (version 4.0.8.1)
########################################################################################################################

java -version
fastqc -version
bwa
samtools
gatk --list

##########################################################
## Download, extract and index the reference chromosome ##
##########################################################

# Download the reference Human chromosome (chromosome 20) from Ensembl
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed reference sequence (.fa.gz)
wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz -O Homo_sapiens.Chr20.fa.gz

# Extract the reference chromosome
# Command: gunzip
# Input: compressed reference sequence (.fa.gz)
# Ouput: reference sequence (remove .gz file)
gunzip Homo_sapiens.Chr20.fa.gz

# Index the reference chromosome
# Command: bwa index
# Input: reference (.fa)
# Ouput: indexed reference (.fa.amb, .fa.ann, .fa.bwt, fa.pac, .fa.sa)
bwa index Homo_sapiens.Chr20.fa

######################################################
## Mapping of a family trio to the reference genome ##
######################################################

# The sequences are from an East Asian (Kinh Vietnamese) family forming a trio : daughter/mother/father
# Data available at http://www.internationalgenome.org/data-portal/sample/HG02024
# Daughter:
#       StudyId: SRP004063
#       SampleName: HG02024
#       Library: Pond-206419
#       ExperimentID: SRX001595
#       RunId: SRR822251
#       PlatformUnit: C1E0PACXX121221.6.tagged_373
#       InstrumentModel: Illumina HiSeq 2000
#       InsertSize: 160
# Mother:
#       StudyId: SRP004063
#       SampleName: HG02025
#       Library: Catch-88584
#       ExperimentID: SRX103760
#       RunId: SRR361100
#       PlatformUnit: BI.PE.110902_SL-HBC_0182_AFCD046MACXX.2.tagged_851
#       InstrumentModel: Illumina HiSeq 2000
#       InsertSize: 96
# Father:
#       StudyId: SRP004063
#       SampleName: HG02026

#############################
## Mapping of the daughter ##
#############################

# Download paired sequencing reads for the daughter (SampleName: HG02024, RunId: SRR822251)
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed sequencing reads (.fastq.gz)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02024/sequence_read/SRR822251_1.filt.fastq.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02024/sequence_read/SRR822251_2.filt.fastq.gz

# Map the paired sequencing reads against the reference Human chromosome 20
# Command: bwa mem
# Input: indexed reference (.fa), and compressed sequencing reads (.fastq.gz)
# Ouput: alignment (.sam)
bwa mem -t 4 -M Homo_sapiens.Chr20.fa SRR822251_1.filt.fastq.gz SRR822251_2.filt.fastq.gz > HG02024_SRR822251.sam

# Compute summary statistics of the alignment
# Command: samtools flagstats
# Input: alignment (.sam)
# Ouput: text file (human and computer readable)
samtools flagstat HG02024_SRR822251.sam > HG02024_SRR822251.sam.flagstats

# Compress the alignment and filter unaligned reads
# Command: samtools view
# Input: alignment (.sam)
# Ouput: compressed alignment (.bam)
samtools view -@ 4 HG02024_SRR822251.sam -Sbh -f 3 > HG02024_SRR822251.bam

# Sort the alignment
# Command: samtools sort
# Input: compressed alignment (.bam)
# Ouput: sorted and compressed alignment (.bam)
samtools sort -@ 4 HG02024_SRR822251.bam > HG02024_SRR822251.sorted.bam

# Add Read group (cf https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups)
# Command: gatk AddOrReplaceReadGroups
# Input: alignment (.bam) and read group (Read group identifier, DNA preparation library identifier, Platform, Platform Unit, Sample)
# Ouput: annotated alignment (.bam)
gatk AddOrReplaceReadGroups -I HG02024_SRR822251.sorted.bam -O daughter.bam \
                            --RGID SRR822251 --RGLB Pond-206419 --RGPL illumina \
                            --RGPU C1E0PACXX121221.6.tagged_373 --RGSM HG02024 --RGPI 160

# Visualize read group
# Command: samtools view -H && grep
# Input: annotated alignment (.bam)
# Ouput: read group
samtools view -H daughter.bam | grep '@RG'

# Compute statistics of the alignment
# Command: samtools-stats
# Input: alignment (.bam)
# Ouput: text file (human and computer readable)
samtools stats daughter.bam > daughter.bam.stats

# Plot statistics of the alignment
# Command: plot-bamstats
# Input: statistics text file (output of samtools-stats)
# Ouput: plots (.png)
plot-bamstats -p ./plots-daughter/ daughter.bam.stats

# Index the alignment
# Command: samtools index
# Input: alignment (.bam)
# Ouput: indexed alignment (.bam.bai)
samtools index daughter.bam


###########################
## Mapping of the mother ##
###########################

# Variables definition
FTP_SEQ_FOLDER=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3 # Ftp folder from 1000Genomes project
RUN_ID=SRR361100 # Read group identifier
SAMPLE_NAME=HG02025 # Sample
INSTRUMENT_PLATFORM=illumina # Platform/technology used to produce the read
LIBRARY_NAME=Catch-88584 # DNA preparation library identifier
RUN_NAME=BI.PE.110902_SL-HBC_0182_AFCD046MACXX.2.tagged_851 # Platform Unit
INSERT_SIZE=96 # Insert size

# Download paired sequencing reads for the mother
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed sequencing reads (.fastq.gz)
wget ${FTP_SEQ_FOLDER}/data/${SAMPLE_NAME}/sequence_read/${RUN_ID}_1.filt.fastq.gz
wget ${FTP_SEQ_FOLDER}/data/${SAMPLE_NAME}/sequence_read/${RUN_ID}_2.filt.fastq.gz

# Map, filter, and sort the paired sequencing reads of the mother against the reference genome
# Command: bwa mem && samtools view && samtools sort
# Input: indexed reference (.fa), and compressed sequencing reads (.fastq.gz)
# Ouput: sorted alignment (.bam)
bwa mem -t 4 -M Homo_sapiens.Chr20.fa ${RUN_ID}_1.filt.fastq.gz ${RUN_ID}_2.filt.fastq.gz | samtools view -@ 4 -bh -f 3 | samtools sort -@ 4 > ${SAMPLE_NAME}.${RUN_ID}.sorted.bam

# Add Read group
# Command: gatk AddOrReplaceReadGroups
# Input: alignment (.bam) and read group
# Ouput: alignment (.bam)
gatk AddOrReplaceReadGroups -I ${SAMPLE_NAME}.${RUN_ID}.sorted.bam -O mother.bam \
                            --RGID ${RUN_ID} --RGLB ${LIBRARY_NAME} --RGPL ${INSTRUMENT_PLATFORM} \
                            --RGPU ${RUN_NAME} --RGSM ${SAMPLE_NAME} --RGPI ${INSERT_SIZE}

# Index the alignment
# Command: samtools index
# Input: alignment (.bam)
# Ouput: indexed alignment (.bam.bai)
samtools index mother.bam

###########################
## Mapping of the father ##
###########################

# Variables definition
SAMPLE_NAME=HG02026 # Sample

# Download index file containing sequencing runs information
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: text file (.index)
wget ${FTP_SEQ_FOLDER}/20130502.phase3.analysis.sequence.index -O 20130502.phase3.index

# Filter paired exome sequencing runs related to father (HG02026)
# Command: grep && grep -v
# Input: tab-separated values file (.index)
# Ouput: filtered comma-separated values file (.index)
grep ${SAMPLE_NAME} 20130502.phase3.index | grep "exome" | grep 'PAIRED' | grep -v 'Solexa' | grep -v 'from blood' | grep -v '_1.filt.fastq.gz' | grep -v '_2.filt.fastq.gz' | sed 's/\t/,/g' > father.index

# File containing the list of alignments (each line is a .bam file)
# This file is necessary to merge multiple alignments into a single alignment.
# Command: touch
# Input: file name
# Ouput: empty file (.bamlist)
touch ${SAMPLE_NAME}.bamlist

# for each sequencing run (the first 10), align to the reference, sort, add read group and index
head -6 father.index | while IFS="," read FASTQ_FILE MD5 RUN_ID STUDY_ID STUDY_NAME CENTER_NAME SUBMISSION_ID SUBMISSION_DATE SAMPLE_ID SAMPLE_NAME POPULATION EXPERIMENT_ID INSTRUMENT_PLATFORM INSTRUMENT_MODEL LIBRARY_NAME RUN_NAME RUN_BLOCK_NAME INSERT_SIZE LIBRARY_LAYOUT PAIRED_FASTQ WITHDRAWN WITHDRAWN_DATE COMMENT READ_COUNT BASE_COUNT ANALYSIS_GROUP
do

    # Variables definition
    FASTQ_FILE_1=${FASTQ_FILE/.filt.fastq.gz/_1.filt.fastq.gz} # Path of the fasta file in the FTP folder
    FASTQ_FILE_2=${FASTQ_FILE/.filt.fastq.gz/_2.filt.fastq.gz} # Path of the fasta file in the FTP folder (pairing file)

    # Download paired sequencing reads for the father
    # Command: wget
    # Input: url (http:// or ftp://)
    # Ouput: compressed sequencing reads (.fastq.gz)
    wget ${FTP_SEQ_FOLDER}/${FASTQ_FILE_1} -O ${SAMPLE_NAME}.${RUN_ID}_1.fastq.gz
    wget ${FTP_SEQ_FOLDER}/${FASTQ_FILE_2} -O ${SAMPLE_NAME}.${RUN_ID}_2.fastq.gz

    # Map, filter, and sort the paired reads of the sequencing run against the reference genome
    # Command: bwa mem && samtools view && samtools sort
    # Input: indexed reference (.fa), and compressed sequencing reads (.fastq.gz)
    # Ouput: sorted alignment (.bam)
    bwa mem -t 4 -M Homo_sapiens.Chr20.fa ${SAMPLE_NAME}.${RUN_ID}_1.fastq.gz ${SAMPLE_NAME}.${RUN_ID}_2.fastq.gz | samtools view -@ 4 -bh -f 3 | samtools sort -@ 4 > ${SAMPLE_NAME}.${RUN_ID}.sorted.bam

    # Add Read group
    # Command: gatk AddOrReplaceReadGroups
    # Input: alignment (.bam) and read group
    # Ouput: alignment (.bam)
    gatk AddOrReplaceReadGroups -I ${SAMPLE_NAME}.${RUN_ID}.sorted.bam -O ${SAMPLE_NAME}.${RUN_ID}.sorted.RG.bam \
                                --RGID ${RUN_ID} --RGLB ${LIBRARY_NAME} --RGPL ${INSTRUMENT_PLATFORM} \
                                --RGPU ${RUN_NAME} --RGSM ${SAMPLE_NAME} --RGPI ${INSERT_SIZE}

    samtools view -H ${SAMPLE_NAME}.${RUN_ID}.sorted.RG.bam | grep '@RG'

    # Index the alignment
    # Command: samtools index
    # Input: alignment (.bam)
    # Ouput: indexed alignment (.bam.bai)
    samtools index ${SAMPLE_NAME}.${RUN_ID}.sorted.RG.bam

    # Append the file name (.bam) to the list of alignments that will be merged
    echo ${SAMPLE_NAME}.${RUN_ID}.sorted.RG.bam >> ${SAMPLE_NAME}.bamlist
done

# Merge the list of alignments into a single file
# Command: samtools merge
# Input: file containing the list of alignments (each line is a .bam file)
# Ouput: alignment (.bam)
samtools merge -b ${SAMPLE_NAME}.bamlist father.bam

# Index the alignment
# Command: samtools index
# Input: alignment (.sam or .bam)
# Ouput: indexed alignment (.sam.bai or .bam.bai)
samtools index father.bam

samtools view -H father.bam | grep '@RG'

#
######################
### Variant Calling ##
######################
#
## GATK best practices pipeline:
## 1. Mark duplicates
## 2. Realigned around known indels
## 3. Recalibrate errors using known SNP
#
## The Ftp folder for the known indels and snps
#ftpMapping=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/
## Use of the Mills and 1000G gold standard dataset for known indels
#known_indels=Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz
## Use of the dbSNP v142 dataset
#known_snps=ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz
#
## Downloading the known indels and SNPs if necessary
#wget ${ftpMapping}/${known_indels} -O ./${known_indels}
#wget ${ftpMapping}/${known_indels}.tbi -O ./${known_indels}.tbi
#wget ${ftpMapping}/${known_snps} -O ./${known_snps}
#wget ${ftpMapping}/${known_snps}.tbi -O ./${known_snps}.tbi
#
#gatk MarkDuplicates \
#    -I father.bam \
#    -O father.markDup.bam \
#    --METRICS_FILE father.markDup_metrics.txt
#
#samtools index father.markDupRG.bam
#
#gatk CreateSequenceDictionary -R Homo_sapiens.Chr20.fa
#
## 2. Realigned around known indels
#gatk  RealignerTargetCreator \
#    -R Homo_sapiens.Chr20.fa \
#    -I father.markDupRG.bam \
#    -o father.intervals \
#    --known ./${known_indels}
#
#gatk IndelRealigner \
#    -R Homo_sapiens.Chr20.fa \
#    -I father.markDupRG.bam \
#    -o father.markDupRG.realigned.bam \
#    --targetIntervals father.intervals \
#    -known ${known_indels}
#
## 3. Recalibrate errors using known SNP
#gatk BaseRecalibrator \
#    -R Homo_sapiens.Chr20.fa \
#    -I father.markDupRG.realigned.bam \
#    -o father.markDupRG.realigned.table \
#    -cov ReadGroupCovariate \
#    -cov QualityScoreCovariate \
#    -cov CycleCovariate -cov ContextCovariate \
#    -knownSites ${known_snps}
#
#gatk PrintReads \
#    -R Homo_sapiens.Chr20.fa \
#    -I father.markDupRG.realigned.bam \
#    -o father.markDupRG.realigned.recal.bam \
#    -BQSR father.markDupRG.realigned.table \
#    --disable_bam_indexing
#
#
#gatk HaplotypeCaller \
#    -R Homo_sapiens.Chr20.fa \
#    -I father.markDupRG.realigned.recal.bam \
#    -o father.gvcf \
#    --genotyping_mode DISCOVERY \
#    -variant_index_type LINEAR \
#    -variant_index_parameter 128000 \
#    --emitRefConfidence GVCF
#
## Converting gVCF to VCF
#gatk GenotypeGVCFs \
#    -R Homo_sapiens.Chr20.fa \
#    --variant father.gvcf \
#    --variant mother.gvcf \
#    --variant daughter.gvcf \
#    -o trio.vcf
#
#
#gatk PhaseByTransmission \
#        -R Homo_sapiens.Chr20.fa \
#        -V trio.vcf \
#        -ped ./20130606_g1k.ped.txt \
#        -o trio.phased.vcf
#
#gatk VariantEval \
#        -R Homo_sapiens.Chr20.fa \
#        -o trio.phased.VE.txt \
#        --eval:set1 trio.vcf \
#        --eval:set2 trio.phased.vcf
#
#gatk GenotypeConcordance \
#        -R Homo_sapiens.Chr20.fa \
#        -eval trio.phased.vcf \
#        -comp trio.vcf \
#        -o trio.GC.txt
