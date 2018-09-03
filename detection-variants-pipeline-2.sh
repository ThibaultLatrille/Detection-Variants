#!/bin/bash
cd ~/DetectionVariants/

########################################################################################################################
# Requirements:
#	Java (version 8)
java -version
#	FastQC (version 0.11.7)
fastqc -version
#	BWA-MEM (version 0.7.17-r1194-dirty)
bwa
#	SAMtools (version 1.9)
samtools
#	IGV (version 2.4.14)
igv.sh &
#	GATK (version 4.0.8.1)
gatk --list
########################################################################################################################


##########################################################
## Download, extract and index the reference chromosome ##
##########################################################

# Download the reference Human chromosome (chromosome 20) from Ensembl
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: downloaded file
wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz -O Homo_sapiens.Chr20.fa.gz

# Extract the reference chromosome
# Command: gunzip
# Input: compressed file (.gz)
# Ouput: extracted file (remove .gz file)
gunzip Homo_sapiens.Chr20.fa.gz

# Index the reference chromosome
# Command: bwa index
# Input: reference (.fa)
# Ouput: indexed reference (.fa.amb, .fa.ann, .fa.bwt, fa.pac, .fa.sa)
bwa index Homo_sapiens.Chr20.fa

# The sequences are from an East Asian (Kinh Vietnamese) family forming a trio : daughter/mother/father
# Data available at http://www.internationalgenome.org/data-portal/sample/HG02024
# Daughter:
#       StudyId: SRP004063
#       SampleName: HG02024
#       Library: Pond-206419
#       ExperimentID: SRX001595
#       RunId: SRR822251
#       InstrumentModel: Illumina HiSeq 2000
# Mother:
#       StudyId: SRP004063
#       SampleName: HG02025
#       Library: Catch-88584
#       ExperimentID: SRX103760
#       RunId: SRR361100
#       InstrumentModel: Illumina HiSeq 2000
# Father:
#       StudyId: SRP004063
#       SampleName: HG02026
#       Library: Catch-111934
#       ExperimentID: SRX111443
#       RunId: SRR389700
#       InstrumentModel: Illumina HiSeq 2000

# Download paired reads for the daughter (SampleName: HG02024, RunId: SRR822251)
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: downloaded file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02024/sequence_read/SRR822251_1.filt.fastq.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02024/sequence_read/SRR822251_2.filt.fastq.gz

# Extract paired reads for the daughter
# Command: gunzip
# Input: compressed file (.gz)
# Ouput: extracted file (remove .gz file)
gunzip SRR822251_1.filt.fastq.gz
gunzip SRR822251_2.filt.fastq.gz

# Map the paired reads against the reference Human chromosome 20
# Command: bwa mem
# Input: indexed reference (.fa), and sequences (.fastq)
# Ouput: alignment (.sam)
bwa mem -t 4 -M Homo_sapiens.Chr20.fa SRR822251_1.filt.fastq SRR822251_2.filt.fastq > HG02024_SRR822251.sam

# Compute summary statistics of the alignment
# Command: samtools flagstats
# Input: alignment (.sam or .bam)
# Ouput: text file (human and computer readable)
samtools flagstat HG02024_SRR822251.sam > HG02024_SRR822251.sam.flagstats

# Compress the alignment and filter unaligned reads
# Command: samtools view
# Input: alignment (.sam)
# Ouput: compressed alignment (.bam)
samtools view -@ 4 HG02024_SRR822251.sam -bh -f 3 > HG02024_SRR822251.bam

# Sort the alignment
# Command: samtools sort
# Input: alignment (.sam)
# Ouput: compressed alignment (.bam)
samtools sort -@ 4 HG02024_SRR822251.bam > daughter.bam

# Compute statistics of the alignment
# Command: samtools-stats
# Input: alignment (.sam or .bam)
# Ouput: text file (human and computer readable)
samtools stats daughter.bam > daughter.bam.stats

# Plot statistics of the alignment
# Command: plot-bamstats
# Input: statistics text file (output of samtools-stats)
# Ouput: plots (.png)
plot-bamstats -p ./plots-daughter/ daughter.bam.stats

# Index the alignment
# Command: samtools index
# Input: alignment (.sam)
# Ouput: indexed alignment (.bam.bai)
samtools index daughter.bam

# Download for HG02025 Father
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02025/sequence_read/SRR361100_1.filt.fastq.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02025/sequence_read/SRR361100_2.filt.fastq.gz
gunzip SRR361100_1.filt.fastq.gz
gunzip SRR361100_2.filt.fastq.gz
bwa mem -t 4 -M Homo_sapiens.Chr20.fa SRR361100_1.filt.fastq SRR361100_2.filt.fastq | samtools view -@ 4 -bh -f 3 | samtools sort -@ 4 > mother.bam
samtools index mother.bam

# Download for HG02026 Father
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02026/sequence_read/SRR389700_1.filt.fastq.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02026/sequence_read/SRR389700_2.filt.fastq.gz
gunzip SRR389700_1.filt.fastq.gz
gunzip SRR389700_2.filt.fastq.gz
bwa mem -t 4 -M Homo_sapiens.Chr20.fa SRR389700_1.filt.fastq SRR389700_2.filt.fastq | samtools view -@ 4 -bh -f 3 | samtools sort -@ 4 > father.bam
samtools index father.bam

gatk MarkDuplicates -I father.bam -O father.markDup.bam -M METRICS_FILE=father.markDup_metrics.txt

gatk AddOrReplaceReadGroups -I father.bam -O father.bam \
                    RGID=${rgid}\
                    RGLB=lib.${rgid}\
                    RGPU=${rgid}\
                    RGPL=illumina \
                    RGSM=${sampleId}