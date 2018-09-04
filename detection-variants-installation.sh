#!/usr/bin/env bash
installationDirectory=~/Shared/Detection-Variants/tools

########################################################################################################################
# JAVA
#   Version: 8
#   Licence: MIT
#   Author: Oracle Corporation
#   URL: https://www.java.com
########################################################################################################################
sudo apt-get install openjdk-8-jdk
sudo apt-get install openjdk-8-jre


########################################################################################################################
# FastQC
#   Version: 0.11.7
#   Licence: BSD, MIT
#   Author: Simon Andrews
#   URL: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#   Citation: Andrews, S. FastQC: a quality control tool for high throughput sequence data. (2010).
########################################################################################################################

# Download and extract
cd ${installationDirectory}
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip

# Export the path of the executable (such that fastqc can be lunched from anywhere)
cd FastQC
chmod 755 fastqc
echo 'export PATH=${PATH}:'$(pwd) >> ~/.profile
source ~/.profile


########################################################################################################################
# BWA-MEM
#   Version: 0.7.17-r1194-dirty
#   Licence: GPLv3
#   Author: Heng Li (lh3@me.com)
#   Repository: https://github.com/lh3/bwa
#   Citation: Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 (2013)
########################################################################################################################

# Download (clone git repository)
cd ${installationDirectory}
git clone https://github.com/lh3/bwa.git

# Create the executable
cd bwa
make

# Export the path of the executable (such that bwa can be lunched from anywhere)
echo 'export PATH=${PATH}:'$(pwd) >> ~/.profile
source ~/.profile


########################################################################################################################
# SAMtools
#   Version: 1.9
#   Licence: BSD, MIT
#   Author: Heng Li
#   URL: http://www.htslib.org/
#   Repository: https://github.com/samtools/samtools
#   Citation: Li H., et al. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9 (2009).
########################################################################################################################

# Download and extract
cd ${installationDirectory}
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar xvjf samtools-1.9.tar.bz2

# Create the executable
cd samtools-1.9
./configure --prefix=$(pwd)
make
make install

# Export the path of the executable (such that samtools can be lunched from anywhere)
echo 'export PATH=${PATH}:'$(pwd)'/bin' >> ~/.profile
source ~/.profile


########################################################################################################################
# Integrative Genomics Viewer (IGV)
#   Version: 2.4.14
#   Licence: MIT
#   Author: James T. Robinson
#   URL: http://software.broadinstitute.org/software/igv/
#   Repository: https://github.com/igvteam/igv/
#   Citation: James T. Robinson, et al. Integrative Genomics Viewer. Nature Biotechnology 29, 24–26 (2011).
########################################################################################################################

# Download and extract
cd ${installationDirectory}
wget http://data.broadinstitute.org/igv/projects/downloads/2.4/IGV_2.4.14.zip
unzip IGV_2.4.14.zip

# Export the path of the executable (such that igv.sh can be lunched from anywhere)
cd IGV_2.4.14
chmod 755 igv.sh
echo 'export PATH=${PATH}:'$(pwd) >> ~/.profile
source ~/.profile


########################################################################################################################
# Genome Analysis ToolKit (GATK)
#   Version: 4.0.8.1
#   Licence: BSD 3-Clause (https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT)
#   Author: Broad Institute, Inc (https://github.com/broadinstitute/gatk/blob/master/AUTHORS.TXT)
#   URL: https://software.broadinstitute.org/gatk/
#   Repository: https://github.com/broadinstitute/gatk
#   Citation: McKenna, Aaron, et al.  The Genome Analysis Toolkit: a MapReduce framework for analyzing next-
#             generation DNA sequencing data. Genome research (2010).
########################################################################################################################

# Download and extract
cd ${installationDirectory}
wget https://github.com/broadinstitute/gatk/releases/download/4.0.8.1/gatk-4.0.8.1.zip
unzip gatk-4.0.8.1.zip

# Export the path of the executable (such that gatk can be lunched from anywhere)
cd gatk-4.0.8.1
chmod 755 gatk
echo 'export PATH=${PATH}:'$(pwd) >> ~/.profile
source ~/.profile


########################################################################################################################
# Picard Suite
#	Version: 2.0
########################################################################################################################

# Download and install
cd ${installationDirectory} 
git clone https://github.com/broadinstitute/picard.git
cd picard/
./gradlew shadowJar

# Export name of the directory
echo 'export PICARDDIR=/pandata/gautier/TP_NGS_THIB/tools/picard/build/libs' >> ~/.profile



