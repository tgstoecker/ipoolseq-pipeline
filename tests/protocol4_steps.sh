#!/bin/bash
# Download and unpack the pipeline
# 1.
VER=latest-release
URL=http://github.com/Cibiv/ipoolseq-pipeline/archive
curl -L -O $URL/$VER.tar.gz
tar xzf $VER.tar.gz
cd ipoolseq-pipeline-$VER

# Install & activate environment
# 2.
./install-environment.sh || exit 1
source ./environment/bin/activate

# No error now
set -e

# Test
# 3.
snakemake data/Uhse_et_al.2018/expA.r1.dv.tab

# Add reference, knockouts, libraries
# 4.
mkdir -p data/your_design
mkdir -p cfg/your_design
# 5.
cp cfg/Uhse_et_al.2018/reference.fa cfg/your_design/reference.fa
cp cfg/Uhse_et_al.2018/knockouts.gff cfg/your_design/knockouts.gff
cp cfg/Uhse_et_al.2018/cassette.fa cfg/your_design/cassette.fa
# 6.
curl -o data/your_design/your_replicate-in.bam ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR219/ERR2190337/r4896/in1.bam
curl -o data/your_design/your_replicate-out.bam ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR219/ERR2190334/r4896/egb73r1.bam 

# Trim
# 7.
snakemake data/your_design/your_replicate-in.trim.1.fq.gz
snakemake data/your_design/your_replicate-out.trim.1.fq.gz
# 8.
snakemake data/your_design/your_replicate-in.fastqc.1.html
snakemake data/your_design/your_replicate-in.fastqc.2.html
snakemake data/your_design/your_replicate-out.fastqc.1.html
snakemake data/your_design/your_replicate-out.fastqc.2.html

# Map & assign
# 9.
snakemake data/your_design/your_replicate-out.assign.bam

# Abundance
# 10.
snakemake data/your_design/your_replicate-in.count.tab
snakemake data/your_design/your_replicate-out.count.tab

# Differential virulence
# 13.
snakemake data/your_design/your_replicate.dv.tab
