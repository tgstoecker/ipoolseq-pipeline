#!/bin/bash
# Enable & update conda
. /opt/conda/etc/profile.d/conda.sh
conda update -n base -c defaults conda
export USER=user

# 1a. Download pipeline
VER=latest-release
URL=http://github.com/Cibiv/ipoolseq-pipeline/archive
curl -L -O $URL/$VER.zip
unzip $VER.zip
cd ipoolseq-pipeline-$VER

# 1b. Create & activate environment
conda env create --file ipoolseq.yaml
conda activate ipoolseq

# No errors now
set -e
set -p pipefail

# 2. Add reference, knockouts, libraries
mkdir -p data/your_design
mkdir cfg/your_design
cp cfg/Uhse_et_al.2018/reference.fa cfg/your_design/reference.fa
cp cfg/Uhse_et_al.2018/knockouts.gff cfg/your_design/knockouts.gff
cp cfg/Uhse_et_al.2018/cassette.fa cfg/your_design/cassette.fa
curl -o data/your_design/your_replicate-in.bam ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR219/ERR2190337/r4896/in1.bam
curl -o data/your_design/your_replicate-out.bam ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR219/ERR2190334/r4896/egb73r1.bam 

# 3. Trim
snakemake data/your_design/your_replicate-in.trim.1.fq.gz
snakemake data/your_design/your_replicate-in.fastqc.1.html
snakemake data/your_design/your_replicate-in.fastqc.2.html

snakemake data/your_design/your_replicate-out.trim.1.fq.gz
snakemake data/your_design/your_replicate-out.fastqc.1.html
snakemake data/your_design/your_replicate-out.fastqc.2.html

# 4. Map & assign
snakemake data/your_design/your_replicate-in.assign.bam
snakemake data/your_design/your_replicate-out.assign.bam

# 5. Abundance
snakemake data/your_design/your_replicate-in.count.tab
snakemake data/your_design/your_replicate-out.count.tab

# 7. Differential virulence
snakemake data/your_design/your_replicate.dv.tab
