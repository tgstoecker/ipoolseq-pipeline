#!/bin/bash
# Activate conda
. /opt/conda/etc/profile.d/conda.sh
conda activate ipoolseq

# Run snakemake in /ipoolseq-pipeline
cd /ipoolseq-pipeline
exec snakemake $*
