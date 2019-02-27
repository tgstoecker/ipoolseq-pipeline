#!/bin/bash

# Enable conda & activate ipoolseq environment
source /opt/conda/etc/profile.d/conda.sh
conda activate ipoolseq

# Run pipeline
cd /ipoolseq-pipeline
snakemake $* || exit 1
