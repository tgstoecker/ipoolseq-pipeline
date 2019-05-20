#!/bin/bash
# Run snakemake in /ipoolseq-pipeline
cd /ipoolseq-pipeline
source ./environment/bin/activate
exec snakemake $*
