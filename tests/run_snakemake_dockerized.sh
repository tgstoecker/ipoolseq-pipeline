#!/bin/bash

# Update the container containing the conda enironment if necessary
. ./update_ipoolseq-environment.sh

# Add the pipeline code to the container created or updated above
mkdir -p ipoolseq-pipeline.ctx
cp ../Snakefile ipoolseq-pipeline.ctx
cp -r ../scripts ipoolseq-pipeline.ctx
cp -r ../cfg ipoolseq-pipeline.ctx
cp run_snakemake.sh ipoolseq-pipeline.ctx
(
	echo "FROM $CONTAINER"
	echo "COPY ipoolseq-pipeline.ctx/ /ipoolseq-pipeline/"
) > ipoolseq-pipeline.Dockerfile
docker build --file=ipoolseq-pipeline.Dockerfile .
PIPELINE_CONTAINER=$(docker build -q --file=ipoolseq-pipeline.Dockerfile . | tail -n1 | sed -n 's/^Successfully built \([a-f0-9]*\)$/\1/p')
echo "*** Container $PIPELINE_CONTAINER contains the up-to-date pipeline code"

# Run the pipeline
docker run -it --rm $PIPELINE_CONTAINER /ipoolseq-pipeline/run_snakemake.sh $*
