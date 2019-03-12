#!/bin/bash
set -e
set -o pipefail

DATA_DIR="$(pwd)/data"
mkdir -p "$DATA_DIR"

# Update the container containing the conda enironment if necessary
./update_ipoolseq-environment.sh

# Add the pipeline code to the container created or updated above
mkdir -p ipoolseq-pipeline.ctx
cp ../Snakefile ../VERSION ipoolseq-pipeline.ctx
cp -r ../scripts ipoolseq-pipeline.ctx
cp -r ../cfg ipoolseq-pipeline.ctx
docker build --tag ipoolseq-pipeline --file=ipoolseq-pipeline.Dockerfile .
echo "*** Container ipoolseq-pipeline contains the up-to-date pipeline code"

# Run the pipeline
docker run --rm \
	--volume "$DATA_DIR:/ipoolseq-pipeline/data" \
	--env DOCKERUSER_UID=$(id -u) --env DOCKERUSER_GID=$(id -g) \
	ipoolseq-pipeline $*
