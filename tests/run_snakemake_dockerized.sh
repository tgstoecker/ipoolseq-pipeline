#!/bin/bash
set -e
set -o pipefail

DATA_DIR="$(pwd)/data"
mkdir -p "$DATA_DIR"

LOG_DIR="$(pwd)/snakemake-log"
mkdir -p "$LOG_DIR"

# Update the container containing the conda enironment if necessary
./update_ipoolseq-environment.sh

# Update the container containing the pipeline code if necessary
./update_ipoolseq-pipeline.sh

# Run the pipeline
docker run --rm \
	--volume "$DATA_DIR:/ipoolseq-pipeline/data" \
	--volume "$LOG_DIR:/ipoolseq-pipeline/.snakemake/log" \
	--env DOCKERUSER_UID=$(id -u) --env DOCKERUSER_GID=$(id -g) \
	ipoolseq-pipeline $*
