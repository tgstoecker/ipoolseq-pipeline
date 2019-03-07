#!/bin/bash

DATA_DIR="$(pwd)/data"
mkdir -p "$DATA_DIR"

# Update the container containing the conda enironment if necessary
. ./update_ipoolseq-environment.sh

# Add the pipeline code to the container created or updated above
mkdir -p ipoolseq-pipeline.ctx
cp ../Snakefile ipoolseq-pipeline.ctx
cp -r ../scripts ipoolseq-pipeline.ctx
cp -r ../cfg ipoolseq-pipeline.ctx
(
	echo "FROM $CONTAINER"
	echo "VOLUME /ipoolseq-pipeline/data"
	echo "COPY entrypoint_root.sh entrypoint_dockeruser.sh /"
	echo "ENTRYPOINT [\"/entrypoint_root.sh\"]"
	echo "CMD [\"--help\"]"
	echo "COPY ipoolseq-pipeline.ctx/ /ipoolseq-pipeline/"
) > ipoolseq-pipeline.Dockerfile
docker build --file=ipoolseq-pipeline.Dockerfile .
PIPELINE_CONTAINER=$(docker build -q --file=ipoolseq-pipeline.Dockerfile . | tail -n1 | sed -n 's/^Successfully built \([a-f0-9]*\)$/\1/p')
echo "*** Container $PIPELINE_CONTAINER contains the up-to-date pipeline code"

# Run the pipeline
docker run --rm \
	--volume "$DATA_DIR:/ipoolseq-pipeline/data" \
	--env DOCKERUSER_UID=$(id -u) --env DOCKERUSER_GID=$(id -g) \
	$PIPELINE_CONTAINER \
	$*
