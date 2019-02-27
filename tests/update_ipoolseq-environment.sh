#!/bin/bash

# This script creates a docker container that contains a conda environment
# that matches the current ipoolseq.yaml. We exploit that docker caches
# containers to avoid re-creating the container if ipoolseq.yaml has not
# changed. The script is meant to be sourced, and sets CONTAINER to the
# hash of the container it created

# Copy environment definition into the build context we'll use
mkdir -p ipoolseq-environment.ctx
cp ../ipoolseq.yaml ipoolseq-environment.ctx

# Build the container
docker build --pull=true --file=ipoolseq-environment.Dockerfile .

# Re-build with "-q" and extract the hash
# Note that due to docker's cache, this shouldn't actually rebuild
CONTAINER=$(docker build -q --file=ipoolseq-environment.Dockerfile . | tail -n1 | sed -n 's/^Successfully built \([a-f0-9]*\)$/\1/p')

echo "*** Container $CONTAINER contains an up-to-date ipoolseq conda environment"
