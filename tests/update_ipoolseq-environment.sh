#!/bin/bash
set -e
set -o pipefail

# This script creates a docker container that contains a conda environment
# that matches the current ipoolseq.yaml. We exploit that docker caches
# containers to avoid re-creating the container if ipoolseq.yaml has not
# changed. The script is meant to be sourced, and sets CONTAINER to the
# hash of the container it created

# Copy environment definition into the build context we'll use
mkdir -p ipoolseq-environment.ctx
cp ../ipoolseq.yaml ipoolseq-environment.ctx

# Build the container
docker build --pull=true --tag ipoolseq-environment --file=ipoolseq-environment.Dockerfile .
echo "*** Container ipoolseq-environment contains an up-to-date ipoolseq conda environment"
