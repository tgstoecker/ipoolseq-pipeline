#!/bin/bash

# Add the pipeline code to the container created or updated above
docker build --tag ipoolseq-protocol4-steps --file=protocol4_steps.Dockerfile .
echo "*** Container ipoolseq-protocol4-steps contains the up-to-date protocol4 steps"

# Run the pipeline
docker run --rm -v conda-pkgs:/conda-pkgs ipoolseq-protocol4-steps
