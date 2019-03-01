#!/bin/bash

# Add the pipeline code to the container created or updated above
docker build --file=protocol4_steps.Dockerfile .
PROT4STEPS_CONTAINER=$(docker build -q --file=protocol4_steps.Dockerfile . | tail -n1 | sed -n 's/^Successfully built \([a-f0-9]*\)$/\1/p')
echo "*** Container $PROT4STEPS_CONTAINER contains the up-to-date protocol4 steps"

# Run the pipeline
docker run --rm $PROT4STEPS_CONTAINER
