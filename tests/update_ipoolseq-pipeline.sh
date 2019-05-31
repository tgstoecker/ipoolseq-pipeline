#!/bin/bash
set -e
set -o pipefail

# Add the pipeline code to the container created or updated above
echo "*** Building container ipoolseq-pipeline containing the current pipeline code"
mkdir -p ipoolseq-pipeline.ctx
cp ../Snakefile ../VERSION ipoolseq-pipeline.ctx
cp -r ../scripts ipoolseq-pipeline.ctx
cp -r ../cfg ipoolseq-pipeline.ctx
docker build --tag ipoolseq-pipeline --file=ipoolseq-pipeline.Dockerfile .

echo "*** Container ipoolseq-pipeline contains the up-to-date pipeline code"
