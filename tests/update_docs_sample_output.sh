#!/bin/bash
set -e
set -o pipefail

# Replicates to create sample output for
REPLICATES=("expA.r1")

# Create temporary directory
DATA_DIR="$(mktemp -d --suffix='-ipoolseq-pipeline-data')"
trap "echo \"Removing $DATA_DIR\"; rm -rf \"$DATA_DIR\"" EXIT
mkdir -p "$DATA_DIR/.snakemake/log"

# Update the container containing the conda enironment if necessary
./update_ipoolseq-environment.sh

# Update the container containing the pipeline code if necessary
./update_ipoolseq-pipeline.sh

# Run the pipeline
echo "*** Running the pipeline for replicates ${REPLICATES[*]}"
docker run --rm \
	--volume "$DATA_DIR:/ipoolseq-pipeline/data" \
	--volume "$DATA_DIR/.snakemake/log:/ipoolseq-pipeline/.snakemake/log" \
	--env DOCKERUSER_UID=$(id -u) --env DOCKERUSER_GID=$(id -g) \
	ipoolseq-pipeline \
		--cores 4 \
		$(printf 'data/Uhse_et_al.2018/%s.dv.html' "${REPLICATES[@]}")

# Copy output
echo "*** Copying sample output to ../docs/sample_output"
for r in "${REPLICATES[@]}"; do
	for f in .dv.{html,tab} {-in,-out}.count.{tab,pdf} {-in,-out}.fastqc.{1,2}.html; do
		cp "$DATA_DIR/Uhse_et_al.2018/$r$f" ../docs/sample_output
	done
done
