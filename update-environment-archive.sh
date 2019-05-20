#!/bin/bash

# Was the script sourced?
if ! (return 0 2>/dev/null); then
	echo "$0 must be sourced, not executed!" >&2
	exit 1
fi

# Check if conda is available
if ! declare -f conda > /dev/null; then
	echo "conda was not enabled properly" >&2
	return
fi

# Check if there are uncommitted changed
if [ $(git status --porcelain | wc -l) != 0 ]; then
	echo "Working copy contains uncommitted changes" >&2
	return
fi

# Create temporary directory
ENVDIR="$(mktemp -d --suffix='-ipoolseq-env')"
trap "echo \"Removing $ENVDIR\"; rm -rf \"$ENVDIR\"" RETURN

echo "*** Creating $ENVDIR ..."
conda create -p "$ENVDIR" --no-default-packages --yes || return

echo "*** Activating $ENVDIR ..."
conda activate "$ENVDIR" || return

echo "*** Making $ENVDIR self-contained"
mkdir -p "$ENVDIR"/pkgs || return
conda config --env --add pkgs_dirs "$ENVDIR"/pkgs || return

echo "*** Installing packages into $ENVDIR ..."
conda env update -f environment.yaml -p "$ENVDIR" --prune || return

echo "*** Packing environment ..."
cp environment.yaml "$ENVDIR" || return
conda pack -p "$ENVDIR" -o environment.tar.gz --format tar.gz --compress-level 9 --n-threads 4 || return

echo "*** Deactivating environment, returning to previously activate environment ..."
conda deactivate || return

echo "*** Comitting environment.tar.gz ..."
git add environment.tar.gz || return
git commit -m "Updated environment.tar.gz based on environment.yaml rev. $(git log -n 1 --pretty=format:%h -- environment.yaml)" || return

echo "*** Comitting updated environent.rev ..."
git log -n 1 --pretty=format:%H -- environment.tar.gz > environment.rev || return
git add environment.rev || return
git commit -m "Updated environment.rev to $(git log -n 1 --pretty=format:%h -- environment.tar.gz) after updating environment.tar.gz based on environment.yaml rev. $(git log -n 1 --pretty=format:%h -- environment.yaml)" || return
