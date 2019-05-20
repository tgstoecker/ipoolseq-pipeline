#!/bin/bash
set -e
set -o pipefail

# This scrips creates a docker container that contains the environment
# from the container references in the current environment.rev. We
# exploit that docker caches containers to avoid re-creating the container
# if environment.rev has not changed. The script is meant to be sourced,
# and sets CONTAINER to the hash of the container it created

# Check that environment.rev matches the revision of environment.tar.gz
if [ "$(git log -n 1 --pretty=format:%H -- ../environment.tar.gz)" != "$(cat ../environment.rev)" ]; then
	echo "environment.ref doesn't match revision of environment.tar.gz" >&2
	exit 1
fi

if [ "$(docker inspect --format '{{ index .Config.Labels "environment-rev"}}' ipoolseq-environment)" == "$(cat ../environment.rev)" ]; then
	echo "*** Container ipoolseq-environment is up-to-date"
	exit
fi

# Hard-link environment archive and install script into the build context we'll use
echo "*** Creating container context ipoolseq-environment.ctx"
mkdir -p ipoolseq-environment.ctx
trap "echo \"*** Cleaning up container context\"; rm -r ipoolseq-environment.ctx" EXIT

echo "*** Linking current environment and install script into container context"
ln ../install-environment.sh ipoolseq-environment.ctx/
ln ../environment.tar.gz ipoolseq-environment.ctx/
ln ../environment.rev ipoolseq-environment.ctx/

# Build the container
echo "*** Building container ipoolseq-environment containing the current environment"
docker build --pull=true --tag ipoolseq-environment --build-arg ENVIRONMENT_REV="$(cat ../environment.rev)" --file=ipoolseq-environment.Dockerfile .
