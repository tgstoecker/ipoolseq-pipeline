#!/bin/bash
set -e
set -o pipefail

# Remove all users & groups, add dockeruser and dockergroup with the invoking user's gid
export DOCKERUSER_UID=${DOCKERUSER_UID:-100} DOCKERUSER_GID=${DOCKERUSER_GID:-100}
head -n1 /etc/passwd > /etc/passwd.new;
mv /etc/passwd.new /etc/passwd
head -n1 /etc/group > /etc/group.new;
mv /etc/group.new /etc/group
groupadd -g $DOCKERUSER_GID dockergroup
useradd --shell /bin/bash -u $DOCKERUSER_UID -g $DOCKERUSER_GID -M -d /tmp dockeruser

# Run /entrypoint_dockeruser.sh as dockeruser
mkdir -p /ipoolseq-pipeline/.snakemake
chown dockeruser.dockergroup /ipoolseq-pipeline/.snakemake
exec runuser \
	-s /bin/bash -l \
	-c "exec /entrypoint_dockeruser.sh \$*" \
	dockeruser -- /bin/bash $*
