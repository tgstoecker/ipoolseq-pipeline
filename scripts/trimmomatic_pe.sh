#!/bin/bash
# trimmomatic_pe.sh, Copyright 2016, 2017, 2018, 2019 Florian G. Pflug
#
# This file is part of the iPool-Seq Analysis Pipeline
#
# The iPool-Seq Analysis Pipeline is free software: you can redistribute it
# and/or modify it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# The iPool-Seq Analysis Pipeline is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with the iPool-Seq Analysis Pipeline.  If not, see
# <http://www.gnu.org/licenses/

INPUT_R1=$1; shift
INPUT_R2=$1; shift
OUTPUT_R1=$1; shift
OUTPUT_R2=$1; shift
TOM_CMDS=$1; shift

set -e; set -o pipefail
SCRATCH="${SCRATCH:-/tmp}"

# Create scratch directory
if ! test -d "$SCRATCH"; then
	echo "*** Creating scratch directory $SCRATCH"
	mkdir -p "$SCRATCH"
	trap "echo \"*** Removing scratch directory $SCRATCH\"; rm -rf \"$SCRATCH\"" EXIT
fi

echo "*** Running Trimm-O-Matic on $INPUT_R1 and $INPUT_R2"

echo "Options: $TOM_OPTS"
echo "Trim Commands: $TOM_CMDS"
trimmomatic PE -threads ${THREADS:-1} \
	"$INPUT_R1" "$INPUT_R2" \
	"$SCRATCH/tom.1p.fq.gz" "$SCRATCH/tom.1u.fq.gz" \
	"$SCRATCH/tom.2p.fq.gz" "$SCRATCH/tom.2u.fq.gz" \
	$TOM_CMDS

echo "*** Synthesiying dummy mates for unpaired output"

(zcat "$SCRATCH/tom.1p.fq.gz"
 zcat "$SCRATCH/tom.1u.fq.gz"
 zcat "$SCRATCH/tom.2u.fq.gz" | sed -n 's|^@\(.*\)/2$$|@\1/1\nN\n+\n!|p'
) | gzip -c > "$OUTPUT_R1" & READ1_PID=$!

(zcat "$SCRATCH/tom.2p.fq.gz"
 zcat "$SCRATCH/tom.1u.fq.gz" | sed -n 's|^@\(.*\)/1$$|@\1/2\nN\n+\n!|p'
 zcat "$SCRATCH/tom.2u.fq.gz"
) | gzip -c > "$OUTPUT_R2" & READ2_PID=$!

wait "$READ1_PID" || exit 1
wait "$READ2_PID" || exit 1

if [ $(zcat "$OUTPUT_R1" | wc -l) -ne $(zcat "$OUTPUT_R2" | wc -l) ]; then
	echo "Output not properly paired" >&2
	exit 1
fi
