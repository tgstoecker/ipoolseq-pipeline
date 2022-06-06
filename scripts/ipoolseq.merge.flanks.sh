#!/bin/bash
# ipoolseq.merge.flanks.sh, Copyright 2020 Florian G. Pflug
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

INPUT_5P=$1; shift
INPUT_3P=$1; shift
OUTPUT=$1; shift

set -e; set -o pipefail

# Create scratch directory
if ! test -d "$SCRATCH"; then
	echo "*** Creating scratch directory $SCRATCH"
	mkdir -p "$SCRATCH"
	trap "echo \"*** Removing scratch directory $SCRATCH\"; rm -rf \"$SCRATCH\"" EXIT
fi

# Copy header fields from 5p input, add read group headers
samtools view -H "$INPUT_5P" > "$SCRATCH/header.txt"
echo -e "@RG\tID:5p\tLB:5p" >> "$SCRATCH/header.txt"
echo -e "@RG\tID:3p\tLB:3p" >> "$SCRATCH/header.txt"

# Merge 5p and 3p libraries, override header since samtools otherwise
# doesn't add any RG header lines, which other tools don't like.
ln -s "$(readlink -e "$INPUT_5P")" "$SCRATCH/5p.bam"
ln -s "$(readlink -e "$INPUT_3P")" "$SCRATCH/3p.bam"
samtools merge -@ $THREADS -r -h "$SCRATCH/header.txt" "$OUTPUT" "$SCRATCH/5p.bam" "$SCRATCH/3p.bam"
