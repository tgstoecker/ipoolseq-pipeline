#!/bin/bash
# read_stats.sh, Copyright 2019 Florian G. Pflug
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

INPUT_RAW_R1="$1"; shift
INPUT_RAW_R2="$1"; shift
INPUT_MAP="$1"; shift
INPUT_ASSIGN="$1"; shift
INPUT_COUNT="$1"; shift
OUTPUT_STATS="$1"; shift

function bam_fragments() {
	file="$1"; shift
	# Count the number of unique read names (i.e. pairs)
	count="$(samtools view $* "$file" | cut -f1 | sort --parallel ${THREADS:-1} -S 1G -u | wc -l)"
	echo $count
}

function fastq_fragments() {
	file="$1"; shift
	# Count the number of unique read names (i.e. pairs)
	count="$[$(zcat "$file" | wc -l) / 4]"
	echo $count
}

function assigned_bam_pattern_fragments() {
	file="$1"; shift
	pattern="$1"; shift
	# Count the number of unique read names (i.e. pairs)
	count="$(samtools view $* "$file" | grep -v "$pattern" | cut -f1 | sort --parallel ${THREADS:-1} -S 1G -u | wc -l)"
	echo $count
}

function assigned_bam_pattern_umis() {
	file="$1"; shift
	pattern="$1"; shift
	# Count the number of unique combinations of XT tag and UMI
	count="$(samtools view $* "$file" | grep -v "$pattern" | sed -n 's/^.*_\([ACGTN]*\)\t.*XT:Z:\([^\t].*\).*$/\2-\1/p' | sort -u --parallel ${THREADS:-1} -S 1G -u | wc -l)"
	echo $count
}

(
	echo -en "After Step\\t#Reads\\t#UMIs\\n";

	echo -en "Sequencing\\t"
	echo -en "$(fastq_fragments "$INPUT_RAW_R1")\\t"
	echo -en "-\\n"

	echo -en "Trimming\\t"
	echo -en "$(bam_fragments "$INPUT_MAP")\\t"
	echo -en "-\\n"

	echo -en "Mapping\\t"
	echo -en "$(bam_fragments "$INPUT_MAP" -F4)\\t"
	echo -en "-\\n"

	echo -en "Assignment\\t"
	echo -en "$(assigned_bam_pattern_fragments "$INPUT_ASSIGN" "XT:Z:(unmatched|ambiguous)" -F4)\\t"
	echo -en "$(assigned_bam_pattern_umis "$INPUT_ASSIGN" "XT:Z:(unmatched|ambiguous)" -F4)\\n"

	echo -en "TRUmiCount\\t"
	echo -en "$(awk 'BEGIN{n=0;IFS="\t"} ($1!="unmatched"||$1!="ambiguous") {n+=$8} END{print n}' < "$INPUT_COUNT")\\t"
	echo -en "$(awk 'BEGIN{n=0;IFS="\t"} ($1!="unmatched"||$1!="ambiguous") {n+=$2} END{print n}' < "$INPUT_COUNT")\\n"
) | tee "$OUTPUT_STATS"
