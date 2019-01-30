#!/bin/bash
INPUT_RAW="$1"; shift
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
	echo -en "$(bam_fragments "$INPUT_RAW")\\t"
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
	echo -en "-\\t"
	echo -en "$(awk 'BEGIN{n=0;IFS="\t"} ($1!="unmatched"||$1!="ambiguous") {n+=$2} END{print n}' < "$INPUT_COUNT")\\n"
) | tee "$OUTPUT_STATS"
