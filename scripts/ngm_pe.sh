#!/bin/bash
INPUT_REF=$1; shift
INPUT_R1=$1; shift
INPUT_R2=$1; shift
OUTPUT_BAM=$1; shift
OUTPUT_BAI=$1; shift
OPTS=$1; shift

set -e; set -o pipefail

# Create scratch directory
if ! test -d "$SCRATCH"; then
	echo "*** Creating scratch directory $SCRATCH"
	mkdir -p "$SCRATCH"
	trap "echo \"*** Removing scratch directory $SCRATCH\"; rm -rf \"$SCRATCH\"" EXIT
fi

echo "*** Running NGM on $INPUT_R1 and $INPUT_R2 against $INPUT_REF"
echo "*** Parameters: $OPTS"
ngm $OPTS -r "$INPUT_REF" -1 "$INPUT_R1" -2 "$INPUT_R2" --bam -o "$SCRATCH/mapped.bam" -t $THREADS

echo "*** Scanning NGM output for invalid reads"
picard ValidateSamFile INPUT="$SCRATCH/mapped.bam" \
	TMP_DIR="$SCRATCH" \
	MAX_OUTPUT=1000000000 IGNORE=RECORD_MISSING_READ_GROUP IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND \
	| sed -n 's/^.*Read name \([^,]*\),.*$/\1/p' \
	| sort \
	| uniq \
	> "$SCRATCH/invalid.txt" \
	|| true

invalid_count=$(wc -l < "$SCRATCH/invalid.txt")
if test $invalid_count -gt 0; then
	echo "*** Filtering the following $invalid_count invalid reads"
	cat "$SCRATCH/invalid.txt"
	picard FilterSamReads INPUT="$SCRATCH/mapped.bam" OUTPUT="$SCRATCH/filtered.bam" SORT_ORDER=queryname \
		READ_LIST_FILE="$SCRATCH/invalid.txt" FILTER=excludeReadList \
		TMP_DIR="$SCRATCH" VALIDATION_STRINGENCY=SILENT WRITE_READS_FILES=false
	rm "$SCRATCH/mapped.bam"
else
	echo "*** No invalid reads, moving mapped.bam to filtered.bam"
	mv "$SCRATCH/mapped.bam" "$SCRATCH/filtered.bam"
fi

echo "*** Sorting by genomic coordinate"
picard SortSam INPUT="$SCRATCH/filtered.bam" OUTPUT="$OUTPUT_BAM" SORT_ORDER=coordinate CREATE_INDEX=true \
	TMP_DIR="$SCRATCH" VALIDATION_STRINGENCY=SILENT
rm "$SCRATCH/filtered.bam"
if ! test -f "$OUTPUT_BAI"; then
	echo "Error: $OUTPUT_BAI was not created successfully"
	exit 1
fi
