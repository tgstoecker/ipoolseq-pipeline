include: "scripts/snakemake.inc"

configfile: "cfg/config.yaml"

rule demultiplex_bam_pe:
	input:
		mpx="data/{sample}/multiplexed.bam",
		libs="data/{sample}/libraries.tab"
	output: dynamic("data/{sample}/{lib}.demux.bam")
	params:
		odir="data/{sample}",
		opts=config_options('deML'),
		scratch=default_scratch
	log:
		"data/{sample}/demultiplex.log"

rule bam_to_fqgz_pe:
	"""Converts a BAM files contained paired-end reads into two (parallel) compressed FASTQ files
	
	From data/{lib}.bam, the two output files data/{lib}.1.fq.fz and data/{lib}.2.fq.gz are created
	"""
	input:
		"data/{dir}/{lib}.bam"
	output:
		r1="data/{dir}/{lib}.1.fq.gz",
		r2="data/{dir}/{lib}.2.fq.gz"
	log:
		"data/{dir}/{lib}.bam2fqgz.log"
	shell:
		"exec > >(tee {log:q}) 2>&1;"
		"set -e; set -o pipefail;"
		"echo \"*** Splitting {input} into {output.r1} and {output.r2}\";"
		"samtools fastq -c1 -1 {output.r1:q} -2 {output.r2:q} {input:q};"

rule bam_idx:
	"""Creates an index for a BAM file
	"""
	input:	"data/{dir}/{lib}.bam"
	output:	"data/{dir}/{lib}.bai"
	shell:	"samtools index {input:q} {output:q}"

rule adapter_readthrough_trim_pe:
	"""Trims read-throughs into the adapter on the other end using Trimmomatic
	
	The adapter sequences and the trimming options are set via cfg/config.yaml
	"""
	input:
		r1="data/{dir}/{lib}.1.fq.gz",
		r2="data/{dir}/{lib}.2.fq.gz"
	output:
		r1="data/{dir}/{lib}.tom.1.fq.gz",
		r2="data/{dir}/{lib}.tom.2.fq.gz"
	params:
		opts=config_options('trimmomatic'),
		scratch=default_scratch
	threads: 16
	log:
		"data/{dir}/{lib}.tom.log"
	shell:
		"exec > >(tee {log:q}) 2>&1;"
		"export SCRATCH={params.scratch:q}; export THREADS={threads};"
		"scripts/trimmomatic_pe.sh {input.r1:q} {input.r2:q} {output.r1:q} {output.r2:q} {params.opts:q}"

rule ipoolseq_transposon_trim_pe:
	"""Trims the iPoolSeq-specific technical sequences (including UMIs), append the UMI to the read name
	
	The read names in the output FASTQ file carry the UMI as a suffic (separated with _), and contains
	only genomic sequences -- the parts overlapping the KO cassette are removed
	"""
	input:
		r1="data/{dir}/{lib}.tom.1.fq.gz",
		r2="data/{dir}/{lib}.tom.2.fq.gz"
	output:
		r1="data/{dir}/{lib}.trim.1.fq.gz",
		r2="data/{dir}/{lib}.trim.2.fq.gz"
	threads: 16
	log:
		"data/{dir}/{lib}.trim.log"
	shell:
		"exec > >(tee {log:q}) 2>&1;"
		"export THREADS={threads};"
		"scripts/ipoolseq.transposon.trim.py {input.r1:q} {input.r2:q} {output.r1:q} {output.r2:q}"

rule map_pe:
	"""Maps the (trimmed) reads to the genome

	The reference genome is set via cfg/config.yaml
	"""
	input:
		r1="data/{dir}/{lib}.trim.1.fq.gz",
		r2="data/{dir}/{lib}.trim.2.fq.gz",
		ref=config_input_file('reference_genome', "data/{dir}/{lib}.map.bam")
	output:
		bam="data/{dir}/{lib}.map.bam",
		bai="data/{dir}/{lib}.map.bai"
	log:    "data/{dir}/{lib}.map.log"
	params:
		opts=config_options('ngm'),
		scratch=default_scratch
	threads: 32
	shell:
		"exec > >(tee {log:q}) 2>&1;"
		"export SCRATCH={params.scratch:q}; export THREADS={threads};"
		"scripts/ngm_pe.sh {input.ref:q} {input.r1:q} {input.r2:q} {output.bam:q} {output.bai:q} {params.opts:q}"

rule ipoolseq_assign_to_knockouts_pe:
	"""Assign the mapped reads to the individual KO strains

	The output reads carry an XT tag that states the name of the KO strain (from the GFF files
	showing the positions of the KO cassette insertions) and the flank (3' or 5') of the KO
	cassette that the read belongs to, in the form '<Name>:{3,5}p'

	The GFF file listing the KO cassette insertion positions is set via cfg/config.yaml
	"""
	input:
		bam="data/{dir}/{lib}.map.bam",
		gff=config_input_file('knockouts', "data/{dir}/{lib}.assign.bam")
	output:	"data/{dir}/{lib}.assign.bam"
	log:	"data/{dir}/{lib}.assign.log"
	shell:
		"exec > >(tee {log:q}) 2>&1;"
		"scripts/ipoolseq.assign.to.knockouts.py {input.gff:q} {input.bam:q} {output:q}"

rule trumicount_pe:
	"""Computes the number of UMIs per flank (5' and 3') of each of the knockouts

	Uses TRUmiCount to count the raw number of UMIs (which the help of UMI-Tools), and
	to correct for UMIs that were not observed due to having too low sequencing coverage

	The TRUmiCount parameters are set via cfg/config.yaml
	"""
	input:
		bam="data/{dir}/{lib}.assign.bam",
		bai="data/{dir}/{lib}.assign.bai"
	output: table="data/{dir}/{lib}.count.tab",
		plot="data/{dir}/{lib}.count.pdf"
	log:	"data/{dir}/{lib}.count.log"
	params:
		opts=config_options('trumicount')
	threads: 32
	shell:
		"exec > >(tee {log:q}) 2>&1;\n"
		"trumicount\\\n"
		"  --input-bam {input.bam:q}\\\n"
		"  --group-per gene\\\n"
		"  --cores {threads}\\\n"
		"  --output-counts {output.table:q}\\\n"
		"  --output-plot {output.plot:q}\\\n"
		"  --umitools-option --per-gene\\\n"
		"  --umitools-option --gene-tag=XT\\\n"
		"  {params.opts}" #params.opts can contain MULTPLE options, hence don't quote
