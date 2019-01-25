include: "scripts/snakemake.inc"

configfile: "cfg/config.yaml"

rule simulate_trumicount_output:
	"""Simulates TRUmiCount-generated count tables
	"""
	input:
		gff=config_input_file('knockouts', "data/Simulation/sim-in.count.tab")
	output:
		pool_in="data/Simulation/sim-in.count.tab",
		pool_out="data/Simulation/sim-out.count.tab",
		truth="data/Simulation/sim.truth.tab"
	script:	"scripts/simulate_trumicount_output.R"

rule download_uhse_et_al:
	output:	"data/Uhse_et_al.2018/exp{experiment}-r{replicate}.{pool}.bam"
	params:
		pool="{pool}",
		experiment="{experiment}",
		replicate="{replicate}"
	shell:
		"if   [[ '{params.pool}' == 'input'  && '{params.experiment}' == 'A' ]]; then ERRID=2190337; FILE='r4896/in{params.replicate}'; "
		"elif [[ '{params.pool}' == 'input'  && '{params.experiment}' == 'B' ]]; then ERRID=2190343; FILE='r5157/in{params.replicate}'; "
		"elif [[ '{params.pool}' == 'output' && '{params.experiment}' == 'A' ]]; then ERRID=2190334; FILE='r4896/egb73r{params.replicate}'; "
		"elif [[ '{params.pool}' == 'output' && '{params.experiment}' == 'B' ]]; then ERRID=2190340; FILE='r5157/od3r{params.replicate}'; "
		"fi;"
		"URL='ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR219/ERR'$[$ERRID+{params.replicate}-1]/\"$FILE\"'.bam';"
		"echo \"Downloading $URL into {output}\";"
		"curl -o {output:q} \"$URL\""

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
		opts=config_options('trimmomatic', required=True),
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
                "  --molecules 1\\\n"
		"  {params.opts}" #params.opts can contain MULTPLE options, hence don't quote

rule differential_abundance:
	"""Compares input (pre-infection) and an output (post-infection) samples for a mutant pool

	Produces an HTML report and a table of log2 fold changes of knockout abundances in
	the output pool relative to a set of known-neutral knockouts (marked with 'Neutral'
	in the knockout GFF file), and normlized for differences in the knockout abundance
	in the input pool
	"""
	input:
		gff=config_input_file('knockouts', "data/Simulation/sim-in.count.tab"),
		pool_in="data/{dir}/{exp}-in.count.tab",
		pool_out="data/{dir}/{exp}-out.count.tab",
		trumicount_pdf_in="data/{dir}/{exp}-in.count.pdf",
		trumicount_pdf_out="data/{dir}/{exp}-out.count.pdf"
	output:
		table="data/{dir}/{exp}.da.tab",
		html="data/{dir}/{exp}.da.html"
	log:	"data/{dir}/{exp}.result.log"
	script:	"scripts/generate_differential_abundance_report.R"
