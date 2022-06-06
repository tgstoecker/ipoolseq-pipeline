# Snakefile, Copyright 2018, 2019, 2020 Florian G. Pflug
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

from snakemake.utils import min_version

min_version("7.0.0")

include: "scripts/snakemake.inc"

configfile: "cfg/config.yaml"

containerized: "docker://tgstoecker/ipoolseq_cbi_transposon:latest"

with open("VERSION") as f:
	VERSION = f.read().strip()

rule help:
	"""
	"""
	run:
		print("This is the iPool-Seq analysis pipeline version %s\n" % VERSION +
		      "Copyright 2017 - 2019 Florian G. Pflug\n"
		      "\n"
		      "Necessary input files\n"
		      "---------------------\n"
		      "\n"
		      "  cfg/your_design/reference.fa                reference genome\n"
		      "  cfg/your_design/cassette.fa                 cassette end sequences\n"
		      "  cfg/your_design/knockouts.gff               knockout cassette locations\n"
		      "\n"
		      "and for each replicate either\n"
		      "\n"
		      "  data/your_design/your_replicate-in.bam      sequenced input pool (both reads)\n"
		      "  data/your_design/your_replicate-out.bam     sequenced output pool (both reads)\n"
		      "\n"
		      "or\n"
		      "\n"
		      "  data/your_design/your_replicate-in.1.fq.gz  sequenced input pool (read 1)\n"
		      "  data/your_design/your_replicate-in.2.fq.gz  sequenced input pool (read 2)\n"
		      "  data/your_design/your_replicate-out.1.fq.gz sequenced output pool (read 1)\n"
		      "  data/your_design/your_replicate-out.2.fq.gz sequenced output pool (read 2)\n"
		      "\n"
		      "Computing abundance tables for your_replicate or your_design\n"
		      "------------------------------------------------------------\n"
		      "\n"
		      "  snakemake data/your_design/your_replicate-in.count.tab\n"
		      "  snakemake data/your_design/your_replicate-out.count.tab\n"
		      "\n"
		      "Running a differential virulence analysis for your_replicate of your_design\n"
		      "---------------------------------------------------------------------------\n"
		      "\n"
		      "  snakemake data/your_design/your_replicate.dv.tab\n"
		      "\n"
		      "This additionally generates the report data/your_design/your_replicate.dv.html\n"
		      "\n"
		      "Additional information\n"
		      "----------------------\n"
		      "\n"
		      "See README.md\n"
		      "\n"
		      "License\n"
		      "-------\n"
		      "This program is distributed in the hope that it will be useful, but WITHOUT ANY\n"
		      "WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n"
		      "PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.\n",
		      file=sys.stderr)

#rule simulate_trumicount_output:
#	"""Simulates TRUmiCount-generated count tables
#	"""
#	input:
#		gff=config_input_file('knockouts', "data/Simulation/sim-in.count.tab")
#	output:
#		pool_in="data/Simulation/sim-in.count.tab",
#		pool_out="data/Simulation/sim-out.count.tab",
#		truth="data/Simulation/sim.truth.tab"
#	script:	"scripts/ipoolseq.simulate.umicounts.R"

#rule download_uhse_et_al:
#	output:	"data/Uhse_et_al.2018/exp{experiment}.r{replicate}-{pool}.bam"
#	params:
#		pool="{pool}",
#		experiment="{experiment}",
#		replicate="{replicate}"
#	wildcard_constraints:
#		experiment="A|B",
#		replicate="1|2|3",
#		pool="in|out"
#	shell:
#		"if   [[ '{params.pool}' == 'in'  && '{params.experiment}' == 'A' ]]; then ERRID=2190337; FILE='r4896/in{params.replicate}';\n"
#		"elif [[ '{params.pool}' == 'in'  && '{params.experiment}' == 'B' ]]; then ERRID=2190343; FILE='r5157/in{params.replicate}';\n"
#		"elif [[ '{params.pool}' == 'out' && '{params.experiment}' == 'A' ]]; then ERRID=2190334; FILE='r4896/egb73r{params.replicate}';\n"
#		"elif [[ '{params.pool}' == 'out' && '{params.experiment}' == 'B' ]]; then ERRID=2190340; FILE='r5157/od3r{params.replicate}';\n"
#		"fi;\n"
#		"URL='ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR219/ERR'$[$ERRID+{params.replicate}-1]/\"$FILE\"'.bam';\n"
#		"echo \"Downloading $URL into {output}\";\n"
#		"curl -o {output:q}\\\n"
#		"  --continue -\\\n"
#		"  --retry 999\\\n"
#		"  --retry-max-time 0\\\n"
#		"  \"$URL\""

#rule bam_to_fqgz_pe:
#	"""Converts a BAM files contained paired-end reads into two (parallel) compressed FASTQ files
#	
#	From data/{lib}.bam, the two output files data/{lib}.1.fq.fz and data/{lib}.2.fq.gz are created
#	"""
#	input:
#		"data/{dir}/{lib}.bam"
#	output:
#		r1="data/{dir}/{lib}.1.fq.gz",
#		r2="data/{dir}/{lib}.2.fq.gz"
#	log:
#		"data/{dir}/{lib}.bam2fqgz.log"
#	shell:
#		"exec > >(tee {log:q}) 2>&1;\n"
#		"set -e; set -o pipefail;\n"
#		"echo \"*** Splitting {input} into {output.r1} and {output.r2}\";\n"
#		"samtools fastq -c1\\\n"
#		"  -1 {output.r1:q}\\\n"
#		"  -2 {output.r2:q}\\\n"
#		"  {input:q};"

rule bam_idx:
	"""Creates an index for a BAM file
	"""
	input:	"data/{dir}/{lib}.bam"
	output:	"data/{dir}/{lib}.bai"
	conda:
		"environment.yaml"
	shell:
		"samtools index\\\n"
		"  {input:q}\\\n"
		"  {output:q}"


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
	threads: 4
	log:
		"data/{dir}/{lib}.tom.log"
	conda:
		"environment.yaml"
	shell:
		"exec > >(tee {log:q}) 2>&1;\n"
		"export SCRATCH={params.scratch:q}; export THREADS={threads};\n"
		"scripts/trimmomatic_pe.sh\\\n"
		"  {input.r1:q}\\\n"
		"  {input.r2:q}\\\n"
		"  {output.r1:q}\\\n"
		"  {output.r2:q}\\\n"
		"  {params.opts:q}"

#ruleorder: adapter_readthrough_trim_pe > bam_to_fqgz_pe

rule ipoolseq_trim_pe:
	"""Trims the iPoolSeq-specific technical sequences (including UMIs), append the UMI to the read name
	
	The read names in the output FASTQ file carry the UMI as a suffic (separated with _), and contains
	only genomic sequences -- the parts overlapping the KO cassette are removed
	"""
	input:
		r1="data/{dir}/{lib}+{flank}.tom.1.fq.gz",
		r2="data/{dir}/{lib}+{flank}.tom.2.fq.gz",
		fa=config_input_file('cassette', "data/{dir}/{lib}.trim.1.fq.gz")
	output:
		r1="data/{dir}/{lib}+{flank}.trim.1.fq.gz",
		r2="data/{dir}/{lib}+{flank}.trim.2.fq.gz"
	params:
		flank="{flank}"
	priority: 1
	threads: 4
	log:
		"data/{dir}/{lib}+{flank}.trim.log"
	conda:
		"environment.yaml"
	shell:
		"exec > >(tee {log:q}) 2>&1;\n"
		"export THREADS={threads};\n"
		"scripts/ipoolseq.trim.py\\\n"
		"  {input.fa:q}\\\n"
		"  {params.flank:q}\\\n"
		"  {input.r1:q}\\\n"
		"  {input.r2:q}\\\n"
		"  {output.r1:q}\\\n"
		"  {output.r2:q}"

rule fastqc_pe:
	input:	"data/{dir}/{lib}.trim.{ri}.fq.gz",
	output: "data/{dir}/{lib}.fastqc.{ri}.html",
	wildcard_constraints:
		ri="1|2"
	priority: 1
	threads: 4
	conda:
		"environment.yaml"
	shell:
		"zcat {input:q}\\\n"
		"| fastqc\\\n"
		"  --threads {threads}\\\n"
		"  --format fastq\\\n"
		"  --outdir \"$(dirname {output:q})\"\\\n"
		"  stdin:\"$(basename {output:q})\";\n"
		"rm {output:q}_fastqc.zip;\n"
		"mv {output:q}_fastqc.html {output:q};\n"

#ruleorder: ipoolseq_trim_pe > bam_to_fqgz_pe

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
	conda:
		"environment.yaml"
	threads: 8
	shell:
		"exec > >(tee {log:q}) 2>&1;\n"
		"export SCRATCH={params.scratch:q}; export THREADS={threads};\n"
		"scripts/ngm_pe.sh\\\n"
		"  {input.ref:q}\\\n"
		"  {input.r1:q}\\\n"
		"  {input.r2:q}\\\n"
		"  {output.bam:q}\\\n"
		"  {output.bai:q}\\\n"
		"  {params.opts:q}"

ruleorder: map_pe > bam_idx

rule ipoolseq_assign_to_isites_pe:
	"""Assign the mapped reads to transposon insertion sites

	The output reads carry an XT tag that contains the chromosome, position and flank of
	the insertion site. Insertion sites are only accepted if both 5' and 3' fragments are
	found
	"""
	input:
		bam_5p="data/{dir}/{lib}+5p.map.bam",
		bam_3p="data/{dir}/{lib}+3p.map.bam"
	output:
		bam_5p="data/{dir}/{lib}+5p.assign.bam",
		bam_3p="data/{dir}/{lib}+3p.assign.bam",
		gff="data/{dir}/{lib}.isites.gff3.gz"
	log:	"data/{dir}/{lib}.assign.log"
	conda:
		"environment.yaml"
	shell:
		"exec > >(tee {log:q}) 2>&1;\n"
		"scripts/ipoolseq.assign.to.isites.py\\\n"
		"  --output-gff {output.gff:q}\\\n"
		"  {input.bam_5p:q}\\\n"
		"  {input.bam_3p:q}\\\n"
		"  {output.bam_5p:q}\\\n"
		"  {output.bam_3p:q}"

rule trumicount_pe:
	"""Computes the number of UMIs per flank (5' and 3') of each of the knockouts

	Uses TRUmiCount to count the raw number of UMIs (which the help of UMI-Tools), and
	to correct for UMIs that were not observed due to having too low sequencing coverage

	The TRUmiCount parameters are set via cfg/config.yaml
	"""
	input:
		bam="data/{dir}/{lib}.assign.bam",
		bai="data/{dir}/{lib}.assign.bai"
	output: counts="data/{dir}/{lib}.count.tab",
		umis="data/{dir}/{lib}.umis.tab",
		readdist="data/{dir}/{lib}.readdist.tab",
		plot="data/{dir}/{lib}.count.pdf"
	log:	"data/{dir}/{lib}.count.log"
	params:
		opts=config_options('trumicount')
	conda:
		"environment.yaml"
	threads: 16
	shell:
		"exec > >(tee {log:q}) 2>&1;\n"
		"trumicount\\\n"
		"  --input-bam {input.bam:q}\\\n"
		"  --group-per gene\\\n"
		"  --include-filter-statistics\\\n"
		"  --output-counts {output.counts:q}\\\n"
		"  --output-umis {output.umis:q}\\\n"
		"  --output-plot {output.plot:q}\\\n"
		"  --output-readdist {output.readdist:q}\\\n"
		"  --paired\\\n"
		"  --umitools-option --per-gene\\\n"
		"  --umitools-option --gene-tag=XT\\\n"
                "  --molecules 1\\\n"
		"  {params.opts}\\\n" #params.opts can contain MULTPLE options, hence don't quote
		"  --cores {threads}"

rule read_stats:
	"""Collects statistics about the number of reads and UMIs remaing after each step
	"""
	input:
		raw_r1="data/{dir}/{lib}.1.fq.gz",
		raw_r2="data/{dir}/{lib}.2.fq.gz",
		map="data/{dir}/{lib}.map.bam",
		assign="data/{dir}/{lib}.assign.bam",
		counts="data/{dir}/{lib}.count.tab"
	output:
		stats="data/{dir}/{lib}.stats.tab"
	conda:
		"environment.yaml"
	threads: 4
	shell:
		"exec > >(tee {log:q}) 2>&1;\n"
		"export THREADS={threads};\n"
		"scripts/read_stats.sh\\\n"
		"  {input.raw_r1:q}\\\n"
		"  {input.raw_r2:q}\\\n"
		"  {input.map:q}\\\n"
		"  {input.assign:q}\\\n"
		"  {input.counts:q}\\\n"
		"  {output.stats:q}"

rule differential_virulence:
	"""Findes KO strains with higher a lower virulence than the wildtype
	
	Produces an HTML report and a table of log2 fold changes of knockout abundances in
	the output pool relative to a set of known-neutral knockouts (marked with 'Neutral'
	in the knockout GFF file), and normlized for differences in the knockout abundance
	in the input pool
	"""
	input:
		ann=config_input_file('reference_annotation', "data/{dir}/{exp}-in+5p.count.tab"),
		essential=config_input_file('essential_genes', "data/{dir}/{exp}-in+5p.count.tab"),
		isites_in="data/{dir}/{exp}-in.isites.gff3.gz",
		isites_out="data/{dir}/{exp}-out.isites.gff3.gz",
		pool_in_5p="data/{dir}/{exp}-in+5p.count.tab",
		pool_in_3p="data/{dir}/{exp}-in+3p.count.tab",
		pool_out_5p="data/{dir}/{exp}-out+5p.count.tab",
		pool_out_3p="data/{dir}/{exp}-out+3p.count.tab",
		stats_in_5p="data/{dir}/{exp}-in+5p.stats.tab",
		stats_in_3p="data/{dir}/{exp}-in+3p.stats.tab",
		stats_out_5p="data/{dir}/{exp}-out+5p.stats.tab",
		stats_out_3p="data/{dir}/{exp}-out+3p.stats.tab",
		fastqc_html_in_5p_r1="data/{dir}/{exp}-in+5p.fastqc.1.html",
		fastqc_html_in_5p_r2="data/{dir}/{exp}-in+5p.fastqc.2.html",
		fastqc_html_in_3p_r1="data/{dir}/{exp}-in+3p.fastqc.1.html",
		fastqc_html_in_3p_r2="data/{dir}/{exp}-in+3p.fastqc.2.html",
		fastqc_html_out_5p_r1="data/{dir}/{exp}-out+5p.fastqc.1.html",
		fastqc_html_out_5p_r2="data/{dir}/{exp}-out+5p.fastqc.2.html",
		fastqc_html_out_3p_r1="data/{dir}/{exp}-out+3p.fastqc.1.html",
		fastqc_html_out_3p_r2="data/{dir}/{exp}-out+3p.fastqc.2.html",
		trumicount_pdf_in_5p="data/{dir}/{exp}-in+5p.count.pdf",
		trumicount_pdf_in_3p="data/{dir}/{exp}-in+3p.count.pdf",
		trumicount_pdf_out_5p="data/{dir}/{exp}-out+5p.count.pdf",
		trumicount_pdf_out_3p="data/{dir}/{exp}-out+3p.count.pdf",
		rmd="scripts/ipoolseq.differential.virulence.Rmd"
	output:
		table="data/{dir}/{exp}.dv.tab",
		html="data/{dir}/{exp}.dv.html",
		zip="data/{dir}/{exp}.dv.zip"
	params:
		version=VERSION,
		dir="{dir}",
		exp="{exp}"
	conda:
		"environment.yaml"
	script:	"scripts/rmarkdown.render.R"
