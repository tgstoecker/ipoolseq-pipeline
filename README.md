# iPool-Seq for Transposons  
[![Snakemake](https://img.shields.io/badge/snakemake->=7.0.0-brightgreen.svg)](https://snakemake.readthedocs.io)  
  
Development of Transposon based iPool-Seq developed in Cthe Crop Bioinformatics group at the University of Bonn - earlier versions can of be found [here](https://github.com/Cibiv/ipoolseq-pipeline/releases).

## Installing the pipeline

Clone this repo:

```
  git clone https://github.com/tgstoecker/ipoolseq-pipeline.git
  cd ipoolseq-pipeline
```

### Option 1 - Installing a Bioconda environment containing all necessary dependencies

The file "environment.yaml" defines a Conda (https://conda.io) environment that
provides all programs necessary for running the iPool-Seq analysis pipeline.  
Installation of the environment with either conda or mamba & followed by activation:  

```
  mamba env create -f environment.yaml
  conda activate ipoolseq-pipeline
```

### Option 2 - Installing software during runtime of the pipeline
Currently, we have one environment.yaml file that isn't split into smaller chunks for individual rules as it is not very large in size.
Nevertheless all rules are linked via the "conda:" directive to this environment.yaml file so that quick installation at runtime is possible.  
For this simply add the `--use-conda` flag to your snakemake command, e.g.:  


### Option 3 (**recommended**) - use our Docker container
To use the container simply add BOTH flags `--use-conda --use-singularity` to snakemake command, e.g.:  

This will pull our ipoolseq_cbi_transposon docker container and then create rule specific (although for now all rules share one) conda environments from within the container.
Take note that the installation of singularity is up to you and can sometimes be fiddly.
The fastest way on a system for which you do not have admin/sudo rights is using conda and specifically requesting the rather old version of `singularity==3.6.1`.

## Running the Pipeline:  

If need be, adjust the config file under `cfg/` and the files under the transposon dir `cfg/{transposon}/`.
The standard currently is `cfg/PiggyBac_2022/`.  
Here you will find 4 files:

```
 - cassete.fa
 - essential.tab
 - reference.fa
 - annotation.gff3.gz
```

`reference.fa` & `annotation.gff3.gz` are the ref assembly and corresponding annotation for the fungus being investigated.  
`cassete.fa` defines the left and right border transposon sequence with which we identify the insertion sites.  
`essential.tab` is a one column info file designating the current list of essential genes known for the fungus. With this we add additional information to the final report.  
The number of cores (XX) should be adjusted to personal compute environments capabilities.

As input data you must also provide paired-end sequencing reads, optimally split into 5' and 3' end libraries for the transposon. If both libraries were sequenced together this is however no problem.  
These seq. reads should be deposited inside a `data/{transposon}` directory - in our standard case `cfg/PiggyBac_2022/`.
Read file names HAVE to adhere o the following schema:  
`{Name}-{in/out}+{3/5}p.{1/2}.fq.gz`  
  
A complete set of typical input files would thus look like this:  

```
  #in 3p
  A_vs_B-in+3p.1.fq.gz
  A_vs_B-in+3p.2.fq.gz
  #in 5p
  A_vs_B-in+5p.1.fq.gz
  A_vs_B-in+5p.2.fq.gz
  # out 3p
  A_vs_B-out+3p.1.fq.gz
  A_vs_B-out+3p.2.fq.gz
  # out 5p
  A_vs_B-out+5p.1.fq.gz
  A_vs_B-out+5p.2.fq.gz
```

Note that if you have **sequenced 3' & 5' together** and thus do not possess the fine-grained split as in the example above, you can simply copy all files and name one set with `3p` and the other with `5p`.  
Based on the input file name the workflow is capable of automatically discarding "wrong" fragments respectively.  
  
FInally, to start the anyalsis enter:

```
  snakemake --cores XX data/{transposon}/{Name}.dv.html
# or
  snakemake --cores XX --use-conda data/{transposon}/{Name}.dv.html
# or
  snakemake --cores XX --use-conda --use-singularity data/{transposon}/{Name}.dv.html

# depending on replicates or multiple experiments consider:
  snakemake --cores XX data/{transposon}/{Name1, Name2}.{rep1, rep2}.dv.html
```

## More Information

See http://www.cibiv.at/software/ipoolseq-pipeline, and our
[publication (Uhse *et al.*, 2019)](http://doi.org/10.1002/cppb.20097) in
*Current Protocols in Plant Biology* that describes both the web-lab and the
data-analysis parts of iPool-Seq in detail, and include a step-by-step
description of how to use this pipeline.

## References

Simon Uhse, Florian G. Pflug, Arndt von Haeseler, Armin Djamei (2019). Insertion pool sequencing
for insertional mutant analysis in complex host-microbe interactions. *Current Protocols in
Plant Biology* 4: e20097. DOI: [http://doi.org/10.1002/cppb.20097](http://doi.org/10.1002/cppb.20097)

Simon Uhse, Florian G. Pflug, Stirnberg Alexandra, Ehrlinger Klaus, Arndt von Haeseler,
Armin Djamei (2018). In vivo insertion pool sequencing identifies virulence factors in
a complex fungal–host interaction. *PLoS Biology* 16(4): e2005129. DOI:
[10.1371/journal.pbio.2005129](https://doi.org/10.1371/journal.pbio.2005129)
