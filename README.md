# iPool-Seq for Transposons
Development of Transposon based iPool-Seq - earlier versions can of be found [here](https://github.com/Cibiv/ipoolseq-pipeline/releases).

## Installing the pipeline

Clone this repo:

```
  git clone https://github.com/tgstoecker/ipoolseq-pipeline.git
  cd ipoolseq-pipeline
```

## Installing a Bioconda environment containing all necessary dependencies

The file "environment.yaml" defines a Conda (https://conda.io) environment that
provides all programs necessary for running the iPool-Seq analysis pipeline.  
Installation of the environment with either conda or mamba & followed by activation:  

```
  mamba env create -f environment.yaml
  conda activate ipoolseq-pipeline
```

## Running the Pipeline on the 12 libraries of Uhse et al.

The following command downloads the raw sequencing reads for the 12 (2
experiments, 3 replicates for each, each replicate consists of an input and
an output pool) libraries from Uhse et al., removes read-throughs and
technical sequences from the reads, maps them to the U. maydis genome, and
counts the number of UMIs per insertional knockout. The number of cores (8)
should be adjusted to the number of cores available.

```
  snakemake --cores XX data/{dir}/exp_name.dv.html
# depending on replicates or multiple experiments consider:
  snakemake --cores XX data/{dir}/{epx}.{rep}.dv.html
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
