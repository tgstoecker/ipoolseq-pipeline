# Installing the pipeline

Download the [https://github.com/Cibiv/ipoolseq-pipeline/archive/latest-release.zip](latest version)
of the iPool-Seq pipeline, and unzip it. On a linux terminal, this is achieved with

```
VER=latest-release
URL=http://github.com/Cibiv/ipoolseq-pipeline/archive
curl -L -O $URL/$VER.zip
unzip $VER.zip
cd ipoolseq-pipeline-$VER
```

# Installing Dependencies

The file "ipoolseq.yaml" defines a conda (https://conda.io) environment that
provides all programs necessary for running the iPool-Seq analysis pipeline.
After installing conda (see https://conda.io), a conda environment called
ipoolseq can be created with

```
  conda env create --file ipoolseq.yaml
```

Remember that (as all conda environments), this environment must, before it
can be used, be activated for the current terminal session by doing

```
  conda activate ipoolseq
```

To prevent python and R from using possibly incompatible package version
installed in your home directory, we recommend you additionally do

```
  export PYTHONNOUSERSITE=1
  export R_LIBS_USER='-'
```

Note that these commands also only affect the current terminal session!

# Running the Pipeline on the 12 libraries of Uhse et al.[1]

The following command downloads the raw sequencing reads for the 12 (2
experiments, 3 replicates for each, each replicate consists of an input and
an output pool) libraries from Uhse et al.[2], removes read-throughs and
technical sequences from the reads, maps them to the U. maydis genome, and
counts the number of UMIs per insertional knockout. The number of cores (8)
should be adjusted to the number of cores available.

```
  snakemake --cores 8 data/Uhse_et_al.2018/exp{A,B}.r{1,2,3}.dv.html
```

# References

[1] Uhse et al., In vivo insertion pool sequencing identifies virulence factors
    in a complex fungal–host interaction. PLoS Biol 16(4), 2018.
    DOI: [10.1371/journal.pbio.2005129](https://doi.org/10.1371/journal.pbio.2005129)
