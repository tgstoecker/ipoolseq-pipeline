# Installing the pipeline

Download the [latest version](https://github.com/Cibiv/ipoolseq-pipeline/archive/latest-release.zip)
of the iPool-Seq pipeline, and unzip it. On a linux terminal, this is achieved with

```
VER=latest-release
URL=http://github.com/Cibiv/ipoolseq-pipeline/archive
curl -L -O $URL/$VER.tar.gz
tar xzf $VER.tar.gz
cd ipoolseq-pipeline-$VER
```

Instead, the desired [release](https://github.com/Cibiv/ipoolseq-pipeline/releases)
can of course also be downloaded and unpacked manually, or cloned using `git clone`.
In that case, all further commands must be entered in a terminal windows whose current
directory is the root directory of the pipeline.

# Installing a Bioconda environment containing all necessary dependencies

The file "environment.yaml" defines a Conda (https://conda.io) environment that
provides all programs necessary for running the iPool-Seq analysis pipeline. To
ensure reproducibility of that environment even if Conda packages are replaced
and removed, our source code repository also contains "environment.tar.gz", a
conda-pack archive of that environent. To unpack that environment into
"./environment" and make it usable, run

```
./install-environment.sh
```

*Note: The script "install-environment.sh" will download "environment.tar.gz" from GitHub
if necessary - as a git LFS object, pristine checkouts or sourcecode archives from
GitHub may contain a pointer file instead of the actual archive. If the download
should fail, you can also manually download
[environment.tar.gz](https://github.com/Cibiv/ipoolseq-pipeline/raw/master/environment.tar.gz)
and place it in the root directory of the pipeline before running "install-environment.sh"*

Remember that (as all conda environments), this environment must, before it
can be used, be activated for the current terminal session by doing

```
source ./environment/bin/activate
```

*This step must be repeated for each new terminal session*

# Running the Pipeline on the 12 libraries of Uhse et al.

The following command downloads the raw sequencing reads for the 12 (2
experiments, 3 replicates for each, each replicate consists of an input and
an output pool) libraries from Uhse et al., removes read-throughs and
technical sequences from the reads, maps them to the U. maydis genome, and
counts the number of UMIs per insertional knockout. The number of cores (8)
should be adjusted to the number of cores available.

```
  snakemake --cores 8 data/Uhse_et_al.2018/exp{A,B}.r{1,2,3}.dv.html
```

# More Information

See http://www.cibiv.at/software/ipoolseq-pipeline, and our
[publication (Uhse *et al.*, 2019)](http://doi.org/10.1002/cppb.20097) in
*Current Protocols in Plant Biology* that describes both the web-lab and the
data-analysis parts of iPool-Seq in detail, and include a step-by-step
description of how to use this pipeline.

# References

Simon Uhse, Florian G. Pflug, Arndt von Haeseler, Armin Djamei (2019). Insertion pool sequencing
for insertional mutant analysis in complex host-microbe interactions. *Current Protocols in
Plant Biology* 4: e20097. DOI: [http://doi.org/10.1002/cppb.20097](http://doi.org/10.1002/cppb.20097)

Simon Uhse, Florian G. Pflug, Stirnberg Alexandra, Ehrlinger Klaus, Arndt von Haeseler,
Armin Djamei (2018). In vivo insertion pool sequencing identifies virulence factors in
a complex fungal–host interaction. *PLoS Biology* 16(4): e2005129. DOI:
[10.1371/journal.pbio.2005129](https://doi.org/10.1371/journal.pbio.2005129)
