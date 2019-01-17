Installing Dependencies
-----------------------

The file "conda.reqs" defines a conda[1] environment that provides all programs
necessary for running the iPool-Seq analysis pipeline. After installing conda
(see https://conda.io), a conda environment can be created with

  conda create -n ipoolseq -c defaults -c conda-forge -c bioconda --file conda.reqs

Remember that (as all conda environments), this environment must than be activ-
ated before it can be used by doing

  conda activate ipoolseq

To prevent python and R from using possibly incompatible package version
installed in your home directory, we recommend you additionally do

  export PYTHONNOUSERSITE=1
  export R_LIBS_USER='-'

Running the Pipeline on the 12 libraries of Uhse et al.[2]
--------------------------------------------------------------

The following command downloads the raw sequencing reads for the 12 (2
experiments, 3 replicates for each, each replicate consists of an input and
an output pool) libraries from Uhse et al.[2], removes read-throughs and
technical sequences from the reads, maps them to the U. maydis genome, and
counts the number of UMIs per insertional knockout. The number of cores (8)
should be adjusted to the number of cores available.

  snakemake --cores 8 data/Uhse_et_al.2018/exp{A,B}-r{1,2,3}.{in,out}.count.tab

More Information
----------------

See Supplementary Information S1 Text "iPool-Seq Analysis Pipeline".

References
----------

[1] https://conda.io/docs/

[2] Uhse et al., In vivo insertion pool sequencing identifies virulence factors
    in a complex fungal–host interaction. PLoS Biol 16(4), 2018.
