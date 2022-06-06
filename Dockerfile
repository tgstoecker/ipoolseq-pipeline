FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="e4ba24f988922046f31fe1170126f80ee848b8262264ead4d24999e37018d6c1"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: environment.yaml
#   prefix: /conda-envs/0a936de3d3823dd271ab824aef854c3c
#   name: cbi_transposon_ipoolseq
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#     - cibiv
#   dependencies:
#     - git
#     - git-lfs
#     - bcbiogff>=0.6.4
#     - biopython>=1.73
#     - bzip2
#     - conda-ecosystem-user-package-isolation
#     - curl
#     - fastqc>=0.11
#     - gawk>=5
#     - namedlist>=1.7
#     - nextgenmap=0.5.5
#     - picard>=2.20
#     - pysam>=0.15
#     - bioconductor-rtracklayer>=1.42
#     - r-data.table>=1.12
#     - r-dt>=0.6
#     - r-gwpcr>=0.9.10
#     - r-rmarkdown>=1.12
#     - r-plotly>=4.8.0
#     - r-venndiagram>=1.6.20
#     - r-rmdformats>=0.3.5
#     - r-zip
#     - samtools>=1.9
#     - snakemake>=5.4.5
#     - trimmomatic>=0.39
#     - trumicount>=0.9.13
#     - umi_tools>=1.1.2
#     - jupyterlab
#     - r-irkernel
RUN mkdir -p /conda-envs/0a936de3d3823dd271ab824aef854c3c
COPY environment.yaml /conda-envs/0a936de3d3823dd271ab824aef854c3c/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/0a936de3d3823dd271ab824aef854c3c --file /conda-envs/0a936de3d3823dd271ab824aef854c3c/environment.yaml && \
    mamba clean --all -y
