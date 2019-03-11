FROM continuumio/miniconda

RUN apt-get install -y unzip && \
    apt-get clean

RUN mkdir /ipoolseq-pipeline-dependencies
COPY ipoolseq-environment.ctx/ipoolseq.yaml /ipoolseq-pipeline-dependencies/ipoolseq.yaml

RUN conda env create --file /ipoolseq-pipeline-dependencies/ipoolseq.yaml && \
    conda clean --all --yes
