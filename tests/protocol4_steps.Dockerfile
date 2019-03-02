FROM continuumio/miniconda

VOLUME  /conda-pkgs

RUN echo "pkgs_dirs:\n  - /conda-pkgs" > /opt/conda/.condarc

RUN apt-get install -y unzip && \
    apt-get clean

COPY protocol4_steps.sh /

ENTRYPOINT ["/protocol4_steps.sh"]
