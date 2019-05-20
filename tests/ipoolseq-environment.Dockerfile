FROM ubuntu:18.04
ARG ENVIRONMENT_REV
LABEL environment-rev=$ENVIRONMENT_REV

RUN mkdir /ipoolseq-pipeline
COPY ipoolseq-environment.ctx/environment.rev /ipoolseq-pipeline/
COPY ipoolseq-environment.ctx/environment.tar.gz /ipoolseq-pipeline/
COPY ipoolseq-environment.ctx/install-environment.sh /ipoolseq-pipeline/

RUN cd /ipoolseq-pipeline && ./install-environment.sh && rm environment.tar.gz install-environment.sh
