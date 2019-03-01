FROM continuumio/miniconda

RUN apt-get install -y unzip && \
    apt-get clean

COPY protocol4_steps.sh /

ENTRYPOINT ["/protocol4_steps.sh"]
