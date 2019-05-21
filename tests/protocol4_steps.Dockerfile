FROM ubuntu:18.04

RUN apt-get update && \
    apt-get install -y curl && \
    apt-get clean

COPY protocol4_steps.sh /

ENTRYPOINT ["/protocol4_steps.sh"]
