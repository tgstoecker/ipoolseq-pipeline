FROM ubuntu:18.04

COPY protocol4_steps.sh /

ENTRYPOINT ["/protocol4_steps.sh"]
