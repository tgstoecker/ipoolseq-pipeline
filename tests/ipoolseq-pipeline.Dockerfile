FROM ipoolseq-environment
VOLUME /ipoolseq-pipeline/data
VOLUME /ipoolseq-pipeline/.snakemake/log
COPY entrypoint_root.sh entrypoint_dockeruser.sh /
ENTRYPOINT ["/entrypoint_root.sh"]
CMD ["--help"]
COPY ipoolseq-pipeline.ctx/ /ipoolseq-pipeline/
