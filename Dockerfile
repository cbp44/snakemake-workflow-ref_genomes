FROM docker.io/cbp44/snakemake-minimal:6.10.0 AS base

WORKDIR /usr/src/code

RUN set -eux; \
  chown snakemake:snakemake /usr/src/code

COPY --chown=snakemake:snakemake . /usr/src/code

CMD ["snakemake", "-c1", "--use-conda", "--cache"]
