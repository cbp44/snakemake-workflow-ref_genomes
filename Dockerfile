############ base ####################
FROM docker.io/debian:bullseye-slim AS base

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Install ca-certificates and force apt to use https
RUN set -eux; \
	apt-get update; \
	apt-get install -q -y --no-install-recommends \
    ca-certificates; \
  # Replace http with https in apt config
  sed -i 's/^deb http:/deb https:/g' /etc/apt/sources.list; \
  # Clean apt files
  apt-get clean; \
	rm -rf /var/lib/apt/lists/* /usr/share/doc/* /usr/share/man/*

RUN set -eux; \
	apt-get update; \
	apt-get install -q -y --no-install-recommends \
    bzip2 \
    curl \
    git \
    gnupg2 \
    gpg \
    squashfs-tools \
    wget; \
  # Clean apt files
  apt-get clean; \
	rm -rf /var/lib/apt/lists/* /usr/share/doc/* /usr/share/man/*
#################### /base ####################

#################### conda ####################
FROM base AS conda

RUN set -eux; \
  # Download anaconda gpg key
  curl https://repo.anaconda.com/pkgs/misc/gpgkeys/anaconda.asc \
    | gpg --dearmor > /tmp/conda.gpg; \
  # Install the file to the right place
  install -o root -g root -m 644 \
    /tmp/conda.gpg /usr/share/keyrings/conda-archive-keyring.gpg; \ 
  rm -f /tmp/conda.gpg; \
  # Check that fingerprint is correct
  gpg --keyring /usr/share/keyrings/conda-archive-keyring.gpg \
    --no-default-keyring \
    --fingerprint 34161F5BF5EB1D4BFBBB8F0A8AEB4F8B29D82806; \
  echo "deb [arch=amd64 signed-by=/usr/share/keyrings/conda-archive-keyring.gpg] https://repo.anaconda.com/pkgs/misc/debrepo/conda stable main" \
    > /etc/apt/sources.list.d/conda.list; \
  chmod 0644 /etc/apt/sources.list.d/conda.list; \
  apt-get update; \
  apt-get install -q -y --no-install-recommends \
    conda=4.9.2-0; \
  ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh; \
  # Delete library files
  find /opt/conda/ -follow -type f -name '*.h' -or -name '*.a' -or -name '*.js.map' -delete; \
  # Clean apt files
  apt-get clean; \
	rm -rf /var/lib/apt/lists/* /usr/share/doc/* /usr/share/man/* /opt/conda/man/*

COPY ./docker/condarc /opt/conda/.condarc

ENV PATH=/opt/conda/bin:$PATH

# Add shell init code to activate conda and load its base environment
RUN set -eux; \
  conda init bash; \
  echo "conda activate base" >> ~/.bashrc

SHELL ["/bin/bash","-c"]
#################### /conda ####################

#################### mamba ####################
FROM conda AS mamba

RUN set -eux; \
  # Install mamba for faster dependency resolution than conda
  conda install -y -c conda-forge \
    conda-forge::mamba==0.19.0; \
  # Clean up anything left over to trim down the image
  conda clean --all -y
#################### /mamba ####################

#################### snakemake ####################
FROM mamba AS snakemake

# Define SNAKEMAKE_VERSION build arg that persists in Docker environment
ARG SNAKEMAKE_VERSION
ENV SNAKEMAKE_VERSION=${SNAKEMAKE_VERSION:-6.12.1}

# Where snakemake conda environments get created
ENV SNAKEMAKE_CONDA_PREFIX=/mnt/snakemake-envs

# Location of snakemake output cache
ENV SNAKEMAKE_OUTPUT_CACHE=/mnt/snakemake-cache

RUN set -eux; \
  mamba create -q -y \
    -c conda-forge -c bioconda \
    -n snakemake \
    bioconda::snakemake=="${SNAKEMAKE_VERSION}" \
    conda-forge::singularity==3.8.5; \    
  mamba clean --all -y

RUN set -eux; \
  sed -i 's/conda activate base/conda activate snakemake/' ~/.bashrc; \
  mkdir -p \
    "$SNAKEMAKE_CONDA_PREFIX" "SNAKEMAKE_OUTPUT_CACHE"; \
  echo "export SNAKEMAKE_OUTPUT_CACHE=$SNAKEMAKE_OUTPUT_CACHE" >> ~/.bashrc; \
  echo "export SNAKEMAKE_CONDA_PREFIX=$SNAKEMAKE_CONDA_PREFIX" >> ~/.bashrc

ENV PATH=/opt/conda/envs/snakemake/bin:$PATH
#################### /snakemake ####################

#################### envs ####################
FROM snakemake AS envs

WORKDIR /tmp/create_envs

COPY ./docker/create_envs.py /tmp/create_envs/

COPY ./workflow /tmp/create_envs/workflow

COPY ./config /tmp/create_envs/config

RUN set -eu; \
  # Activate snakemake environment
  source ~/.bashrc; \
  conda activate snakemake; \
  # Run script to create environments
  chmod +x /tmp/create_envs/create_envs.py; \
  /tmp/create_envs/create_envs.py; \
  # Have snakemake create environments
  snakemake --cores=all \
    --use-conda --conda-create-envs-only --conda-cleanup-pkgs=cache; \
  # Delete the tmp workspace
  rm -rf /tmp/create_envs
#################### /envs ####################

#################### workflow ####################
FROM envs AS workflow

WORKDIR /usr/src/code

COPY . /usr/src/code

VOLUME ["$SNAKEMAKE_OUTPUT_CACHE", \
        "$SNAKEMAKE_CONDA_PREFIX"]

CMD ["snakemake", "--cores=all", "--use-conda", "--cache"]
#################### workflow ####################
