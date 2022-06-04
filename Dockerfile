FROM docker.io/snakemake/snakemake:v7.8.1 AS base

#################### snakemake ####################
FROM base AS snakemake

# Important to set this for using sort commands
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Add useful aliases for use inside the container
RUN echo "alias ls='/bin/ls -lh --color=auto'" >> /etc/bash.bashrc
#################### /snakemake ####################


#################### workflow ####################
FROM snakemake AS workflow

# Where the workflow code will be stored
ENV WORKFLOW_PATH=/mnt/workflow

# Location of snakemake output cache
ENV SNAKEMAKE_OUTPUT_CACHE=/mnt/snakemake-cache

# Where snakemake conda environments get created
ENV SNAKEMAKE_CONDA_PREFIX=/conda-envs

VOLUME ["${SNAKEMAKE_OUTPUT_CACHE}", "${SNAKEMAKE_CONDA_PREFIX}"]

WORKDIR "${WORKFLOW_PATH}"

COPY ./ "${WORKFLOW_PATH}"

# Define volume for the workflow results
VOLUME ["${WORKFLOW_PATH}/resources"]

# Default command runs 3 rules at a time, downloads at most 2 at the same time
#   and uses conda for software environments.
CMD ["snakemake","--cores=1", "--printshellcmds", "--use-conda"]
#################### /workflow ###################
