#!/bin/bash
set -euxo pipefail

# Configuration
CONDA_VERSION=py39_4.10.3
MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh"
SHA256SUM="1ea2f885b4dbc3098662845560bc64271eb17085387a70c2ba3f29fff6f8d52f"
MAMBA_VERSION=0.19.1
SNAKEMAKE_VERSION=6.12.3
SINGULARITY_VERSION=3.8.5
COPIER_VERSION=5.1.0
SNAKEFMT_VERSION=0.4.4
PYDANTIC_VERSION=1.9.0
JSONSCHEMA_VERSION=4.3.3
SNAKEMAKE_CONDA_PREFIX=/mnt/snakemake-envs
SNAKEMAKE_OUTPUT_CACHE=/mnt/snakemake-cache


# Change to https for apt requests
sed -i 's/^deb http:/deb https:/g' /etc/apt/sources.list
apt-get update


# Install dependencies
apt-get --no-install-recommends install -qq -y \
    curl \
    git


## Install miniconda3
workdir=$(mktemp -d)
cd $workdir
wget --quiet "${MINICONDA_URL}" -O miniconda.sh
echo "${SHA256SUM} miniconda.sh" > conda_shasum
mkdir -p /opt
sh miniconda.sh -b -p /opt/conda
sha256sum --check conda_shasum
cd /
rm -rf "$workdir"
ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh
echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/vagrant/.bashrc
echo "conda activate base" >> /home/vagrant/.bashrc
/opt/conda/bin/conda init bash
source /root/.bashrc
conda activate base


## Install mamba for faster dependency resolution
conda install -y -c conda-forge \
    -n base \
    conda-forge::mamba=="${MAMBA_VERSION}"


## Install Snakemake
mamba create -y -c conda-forge -c bioconda \
    -n snakemake \
    bioconda::snakemake-minimal=="${SNAKEMAKE_VERSION}" \
    conda-forge::singularity=="${SINGULARITY_VERSION}"


## Create dev environment with useful tools
mamba create -y -n dev \
    -c conda-forge -c bioconda -c defaults \
    bioconda::snakefmt=="${SNAKEFMT_VERSION}" \
    conda-forge::pydantic=="${PYDANTIC_VERSION}" \
    conda-forge::jsonschema=="${JSONSCHEMA_VERSION}"


## Install copier to copier conda env
apt-get install -y --no-install-recommends --no-install-suggests libc6-dev gcc
mamba create -y -n copier pip
conda activate copier
pip install --no-cache-dir copier=="${COPIER_VERSION}"
conda activate base
apt-get purge -y libc6-dev gcc
apt autoremove -y


# Clean conda and apt-get cache for space
conda clean -afy
apt-get clean all


# Setup snakemake conda env and output cache
mkdir -p "${SNAKEMAKE_CONDA_PREFIX}" "${SNAKEMAKE_OUTPUT_CACHE}"
chown vagrant:vagrant "${SNAKEMAKE_CONDA_PREFIX}" "${SNAKEMAKE_OUTPUT_CACHE}"
echo "export SNAKEMAKE_OUTPUT_CACHE=${SNAKEMAKE_OUTPUT_CACHE}" >> /home/vagrant/.bashrc
echo "export SNAKEMAKE_CONDA_PREFIX=${SNAKEMAKE_CONDA_PREFIX}" >> /home/vagrant/.bashrc
su -c "conda config --append envs_dirs /opt/conda/envs" vagrant
su -c "conda config --append envs_dirs /mnt/snakemake-envs" vagrant


# Create local Snakemake profile --- 'snakemake --profile local'
mkdir -p /etc/xdg/snakemake/local
cat << EOF > /etc/xdg/snakemake/local/config.yaml
conda-cleanup-pkgs: cache
conda-frontend: mamba
conda-prefix: /mnt/snakemake-envs
cores: all
keep-going: false
printshellcmds: true
rerun-incomplete: true
use-conda: true
EOF


# Setup some aliases
cat << EOF >> /home/vagrant/.bashrc
alias ls='/bin/ls --color=auto -lha'
alias sml='snakemake --profile local'

conda activate snakemake
complete -o bashdefault -C snakemake-bash-completion snakemake
EOF
