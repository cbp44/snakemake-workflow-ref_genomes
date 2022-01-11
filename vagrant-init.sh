#!/bin/bash
set -euxo pipefail


# Configuration
CONDA_VERSION=py39_4.10.3
COOKIECUTTER_VERSION=1.7.3
COPIER_VERSION=5.1.0
JSONSCHEMA_VERSION=4.3.3
MAMBA_VERSION=0.19.1
PYDANTIC_VERSION=1.9.0
SINGULARITY_VERSION=3.8.5
SNAKEDEPLOY_VERSION=0.3.0
SNAKEFMT_VERSION=0.4.4
SNAKEMAKE_VERSION=6.12.3
SNAKEMAKE_CONDA_PREFIX=/mnt/snakemake-envs
SNAKEMAKE_OUTPUT_CACHE=/mnt/snakemake-cache
MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh"
MINICONDA_SHA256SUM="1ea2f885b4dbc3098662845560bc64271eb17085387a70c2ba3f29fff6f8d52f"


# Setup some bash aliases
setup_aliases() {
    cat << EOF | tee /root/.bashrc /home/vagrant/.bash_aliases
alias ls='/bin/ls --color=auto -lha'
EOF
}


# Change to https for apt requests
enable_apt_https() {
    sed -i 's/^deb http:/deb https:/g' /etc/apt/sources.list
}


# Install system dependencies and development tools
install_system_deps() {
    apt-get --no-install-recommends install -qq -y \
        curl \
        git
}


# Install miniconda3
install_miniconda() {
    # Install mamba for faster dependency resolution
    install_mamba() {
        conda install -y -c conda-forge \
            -n base \
            conda-forge::mamba=="${MAMBA_VERSION}"
    }

    local prevdir=$(pwd)
    local workdir=$(mktemp -d)

    cd "${workdir}"
    
    # Download installer
    wget --quiet "${MINICONDA_URL}" -O miniconda.sh

    # Check hash
    # NOTE: sha256 hash at top not 100% verified
    echo "${MINICONDA_SHA256SUM} miniconda.sh" > conda_shasum
    sha256sum --check conda_shasum

    # Run installer
    sh miniconda.sh -b -p /opt/conda
    
    # Setup conda auto-activation
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh
    
    # Setup for root and vagrant users
    /opt/conda/bin/conda init bash
    su -c "/opt/conda/bin/conda init bash" vagrant

    # Cleanup files
    cd "${prevdir}"
    rm -rf "$workdir"

    # Load conda for rest of Vagrantfile
    source /root/.bashrc
    install_mamba
}


# Install Snakemake workflow engine along with useful dev tools
install_snakemake() {
    mamba create -y -n snakemake \
        -c conda-forge -c bioconda \
        conda-forge::cookiecutter=="${COOKIECUTTER_VERSION}" \
        conda-forge::jsonschema=="${JSONSCHEMA_VERSION}" \
        conda-forge::pydantic=="${PYDANTIC_VERSION}" \
        conda-forge::singularity=="${SINGULARITY_VERSION}" \
        bioconda::snakedeploy="${SNAKEDEPLOY_VERSION}" \
        bioconda::snakefmt=="${SNAKEFMT_VERSION}" \
        bioconda::snakemake=="${SNAKEMAKE_VERSION}"

    # Setup paths for Snakemake cache and conda env storage
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
use-conda: true
EOF

    # Activate snakemake env & auto-completion by default for vagrant user
    cat << EOF >> /home/vagrant/.bashrc
conda activate snakemake
complete -o bashdefault -C snakemake-bash-completion snakemake

cd /vagrant
EOF

    # Add an alias to make it easier to run snakemake
    cat << EOF >> /home/vagrant/.bash_aliases
alias sml='snakemake --profile local'
EOF
}


# Install copier template engine
install_copier() {
    apt-get install -y --no-install-recommends --no-install-suggests libc6-dev gcc
    mamba create -y -n copier pip
    conda activate copier
    pip install --no-cache-dir copier=="${COPIER_VERSION}"
    conda activate base
    apt-get purge -y libc6-dev gcc
    apt autoremove -y
}


# Clean conda and apt-get cache to save space
clean_cache() {
    conda clean -afy
    apt-get clean all
}


# Run required setup functions
setup_aliases
enable_apt_https
install_system_deps
install_miniconda
install_snakemake
install_copier
clean_cache
