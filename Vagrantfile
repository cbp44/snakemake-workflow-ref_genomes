# -*- mode: ruby -*-
# vi: set ft=ruby :

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.
Vagrant.configure("2") do |config|
    # For a complete reference, please see the online documentation at
    # https://docs.vagrantup.com.

    config.vm.define :devbox do |devbox|
        devbox.vm.box = "debian/bullseye64"
        devbox.vm.box_version = "11.20211018.1"
    
        # Disable automatic box update checking. If you disable this, then
        # boxes will only be checked for updates when the user runs
        # `vagrant box outdated`. This is not recommended.
        devbox.vm.box_check_update = true
        
        devbox.vm.synced_folder "./", "/vagrant"

        devbox.vm.provider :libvirt do |domain|
            domain.title = "snakemake-workflow-ref_genomes"

            domain.memory = 4096
            domain.cpus = 2
            domain.cpu_mode = "host-model"

            domain.graphics_type = "vnc"
            domain.graphics_autoport = true
            domain.graphics_ip = "127.0.0.1"
            domain.video_type = "qxl"

            domain.volume_cache = "none"
        end

        devbox.vm.provision :shell, inline: <<-SHELL
            # Change to https for apt requests
            sed -i 's/^deb http:/deb https:/g' /etc/apt/sources.list
            apt-get update
            

            # Install git
            apt-get install -y git
            

            # Install curl
            apt-get install -y curl


            # Install docker-compose
            curl --silent -L "https://github.com/docker/compose/releases/download/1.29.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
            chmod +x /usr/local/bin/docker-compose

            # Install miniconda3
            CONDA_VERSION=py39_4.10.3
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh"
            SHA256SUM="1ea2f885b4dbc3098662845560bc64271eb17085387a70c2ba3f29fff6f8d52f"

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


            # Install mamba for faster dependency resolution
            MAMBA_VERSION=0.19.1
            conda install -y -c conda-forge conda-forge::mamba=="${MAMBA_VERSION}"


            # Install Snakemake
            SNAKEMAKE_VERSION=6.12.3
            mamba create -y -c conda-forge -c bioconda -n snakemake bioconda::snakemake-minimal=="${SNAKEMAKE_VERSION}" conda-forge::singularity==3.8.5

           
            # Install snakefmt into snakefmt conda env
            mamba create -y -n snakefmt pip
            conda activate snakefmt
            pip install --no-cache-dir snakefmt
            conda activate base


            # Install copier
            mamba create -y -n copier pip
            conda activate copier
            apt-get install -y --no-install-recommends --no-install-suggests libc6-dev gcc
            pip install --no-cache-dir copier==5.1.0
            apt-get purge -y libc6-dev gcc
            apt autoremove -y
            conda activate base


            # Clean conda and apt-get cache for space
            conda clean -afy
            apt-get clean all


            # Setup snakemake conda env and output cache
            SNAKEMAKE_CONDA_PREFIX=/mnt/snakemake-envs
            SNAKEMAKE_OUTPUT_CACHE=/mnt/snakemake-cache
            mkdir -p "${SNAKEMAKE_CONDA_PREFIX}" "${SNAKEMAKE_OUTPUT_CACHE}"
            chown vagrant:vagrant "${SNAKEMAKE_CONDA_PREFIX}" "${SNAKEMAKE_OUTPUT_CACHE}"
            echo "export SNAKEMAKE_OUTPUT_CACHE=${SNAKEMAKE_OUTPUT_CACHE}" >> /home/vagrant/.bashrc
            echo "export SNAKEMAKE_CONDA_PREFIX=${SNAKEMAKE_CONDA_PREFIX}" >> /home/vagrant/.bashrc
            su -c "conda config --append envs_dirs /opt/conda/envs" vagrant
            su -c "conda config --append envs_dirs /mnt/snakemake-envs" vagrant
        SHELL

        # devbox.vm.provision :docker do |docker|
        #     # docker.pull_images "docker.io/continuumio/miniconda3:4.10.3"
        # end 
    end
end
