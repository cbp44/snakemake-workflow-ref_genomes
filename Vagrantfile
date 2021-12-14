# -*- mode: ruby -*-
# vi: set ft=ruby :

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.
Vagrant.configure("2") do |config|
    # For a complete reference, please see the online documentation at
    # https://docs.vagrantup.com.

    # Create a forwarded port mapping which allows access to a specific port
    # within the machine from a port on the host machine. In the example below,
    # accessing "localhost:8080" will access port 80 on the guest machine.
    # NOTE: This will enable public access to the opened port
    # config.vm.network "forwarded_port", guest: 80, host: 8080

    # Create a forwarded port mapping which allows access to a specific port
    # within the machine from a port on the host machine and only allow access
    # via 127.0.0.1 to disable public access
    # config.vm.network "forwarded_port", guest: 80, host: 8080, host_ip: "127.0.0.1"


    config.vm.define :devbox do |devbox|
        devbox.vm.box = "debian/bullseye64"
        devbox.vm.box_version = "11.20211018.1"
    
        # Disable automatic box update checking. If you disable this, then
        # boxes will only be checked for updates when the user runs
        # `vagrant box outdated`. This is not recommended.
        devbox.vm.box_check_update = false
        
        devbox.vm.synced_folder "./", "/vagrant"
        # devbox.vm.synced_folder "./", "/vagrant", type: '9p', disabled: false, accessmode: "mapped", owner: "vagrant", mount: true

        devbox.vm.network :private_network,
            :type => "dhcp",
            :libvirt__network_address => "10.44.44.0",
            :libvirt__forward_mode => "nat",
            :libvirt__guest_ipv6 => "no"

        devbox.vm.provider :libvirt do |domain|
            domain.title = "snakemake-workflow-ref_genomes"

            domain.memory = 4096
            domain.cpus = 2
            domain.cpu_mode = "host-model"
            domain.random_hostname = true
            # domain.cpu_mode = "host-passthrough"
            # domain.nested = true

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
            apt-get install -y ca-certificates curl

            # Install docker-compose
            curl -L "https://github.com/docker/compose/releases/download/1.29.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
            chmod +x /usr/local/bin/docker-compose

            apt-get clean all
        SHELL

        devbox.vm.provision :docker do |docker|
            docker.pull_images "docker.io/continuumio/miniconda3:4.10.3"
            # docker.build_image "/vagrant"
        end 
    end
end
