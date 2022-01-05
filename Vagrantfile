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

            domain.memory = 6144
            domain.cpus = 2
            # domain.cpu_mode = "host-model"

            domain.graphics_type = "vnc"
            domain.graphics_autoport = true
            domain.graphics_ip = "127.0.0.1"
            domain.video_type = "qxl"

            # domain.graphics_type = "none"
            
            domain.volume_cache = "none"
        end

        devbox.vm.provision :shell, path: "vagrant-init.sh"
    end
end
