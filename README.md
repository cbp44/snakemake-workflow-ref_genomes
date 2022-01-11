# ngs-workflow-ref_genomes

Snakemake workflow for downloading Ensembl reference genomes for use in downstream bioinformatics analyses.

## Quickstart

### Vagrant

```shell
# Start and provision the vagrant box
host> vagrant up

# Connect to the vagrant box through SSH
host> vagrant ssh

# Run copier and answer the questions
vagrant> conda activate copier
vagrant> copier /vagrant ~/test
vagrant> cd ~/test
vagrant> conda activate snakemake
vagrant> snakemake --profile local
```

#### Stop VM

To stop the VM use the following command. You will later be able to start the same VM again with `vagrant up`.

```shell
vagrant halt
```

#### Destroy VM

This will permanently halt and delete the VM.

```shell
vagrant destroy
```

#### Firewall issues

For Vagrant to work, SSH and NFS connections need to be allowed. Here are example `ufw` rules to allow them on the bridge interface. You only need this if you are running a restrictive firewall that denys outbound traffic by default.

```shell
host> ufw allow in on virbr1 from 192.168.121.0/24 to any port 2049 proto tcp comment 'NFS'
host> ufw allow out on virbr1 to any port 22 proto tcp comment 'SSH'
```
