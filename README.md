# ngs-workflow-ref_genomes

Snakemake workflow for downloading Ensembl reference genomes for use in downstream bioinformatics analyses.

## Quickstart

### Vagrant

```shell
# Allow connections for SSH and NFS
# Only needed if running ufw and you have outbound traffic blocked by default
host> ufw allow in on virbr1 from 192.168.121.0/24 to any port 2049 proto tcp comment 'NFS'
host> ufw allow out on virbr1 to any port 22 proto tcp comment 'SSH'

# Start and provision the vagrant box
host> vagrant up

# Connect to the vagrant box through SSH
host> vagrant ssh

vagrant> cd /vagrant
vagrant> conda activate copier
# Run copier and answer the questions
vagrant> copier /vagrant ~/test
vagrant> cd ~/test
vagrant> conda activate snakemake
vagrant> snakemake --use-conda -c 1 -p
```

#### Destroy VM

```shell
vagrant destroy
```
