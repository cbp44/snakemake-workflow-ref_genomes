# ngs-workflow-ref_genomes

Snakemake workflow for downloading Ensembl reference genomes for use in downstream bioinformatics analyses.

## Quickstart

### Docker

```shell
docker build -t ref_genomes .
docker volume create snakemake-cache
docker run -v ${PWD}/workflow:/usr/src/code/workflow \
  -v ${PWD}/config:/usr/src/code/config \
  -v snakemake-cache:/mnt/snakemake-cache \
  ref_genomes
```

### Docker Compose

```shell
docker volume create snakemake-cache
docker-compose build
docker-compose up
```

### Vagrant

```shell
# Allow connections for SSH and NFS
# Only needed if running ufw and you have outbound traffic blocked by default
host> ufw allow in on virbr1 from 192.168.121.0/24 to any port 2049 proto tcp comment 'NFS'
host> ufw allow out on virbr1 to any port 22 proto tcp comment 'SSH'

# Start and provision the vagrant box
host> vagrant up
host> vagrant provision

# Connect to the vagrant box through SSH
host> vagrant ssh

vagrant> cd /vagrant
vagrant> sudo docker volume create snakemake-cache
vagrant> sudo docker-compose build
vagrant> sudo docker-compose run --rm workflow bash

docker> copier /usr/src/code /usr/src/test
docker> cd /usr/src/test
docker> snakemake --use-conda -p --cores=all --cache
```
