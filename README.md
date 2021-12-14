# ngs-workflow-reference_genomes

Reproducible workflow for preparing reference genomes for use in bioinformatics analyses using Snakemake.

## Quickstart

```shell
docker-compose build
docker-compose up
```

## Vagrant

```shell
host> ufw allow in on virbr1 from 192.168.121.0/24 to any port 2049 proto tcp comment 'nfs'
host> ufw allow out on virbr1 to any port 22 proto tcp comment 'ssh'
host> vagrant up
host> vagrant provision
host> vagrant ssh

vagrant> cd /vagrant
vagrant> sudo docker-compose build
vagrant> sudo docker-compose run --rm workflow bash

docker> copier /usr/src/code /usr/src/test
docker> cd /usr/src/test
docker> snakemake --use-conda -p --cores=all --cache
```
