# snakemake-workflow-ref_genomes

Snakemake workflow for downloading Ensembl reference genomes for use in downstream bioinformatics analyses.

## Getting Started

### Docker

```shell
# Create conda-env Docker volume to share conda envs among workflows
docker volume create conda-envs
# docker volume create --opt type=none --opt o=bind --opt device=/media/user/data/conda-envs conda-envs

# Create named bind mount to hold genome data
docker volume create --opt type=none --opt o=bind --opt device=/media/user/data/genomes genomes

# Build copier
docker compose build copier

# Run copier to initialize a new workflow, answer the questions and pick human
docker compose run --rm -it copier copy /mnt/workflow /mnt/genomes/homo_sapiens

cd /media/user/data/genomes/homo_sapiens

# Build the workflow image
docker compose build workflow

# Run Snakemake within the workflow service to automatically run the workflow
docker compose run --rm -it workflow snakemake --use-conda -c2 -p -n
```
