#!/usr/bin/env python3
# Script to create cached conda environments for Snakemake to use.
# It only needs a valid Snakefile and the conda envs so it is a little
# nicer to use in containers than running Snakemake in create conda envs mode.

import shutil
import os
import subprocess
from pathlib import Path

from snakemake.deployment.conda import Env
from snakemake.workflow import Workflow

CODE_PATH = "/tmp/create_envs"
# Where the conda environments will be created
CONDA_ENV_PATH = "/mnt/snakemake-envs"

def create_environment(env_file, env_dir):
  """
  Create a conda environment defined in env_file, place it in env_dir
  """
  workflow = Workflow(os.path.join(CODE_PATH, "workflow", "Snakefile"), 
                      use_conda=True, 
                      conda_cleanup_pkgs="cache", 
                      conda_prefix=env_dir,
                      overwrite_configfiles=[])

  conda_env = Env(env_file,
    workflow=workflow,
    env_dir=env_dir)

  # Copy the conda environment file to the right place for Snakemake
  env_file_copy = "{}.yaml".format(conda_env.path)
  os.makedirs(conda_env.path, exist_ok=True)
  shutil.copy(env_file, env_file_copy)

  # Build the command to create the environment
  cmd = ["mamba","env","create",
        "--quiet",
        '-f={}'.format(env_file_copy),
        "--force",
        '-p {}'.format(conda_env.path)]

  # Run command to build environment
  subprocess.run(" ".join(cmd), shell=True, check=True)
  

def clean_conda_artifacts():
  """
  Clean up the conda artifacts and cache files.
  """
  subprocess.run(" ".join(["mamba","clean","--quiet","-y","--all"]), 
    shell=True, 
    check=True)

# Create all environment definitions in the envs directory
envs_path = Path(os.path.join(CODE_PATH, "workflow", "envs"))
conda_env_yamls = []
conda_env_yamls.extend([str(f.resolve()) for f in envs_path.rglob("*.yml") if f.is_file()])
conda_env_yamls.extend([str(f.resolve()) for f in envs_path.rglob("*.yaml") if f.is_file()])

# Loop and create environments
for env_file in conda_env_yamls:
  create_environment(
      env_file=env_file,
      env_dir=CONDA_ENV_PATH
    )

print("Created Snakemake conda environments for: {}".format(", ".join(conda_env_yamls)))

# Remove conda artifacts
clean_conda_artifacts()
