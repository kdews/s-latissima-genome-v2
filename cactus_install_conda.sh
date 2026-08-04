#!/bin/bash

# Create cactus_gpu env with KegAlign (+ additional tools)
module purge # clear conflicting environments/modules
mamba create -n cactus_gpu -c bioconda kegalign-full
conda activate cactus_gpu

# Patch KegAlign scripts (see https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/Dockerfile.kegalign)
env_bin="$(dirname "$(which run_kegalign)")"
dos2unix "$env_bin"/run_kegalign "$env_bin"/*.py

# Download and extract Cactus pre-compiled binary
cd ~/bin
wget https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v3.2.1/cactus-bin-v3.2.1.tar.gz
tar -xzf cactus-bin-v3.2.1.tar.gz
cd cactus-bin-v3.2.1

# Cactus requires Python >3.9
which python
python --version

# Install virtualenv as Python module (if needed)
python3 -m pip install --user virtualenv

# Create virtual environment for Cactus
vnv="venv-cactus-v3.2.1"
python3 -m virtualenv "$vnv"
printf "export PATH=$(pwd)/bin:\$PATH\nexport PYTHONPATH=$(pwd)/lib:\$PYTHONPATH\nexport LD_LIBRARY_PATH=$(pwd)/lib:\$LD_LIBRARY_PATH\n" >> "$vnv"/bin/activate

# Activate Cactus virtual env
source "$vnv"/bin/activate

# Install Python module dependencies
python3 -m pip install -U setuptools pip wheel
python3 -m pip install -U .
python3 -m pip install -U -r ./toil-requirement.txt

# Fix KegAlign dependencies
python3 -m pip install -U bashlex mypy nvidia-ml-py flake8
mamba install -y -c conda-forge black zlib gxx=13 cmake

# Test Cactus install
cactus ./jobstore ./examples/evolverMammals.txt ./evolverMammals.hal
sbatch \
  /project2/noujdine_61/kdeweese/latissima/corteva_genome/s-latissima-genome-v2/prog_cactus.sbatch \
  ./examples/evolverMammals.txt

# Test Cactus Pangenome pipeline
cactus-pangenome ./js ./examples/evolverPrimates.txt --outDir primates-pg --outName primates-pg --reference simChimp --vcf --giraffe --gfa --gbz
sbatch \
  /project2/noujdine_61/kdeweese/latissima/corteva_genome/s-latissima-genome-v2/cactus_pangenome.sbatch \
  ./examples/evolverPrimates.txt simChimp
