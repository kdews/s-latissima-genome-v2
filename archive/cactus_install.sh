#!/bin/bash

# Download and extract Cactus pre-compiled binary
cd ~/bin
wget https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v3.1.3/cactus-bin-v3.1.3.tar.gz
tar -xvf cactus-bin-v3.1.3.tar.gz
cd cactus-bin-v3.1.3

# Clear conflicting environments/modules
conda deactivate
module purge

# Load Python 3.11 (Cactus requires Python >3.9)
module load ver/2506 gcc/14.3.0 python/3.11.14
# Install virtualenv as Python module (if needed)
python3 -m pip install --user virtualenv

# Create virtual environment for Cactus
python3 -m virtualenv venv-cactus-v3.1.3
printf "export PATH=$(pwd)/bin:\$PATH\nexport PYTHONPATH=$(pwd)/lib:\$PYTHONPATH\nexport LD_LIBRARY_PATH=$(pwd)/lib:\$LD_LIBRARY_PATH\n" >> venv-cactus-v3.1.3/bin/activate

# Activate Cactus virtual env
source venv-cactus-v3.1.3/bin/activate

# Install Python module dependencies
python3 -m pip install -U setuptools pip wheel
python3 -m pip install -U .
python3 -m pip install -U -r ./toil-requirement.txt

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
