#!/bin/bash
#SBATCH -J hal_extract
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --mem=15gb
#SBATCH --time=05:00:00

# Script variables
# Extract job name from script
if [[ -z "${SLURM_JOB_NAME}" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
outdir="synteny_PAFs" # output directory name
mergePaf="merged_synteny.paf"

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Extracts individual and merged PAF files from chromosome alignments.
Requires output directory from cactus-pangenome (cactusDir) as input.
Usage: sbatch hal_extract.sbatch <cactusDir> <reference> <query>

Required input:
  cactusDir         path to directory containing output of cactus-pangenome
  reference         assembly ID to be used as reference
  query             assembly ID to be used as query

Output:
  <cactusDir>/$outdir/  output directory for PAFs

Requires:
 - virtualenv (https://virtualenv.pypa.io/en/latest)
 - Toil (https://toil.readthedocs.io/en/latest)
 - Cactus (https://github.com/ComparativeGenomicsToolkit/cactus)"
  exit 0
fi

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"

# Load Python module and Cactus virtualenv
module purge
module load ver/2506 gcc/14.3.0 python/3.11.14
env_path="/home1/kdeweese/bin/cactus-bin-v3.1.3/venv-cactus-v3.1.3/bin/activate"
if [[ -e "$env_path" ]]
then
  source "$env_path"
else
  echo "Error on source of $env_path"
  exit 1
fi

# # Test hal2maf
# which hal2maf
# echo
# Test halSynteny
which halSynteny
echo
# Test Python (for psl2paf conversion)
which python3
python3 --version
echo

# Input
# cactusDir=cactus_pangenome_01-28-26_6039802
# reference=corteva_v2
# query=roscoff_v2
if [[ "$#" -eq 3 ]]
then
  cactusDir="$1"
  reference="$2"
  query="$3"
  if [[ -d "$cactusDir" ]]
  then
    # Redirect output
    echo "Using results in: $cactusDir"
    outdir="$cactusDir/$outdir"
    mergePaf="$outdir/$mergePaf"
    mkdir -p "$outdir"
    echo "Directing output to: $outdir"
  else
    echo "Error: no results directory found ($cactusDir)."
    exit 1
  fi
else
  echo "Error: 3 positional arguments expected."
  exit 1
fi

# Input HAL files
mapfile -t inHals < <(ls "$cactusDir/chrom-alignments/"*".hal")
outPafs=()
for inHal in "${inHals[@]}"
do
  # Prefix output files by inHal name
  outname="$(basename "${inHal%.*}")" # output file prefix
  outPsl="$outdir/${outname}_synteny.psl"
  outPaf="$outdir/${outname}_synteny.paf"

  # # Run hal2maf
  # cmd=(
  #   hal2maf
  #   --hdf5InMemory
  #   "$inHal"
  #   "$outMaf"
  # )
  # echo "${cmd[*]}"
  # "${cmd[@]}"

  # Run halSynteny
  cmd=(
    halSynteny
    --targetGenome "$reference"
    --queryGenome "$query"
    "$inHal"
    "$outPsl"
  )
  echo "${cmd[*]}"
  "${cmd[@]}"
  
  # Convert PSL to PAF
  if [[ -f "$outPsl" ]]
  then
    echo "python3 psl2paf.py $outPsl $outPaf"
    python3 psl2paf.py "$outPsl" "$outPaf"
    if [[ -f "$outPaf" ]]
    then
      outPafs+=("$outPaf")
    fi
  fi
  echo
done

# Merge per-chromosome PAFs
if [[ -n "${outPafs[*]}" ]] 
then
  echo "Merging synteny PAFs into: $mergePaf"
  echo "cat ${outPafs[*]} > $mergePaf"
  cat "${outPafs[@]}" > "$mergePaf"
else
  echo "Error: no PAFs created. Check halSynteny + psl2paf steps."
  exit 1
fi

# Print date and time
echo
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
