#!/bin/bash
#SBATCH -J hal_extract
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --mem=15gb
#SBATCH --time=05:00:00

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
outdir="synteny_PAFs" # output directory name
mergePaf="merged_synteny.paf"

# Help message
show_help () {
  echo "\
Extracts individual and merged PAF files from chromosome alignments.
Requires output directory from cactus-pangenome (cactusDir) as input.
Usage: sbatch $job_name.sbatch <cactusDir> <ref_label> <query_label>

Required input:
  cactusDir         path to directory containing output of cactus-pangenome
  ref_label         assembly label to be used as reference
  query_label             assembly label to be used as query_label

Output:
  <cactusDir>/$outdir/  output directory for PAFs

Requires:
 - virtualenv (https://virtualenv.pypa.io/en/latest)
 - Toil (https://toil.readthedocs.io/en/latest)
 - Cactus (https://github.com/ComparativeGenomicsToolkit/cactus)
   - contains HALtools distribution (http://github.com/glennhickey/hal)"
}
# Catch help flag and missing argument(s)
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  show_help
  exit 0
elif (( $# < 3 ))
then
  echo "Error: 3 positional arguments required."
  echo
  show_help
  exit 1
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
cactusDir="${1:-cactus_pangenome_results/01-28-26_6039802}"
ref_label="${2:-corteva_v2}"
query_label="${3:-roscoff_v2}"
if [[ -d "$cactusDir" ]]
then
  # Redirect output
  echo "Using results in: $cactusDir"
  outdir="$cactusDir/$outdir"
  mergePaf="$outdir/$mergePaf"
  mkdir -p "$outdir"
  echo "Directing output to: $outdir"
else
  echo "Error: No results directory found in $cactusDir"
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
    --targetGenome "$ref_label"
    --queryGenome "$query_label"
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
