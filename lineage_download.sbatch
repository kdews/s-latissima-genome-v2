#!/bin/bash
#SBATCH -J lineage_download
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p main
#SBATCH --mem=4g

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

# Load conda
cond="$HOME/.conda_for_sbatch.sh"
if [[ -e "$cond" ]]
then
  source "$cond"
else
  echo "Error sourcing $cond"
  exit 1
fi

# Load BUSCO env
module purge
conda activate busco
busco --version || exit 1

# Download required BUSCO lineage datasets
lineages=(
  "eukaryota"
  "stramenopiles"
)
for i in "${!lineages[@]}"
do
  odb="$(
    busco --list-datasets \
      | grep "${lineages[i]}" \
      | sed -r "s/.*(${lineages[i]}_odb.+) .*/\1/g"
  )"
  busco --download "$odb"
done
