#!/bin/bash
#SBATCH -J sdr_annot
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p htcondor
#SBATCH --time=01:00:00
#SBATCH --mem=40g
#SBATCH --cpus-per-task=16

