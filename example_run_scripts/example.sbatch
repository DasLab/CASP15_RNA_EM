#!/bin/bash
#SBATCH -J XXX
#SBATCH -o run_files/XXX.out
#SBATCH -e run_files/XXX.err
#SBATCH -p biochem
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -N 1

ulimit -s unlimited; lga
module load gcc/12.1.0

