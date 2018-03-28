#!/bin/bash -l
#SBATCH -p Lewis
#SBATCH -o $3.out
#SBATCH -e $3.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 60G
#SBATCH -t 2:00:00

Rscript de_analysis.R $1 $2