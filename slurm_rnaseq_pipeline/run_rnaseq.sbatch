#!/bin/bash -l
#SBATCH -p Lewis
#SBATCH -o %J.out
#SBATCH -e %J.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 50G
#SBATCH -t 16:00:00

Rscript rnaseq_pipeline.R $1