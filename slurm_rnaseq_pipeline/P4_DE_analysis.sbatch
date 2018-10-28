#!/bin/bash -l
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 60G
#SBATCH -t 2:00:00

Rscript /storage/htc/bdm/ccm3x/rnaseq_pipeline/de_analysis.R $1 $2 