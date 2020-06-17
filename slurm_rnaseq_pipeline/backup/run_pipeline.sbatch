#!/bin/bash -l
#SBATCH -p Lewis
#SBATCH -o %J.out
#SBATCH -e %J.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 23:00:00

module load star/star-2.5.2b
Rscript /storage/htc/bdm/ccm3x/rnaseq_pipeline/rnaseq_pipeline.R $1