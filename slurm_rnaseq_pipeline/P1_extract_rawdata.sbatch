#!/bin/bash -l
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 80G
#SBATCH -t 2:00:00

gzip -d $1