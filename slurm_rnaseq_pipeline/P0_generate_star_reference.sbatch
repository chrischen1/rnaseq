#!/bin/bash -l
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 110G
#SBATCH -t 10:00:00

STAR --runThreadN 1 --runMode genomeGenerate --genomeDir $1 --genomeFastaFiles $2 --sjdbGTFfile &3 --sjdbOverhang 100 