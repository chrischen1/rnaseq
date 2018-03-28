#!/bin/bash -l
#SBATCH -p Lewis
#SBATCH -o $3.out
#SBATCH -e $3.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 80G
#SBATCH -t 2:00:00

/group/birchler-cheng/tools/FASTX_Toolkit/fastq_quality_filter -Q 33  -q 20 -p 80  -i $1 -o $2