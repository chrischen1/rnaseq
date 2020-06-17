#! /bin/bash

#SBATCH -p Interactive  # use the Lewis partition
#SBATCH -J cpu_sub_8n_22h  # give the job a custom name
#SBATCH -o output.fa  # give the job output a custom name
#SBATCH -t 0-1:00  # two hour time limit
#SBATCH --mem=24G
#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $1
