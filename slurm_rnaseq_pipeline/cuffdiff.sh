#! /bin/bash

#SBATCH -p Interactive  # use the Lewis partition
#SBATCH -J cpu_sub_8n_22h  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 0-4:00  # two hour time limit
#SBATCH --mem=94G
#SBATCH -N 1  # number of nodes
#SBATCH -n 16  # number of cores (AKA tasks)
export  HDF5_USE_FILE_LOCKING='FALSE'
# Commands here run only on the first core
echo "$*"
module load cufflinks/cufflinks-2.2.1



cuffdiff -p 16 -o $1 --use-sample-sheet $2 $3