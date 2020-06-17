#! /bin/bash

#SBATCH -p Interactive  # use the Lewis partition
#SBATCH -J index  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 0-2:00  # two hour time limit
#SBATCH --mem=94G
#SBATCH -N 1  # number of nodes
#SBATCH -n 16  # number of cores (AKA tasks)
export  HDF5_USE_FILE_LOCKING='FALSE'
# Commands here run only on the first core
echo "$*"

module load star/star-2.7.0e
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir $1 --genomeFastaFiles $2 --sjdbGTFfile $3 --sjdbOverhang 100 