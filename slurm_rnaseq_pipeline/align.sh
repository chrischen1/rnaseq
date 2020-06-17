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
module load star/star-2.7.0e
# Commands with srun will run on all cores in the allocation
STAR --runThreadN 16 --genomeDir $1 --readFilesIn $2 --outFileNamePrefix $3 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical