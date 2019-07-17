# This is for single-end and adapter trimmed reads, just provide the project folder
# The raw fastq files should be put at project_folder/raw

# Usage: Rscript rnaseq_pipeline.R project_folder [gtf_file] [meta_file]

args <- commandArgs(TRUE)
# args = c('/storage/htc/bdm/Collaboration/Will_rnaseq/',
#          '/storage/htc/bdm/ref_genomes/mouse/Mus_musculus.GRCm38.97.gtf',
#          '/storage/htc/bdm/Collaboration/Will_rnaseq/group.tsv')

job_prefix = '#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 0-8:00  # two hour time limit
#SBATCH --mem=48G
#SBATCH -N 1  # number of nodes
#SBATCH -n 8  # number of cores (AKA tasks)

'
qc_commands = '
# QC with fastqc
module load fastqc/fastqc-0.11.7
fastqc '

trim_commands = '
# Trimming
module load trimmomatic/trimmomatic-0.38
trimmomatic SE -phred33 -threads 8'
trim_configs = 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 '

alignment_commands = '
# Alignment
module load star/star-2.7.0e
STAR --runThreadN 8 --genomeDir /storage/htc/bdm/ref_genomes/mouse/star_ref/ --quantMode GeneCounts --readFilesCommand zcat '

project_dir <- args[1]
gtf_file <- args[2]
meta_file <- args[3]

rawdata_path <- paste0(project_dir,'/raw/')
qc_data_path <- paste0(project_dir,'/fastqc/')
trim_data_path <- paste0(project_dir,'/trim/')
alignment_result <- paste0(project_dir,'/alignment/')
script_path <- paste0(project_dir,'/scripts/')
log_path <- paste0(project_dir,'/logs/',sep = '')
output_result <- paste(project_dir,'final/',sep = '')

dir.create(project_dir,showWarnings = F)
dir.create(qc_data_path,showWarnings = F)
dir.create(trim_data_path,showWarnings = F)
dir.create(alignment_result,showWarnings = F)
dir.create(script_path,showWarnings = F)
dir.create(output_result,showWarnings = F)
dir.create(log_path,showWarnings = F)

setwd(log_path)

# for(i in list.files(rawdata_path,pattern = '.+.trimmed.fastq.gz')){
#   sample_name <- gsub('_R1_001.trimmed.fastq.gz','',i)
#   command_qc <- paste0(job_prefix,qc_commands,' ',rawdata_path,i,' -o ',qc_data_path)
#   writeChar(command_qc,paste0(script_path,sample_name,'.sh'),eos = '\n')
# }
  
# After qc, trim results and perform alignment
for(i in list.files(rawdata_path,pattern = '.+.trimmed.fastq.gz')){
  sample_name <- gsub('_R1_001.trimmed.fastq.gz','',i)
  command_trim <- paste0(trim_commands,' ',rawdata_path,i,' ',trim_data_path,i,' ',trim_configs)
  align_out_path <- paste0(alignment_result,'/',sample_name,'/')
  dir.create(align_out_path,showWarnings = F)
  command_align <- paste0(alignment_commands,' --readFilesIn ',trim_data_path,i,' --outFileNamePrefix ',align_out_path)
  writeLines(paste(job_prefix,command_trim,command_align,sep = '\n'),paste0(script_path,sample_name,'.sh'))
}
for(i in list.files(script_path,full.names = T)){
  system(paste('sbatch',i))
}

# combine alignment qc and counts
qc_table <- cnt <- NULL
cnt_colnames <- NULL
all_results <- list.dirs(alignment_result,recursive = F)
for(i in all_results){
  new_cnt <- read.delim(paste(i,'/ReadsPerGene.out.tab',sep = ''),row.names = 1)
  new_qc <- read.delim(paste(i,'/Log.final.out',sep = ''),header = F,sep = '\t')
  if(is.null(qc_table)){
    qc_table <- new_qc
    cnt <- new_cnt[-(1:4),1,drop=FALSE]
  }else{
    qc_table <- cbind(qc_table,new_qc$V2)
    cnt <- cbind(cnt,new_cnt[rownames(cnt),1])
  }
  cnt_colnames <- c(cnt_colnames,gsub('.+//','',i))
}
colnames(qc_table) <- cnt_colnames
colnames(cnt) <- cnt_colnames

write.csv(qc_table,paste(output_result,'alignment_qc.csv',sep = ''))
write.csv(cnt,paste(output_result,'counts.csv',sep = ''))

library(edgeR)
source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')
# get RPKM
if(gtf_file != ''){
  rpkm_data <- get_rpkm(cnt,gtf_file)
  write.csv(rpkm_data,paste(output_result,'rpkm.csv',sep = ''))
}

# differential expression analysis
if(meta_file!=''){
  grp <- read.csv(meta_file,as.is=T,row.names=1)
  de_result <- edgeR_wrapper(cnt,grp)
  
  de_results_path <- paste0(output_result,'/Differential_expression/')
  dir.create(de_results_path)
  for(i in names(de_result)){
    write.csv(de_result[[i]],paste0(de_results_path,i,'.csv'))
  }
}

#RNASeq pipeline finished
print(paste(Sys.time(),'RNASeq pipeline finished'))

