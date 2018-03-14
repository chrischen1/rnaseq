download_script = '/storage/htc/bdm/Collaboration/liu_rnaseq/src/download_fastq.sh'
ref_genome = '/storage/htc/bdm/ref_genomes/mouse/star_ref'
project_dir = '/storage/htc/bdm/Collaboration/liu_rnaseq/'
meta_file = '/storage/htc/bdm/Collaboration/liu_rnaseq/src/meta_info.csv'

rawdata_path <- paste(project_dir,'rawdata/',sep = '')
trim_data_path <- paste(project_dir,'trimdata/',sep = '')
alignment_result <- paste(project_dir,'alignment_result/',sep = '')
output_result <- paste(project_dir,'final_result/',sep = '')

dir.create(rawdata_path,showWarnings = F)
dir.create(trim_data_path,showWarnings = F)
dir.create(alignment_result,showWarnings = F)

source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')
# download fastq
setwd(rawdata_path)
system(paste('cd',rawdata_path))
system(paste('bash',download_script))

# After downloading, extract results
for(i in list.files(rawdata_path,pattern = '.+.fastq.gz')){
  system(paste('srun -n1 -p Lewis -t 2:00:00 --mem 40G gzip -d ',rawdata_path,i,' &',sep = ''))
}
Sys.sleep(60)
while(slurm_running()) {
  Sys.sleep(180)
}

# After extraction, trim results
for(i in list.files(rawdata_path,pattern = '.+.fastq$')){
  system(paste('srun -n1 -p Lewis -t 2:00:00 --mem 40G /group/birchler-cheng/tools/FASTX_Toolkit/fastq_quality_filter -Q 33  -q 20 -p 80  -i ',rawdata_path,i,' -o ',trim_data_path,i,' &',sep = ''))
}
Sys.sleep(60)
while(slurm_running()) {
  Sys.sleep(180)
}

# After trimming, start alignment
system('module load star/star-2.5.2b')
for(i in list.files(trim_data_path,pattern = '.+.fastq$')){
  out_path_i <- paste(alignment_result,gsub('.fastq','',i),sep = '')
  dir.create(out_path_i,showWarnings = F)
  system(paste('srun -n1 -p Lewis -t 4:00:00 --mem 80G STAR --runThreadN 1 --quantMode GeneCounts --genomeDir ',ref_genome,' --readFilesIn ',trim_data_path,i,' --outFileNamePrefix ',out_path_i,'/ &',sep = ''))
}
Sys.sleep(60)
while(slurm_running()) {
  Sys.sleep(300)
}

#After alignment, start differential expression analysis
cnt <- NULL

grp <- read.csv(meta_file,as.is=T,row.names=1)
de_result <- edgeR_wrapper(cnt,grp)