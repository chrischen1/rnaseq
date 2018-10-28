args <- commandArgs(TRUE)
source(args[1])

script_dir = '/storage/htc/bdm/ccm3x/rnaseq_pipeline/'
source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')
library(edgeR)

rawdata_path <- fastq_dir
trim_data_path <- paste(project_dir,'trimdata/',sep = '')
alignment_result <- paste(project_dir,'alignment_result/',sep = '')
log_path <- paste(project_dir,'logs/',sep = '')
output_result <- paste(project_dir,'final_result/',sep = '')

dir.create(project_dir,showWarnings = F)
dir.create(trim_data_path,showWarnings = F)
dir.create(alignment_result,showWarnings = F)
dir.create(output_result,showWarnings = F)
dir.create(log_path,showWarnings = F)
job_name <- substr(gsub('.+/(.+)/','\\1',project_dir),1,8)

setwd(log_path)
# After downloading, extract results
print(paste(Sys.time(),'extracting results'))
for(i in list.files(rawdata_path,pattern = '.+.fastq.gz')){
  file2extract <- paste(rawdata_path,i,sep = '')
  system(paste('sbatch -J ',job_name,' ',script_dir,'P1_extract_rawdata.sbatch ',file2extract,sep = ''))
}
Sys.sleep(60)
while(slurm_running(job_name)) {
  Sys.sleep(180)
}

# After extraction, trim results
print(paste(Sys.time(),'trimming results'))
for(i in list.files(rawdata_path,pattern = '.+.fastq$')){
  infile <- paste(rawdata_path,i,sep = '')
  outfile <- paste(trim_data_path,i,sep = '') 
  system(paste('sbatch -J ',job_name,' ',script_dir,'P2_trim_data.sbatch ',infile,' ',outfile,sep = ''))
}
Sys.sleep(60)
while(slurm_running(job_name)) {
  Sys.sleep(180)
}

# After trimming, start alignment
print(paste(Sys.time(),'start alignment'))
for(i in list.files(trim_data_path,pattern = '.+.fastq.*')){
  out_path_i <- paste(alignment_result,gsub('.fastq.*','',i),'/',sep = '')
  dir.create(out_path_i,showWarnings = F)
  infile <- paste(trim_data_path,i,sep = '')
  system(paste('sbatch -J ',job_name,' ',script_dir,'P3_alignment_with_star.sbatch ',ref_genome,' ',infile,' ',out_path_i,sep = ''))
}
Sys.sleep(60)
while(slurm_running(job_name)) {
  Sys.sleep(300)
}

#After alignment, start differential expression analysis
# combine counts
cnt <- NULL
cnt_colnames <- NULL
all_results <- list.dirs(alignment_result,recursive = F)
for(i in all_results){
  new_cnt <- read.delim(paste(i,'/ReadsPerGene.out.tab',sep = ''),row.names = 1)
  if(is.null(cnt)){
    cnt <- new_cnt[-(1:3),1,drop=FALSE]
  }else{
    cnt <- cbind(cnt,new_cnt[rownames(cnt),1])
  }
  
  cnt_colnames <- c(cnt_colnames,gsub('.+//','',i))
}
colnames(cnt) <- cnt_colnames
write.csv(cnt,paste(output_result,'counts.csv',sep = ''))

# get RPKM
if(gtf_file != ''){
  gtf <- read.delim(gtf_file,as.is = T,comment.char = '#',header = F,sep = '\t')
  gtf_gene <- gtf[gtf$V3=='gene',]
  gene_length <- as.numeric(gtf_gene$V5)-as.numeric(gtf_gene$V4)
  names(gene_length) <- gsub('gene_id ([A-Za-z0-9]+); .+','\\1',gtf_gene$V9)
  rpkm_data <- rpkm(cnt,gene.length = gene_length[rownames(cnt)])
  write.csv(rpkm_data,paste(output_result,'rpkm.csv',sep = ''))
}

# differential expression analysis
if(meta_file!=''){
  print(paste(Sys.time(),'start differential expression analysis'))
  system(paste('sbatch -J ',job_name,' ',script_dir,'P4_DE_analysis.sbatch ',output_result,' ',meta_file,sep = ''))
  Sys.sleep(60)
  while(slurm_running(job_name)) {
    Sys.sleep(300)
  }
}

#RNASeq pipeline finished
print(paste(Sys.time(),'RNASeq pipeline finished'))

