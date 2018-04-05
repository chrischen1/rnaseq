script_dir = '/storage/htc/bdm/ccm3x/rnaseq_pipeline/'
args <- commandArgs(TRUE)
source(args[1])
source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')

rawdata_path <- paste(project_dir,'rawdata/',sep = '')
trim_data_path <- paste(project_dir,'trimdata/',sep = '')
alignment_result <- paste(project_dir,'alignment_result/',sep = '')
log_path <- paste(project_dir,'logs/',sep = '')
output_result <- paste(project_dir,'final_result/',sep = '')

dir.create(project_dir,showWarnings = F)
dir.create(rawdata_path,showWarnings = F)
dir.create(trim_data_path,showWarnings = F)
dir.create(alignment_result,showWarnings = F)
dir.create(output_result,showWarnings = F)
dir.create(log_path,showWarnings = F)
job_name <- substr(gsub('.+/(.+)/','\\1',project_dir),1,8)

# download fastq
setwd(rawdata_path)
system(paste('cd',rawdata_path))
if(download_script!=''){
  print(paste(Sys.time(),'downloading fastq files'))
  system(paste('bash',download_script))
}else{
  print(paste(Sys.time(),'copying fastq files'))
  system(paste('cp ',fastq_dir,'*.fastq.gz ./',sep = ''))
}

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
system('module load star/star-2.5.2b')
for(i in list.files(trim_data_path,pattern = '.+.fastq$')){
  out_path_i <- paste(alignment_result,gsub('.fastq','',i),'/',sep = '')
  dir.create(out_path_i,showWarnings = F)
  infile <- paste(trim_data_path,i,sep = '')
  system(paste('sbatch -J ',job_name,' ',script_dir,'P3_alignment_with_star.sbatch ',ref_genome,' ',infile,' ',out_path_i,sep = ''))
}
Sys.sleep(60)
while(slurm_running(job_name)) {
  Sys.sleep(300)
}

#After alignment, start differential expression analysis
if(meta_file!=''){
  print(paste(Sys.time(),'start differential expression analysis'))
  system(paste('sbatch -J ',job_name,' ',script_dir,'P4_DE_analysis.sbatch ',alignment_result,' ',output_result,' ',meta_file,sep = ''))
  Sys.sleep(60)
  while(slurm_running(job_name)) {
    Sys.sleep(300)
  }
}

#RNASeq pipeline finished
print(paste(Sys.time(),'RNASeq pipeline finished'))

