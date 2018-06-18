args <- commandArgs(TRUE)
source(args[1])

source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')
script_dir = '/storage/htc/bdm/ccm3x/rnaseq_pipeline/'

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

# download fastq
setwd(rawdata_path)
system(paste('cd',rawdata_path))
if(download_script!=''){
  print(paste(Sys.time(),'downloading fastq files'))
  system(paste('bash',download_script))
}else{
  print(paste(Sys.time(),'copying fastq files'))
  system(paste('cp ',fastq_dir,'*.fastq.gz ./',sep = ''))
  system(paste('cp ',fastq_dir,'*.fastq ./',sep = ''))
}

setwd(log_path)
# After downloading, extract results
print(paste(Sys.time(),'extracting results'))
for(i in list.files(rawdata_path,pattern = '.+.fastq.gz')){
  file2extract <- paste(rawdata_path,i,sep = '')
  system(paste('gzip -d',file2extract))
}

# After extraction, trim results
print(paste(Sys.time(),'trimming results'))
for(i in list.files(rawdata_path,pattern = '.+.fastq$')){
  infile <- paste(rawdata_path,i,sep = '')
  outfile <- paste(trim_data_path,i,sep = '') 
  system(paste('/group/birchler-cheng/tools/FASTX_Toolkit/fastq_quality_filter -Q 33  -q 20 -p 80  -i',infile,'-o',outfile))
}

# After trimming, start alignment
print(paste(Sys.time(),'start alignment'))
system('module load star/star-2.5.2b')
for(i in list.files(trim_data_path,pattern = '.+.fastq$')){
  out_path_i <- paste(alignment_result,gsub('.fastq','',i),'/',sep = '')
  dir.create(out_path_i,showWarnings = F)
  infile <- paste(trim_data_path,i,sep = '')
  system(paste('STAR --runThreadN 1 --quantMode GeneCounts --genomeDir',ref_genome,'--readFilesIn',infile,'--outFileNamePrefix',out_path_i))
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


if(meta_file!=''){
  print(paste(Sys.time(),'start differential expression analysis'))
  system(paste('Rscript /storage/htc/bdm/ccm3x/rnaseq_pipeline/de_analysis.R',alignment_result,output_result,meta_file))
}

#RNASeq pipeline finished
print(paste(Sys.time(),'RNASeq pipeline finished.'))

