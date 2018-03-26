rawdata_path <- paste(project_dir,'rawdata/',sep = '')
trim_data_path <- paste(project_dir,'trimdata/',sep = '')
alignment_result <- paste(project_dir,'alignment_result/',sep = '')
output_result <- paste(project_dir,'final_result/',sep = '')

dir.create(rawdata_path,showWarnings = F)
dir.create(trim_data_path,showWarnings = F)
dir.create(alignment_result,showWarnings = F)
dir.create(output_result,showWarnings = F)
job_name <- gsub('/','_',project_dir)
source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')

# download fastq
setwd(rawdata_path)
system(paste('cd',rawdata_path))
system(paste('bash',download_script))

# After downloading, extract results
for(i in list.files(rawdata_path,pattern = '.+.fastq.gz')){
  system(paste('srun -n1 -p Lewis -t 2:00:00 --mem 40G -J ',job_name,' gzip -d ',rawdata_path,i,' &',sep = ''))
}
Sys.sleep(60)
while(slurm_running(job_name)) {
  Sys.sleep(180)
}

# After extraction, trim results
for(i in list.files(rawdata_path,pattern = '.+.fastq$')){
  system(paste('srun -n1 -p Lewis -t 2:00:00 --mem 40G -J ',job_name,' /group/birchler-cheng/tools/FASTX_Toolkit/fastq_quality_filter -Q 33  -q 20 -p 80  -i ',rawdata_path,i,' -o ',trim_data_path,i,' &',sep = ''))
}
Sys.sleep(60)
while(slurm_running(job_name)) {
  Sys.sleep(180)
}

# After trimming, start alignment
system('module load star/star-2.5.2b')
for(i in list.files(trim_data_path,pattern = '.+.fastq$')){
  out_path_i <- paste(alignment_result,gsub('.fastq','',i),sep = '')
  dir.create(out_path_i,showWarnings = F)
  system(paste('srun -n1 -p Lewis -t 6:00:00 --mem 110G -J ',job_name,' STAR --runThreadN 1 --quantMode GeneCounts --genomeDir ',ref_genome,' --readFilesIn ',trim_data_path,i,' --outFileNamePrefix ',out_path_i,'/ &',sep = ''))
}
Sys.sleep(60)
while(slurm_running(job_name)) {
  Sys.sleep(300)
}

#After alignment, start differential expression analysis
cnt <- NULL
cnt_colnames <- NULL
all_results <- list.dirs(alignment_result)[-1]
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
# cnt <- gene_matrix_conversion(cnt,dataset = 'mmusculus_gene_ensembl',gene_id2 = 'mgi_symbol')

grp <- read.csv(meta_file,as.is=T,row.names=1)
cnt <- cnt[,rownames(grp)]
de_result <- edgeR_wrapper(cnt,grp)
if(is.list(de_result)){
  write.csv(de_result$pmat,paste(output_result,'pval.csv',sep = ''))
  write.csv(de_result$fdr_mat,paste(output_result,'fdr.csv',sep = ''))
  write.csv(de_result$logFC,paste(output_result,'logFC.csv',sep = ''))
}else{
  write.csv(de_result,paste(output_result,'de_results.csv',sep = ''))
}
colnames(cnt) <- gsub('^D','',colnames(cnt) )
write.csv(cnt,paste(output_result,'counts.csv',sep = ''))
