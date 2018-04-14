args <- commandArgs(TRUE)

alignment_result <- args[1]
output_result <- args[2]
meta_file <- args[3]

source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')

cnt <- read.csv(paste(output_result,'counts.csv',sep = ''),row.names = 1)
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
