args <- commandArgs(TRUE)

alignment_result <- args[1]
output_result <- args[2]
meta_file <- args[3]

source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')

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
