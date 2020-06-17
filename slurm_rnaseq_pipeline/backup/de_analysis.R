args <- commandArgs(TRUE)

output_result <- args[1]
meta_file <- args[2]

source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')

cnt <- read.csv(paste(output_result,'counts.csv',sep = ''),row.names = 1,check.names = F)
grp <- read.csv(meta_file,as.is=T,row.names=1)

cnt <- cnt[,rownames(grp)]
cnt <- cnt[apply(cnt,1,sum)>0,]
de_result <- edgeR_wrapper(cnt,grp)

for(i in 1:length(de_result)){
  write.csv(de_result[[i]],paste(output_result,'DE_',names(de_result)[i],'.csv',sep = ''))
}


