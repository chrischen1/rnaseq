final_path = '/storage/htc/bdm/Collaboration/David_rnaseq/final_result'

all_tfs <- read.csv('/storage/htc/bdm/Collaboration/David_rnaseq/tf_list.csv',header = F,as.is = T)$V1
sig_genes <- c()
for(i in list.files(final_path,pattern = 'DE_.+',full.names = T)){
  de_result <- read.csv(i,as.is = T,row.names = 1)
  sig_genes <- c(sig_genes,rownames(de_result)[de_result$FDR<0.05])
}
sig_genes <- unique(sig_genes)
sig_tfs <- intersect(all_tfs,sig_genes)
write.table(sig_genes,paste0(final_path,'sig_genes.csv'),row.names = F,col.names = F,quote = F)
write.table(sig_tfs,paste0(final_path,'sig_tfs.csv'),row.names = F,col.names = F,quote = F)

# rpkm_data <- read.csv('/storage/htc/bdm/Collaboration/David_rnaseq/final_result/rpkm.csv',row.names = 1,check.names = F)
# rpkm_sig <- rpkm_data[sig_genes,]
# log2_rpkm_sig <- log2(rpkm_sig+0.001)



