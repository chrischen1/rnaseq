cnt <- read.csv('~/Dropbox/MU/workspace/will_rnaseq/final/counts.csv',row.names = 1,check.names = F)
grp <- read.csv('~/Dropbox/MU/workspace/will_rnaseq/group.csv',row.names = 1,check.names = F)
rpkm_data <- read.csv('~/Dropbox/MU/workspace/will_rnaseq/final/rpkm.csv',row.names = 1,check.names = F)

source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')

rownames(grp) <- paste0('S',rownames(grp))
colnames(rpkm_data) <- colnames(cnt) <- gsub('.+_','',colnames(cnt))
cnt <- cnt[,rownames(grp)]
rpkm_data <- rpkm_data[,rownames(grp)]

# heatmap
rpkm_data1 <- rpkm_data[apply(rpkm_data, 1, function(x)(sum(x>0)))>22,]
my_sample_col <- data.frame(group = grp$group,condition = grp$condition)
row.names(my_sample_col) <- rownames(grp)
library(pheatmap)
pheatmap(log2(rpkm_data1+0.01),show_rownames = F,scale = 'column', annotation_col = my_sample_col)

# differential expression analysis
de_result <- edgeR_wrapper(cnt,grp)
de_results_path <- '~/Dropbox/MU/workspace/will_rnaseq/final/DE/'
for(i in names(de_result)){
  write.csv(de_result[[i]],paste0(de_results_path,i,'.csv'))
}

grp1 <- grp[(grp$condition=='control') | (grp$condition=='control_TMT'),]
grp1$control[grp1$condition=='control_TMT'] <- FALSE
grp1$group[grp1$condition=='control_TMT'] <- 'normal'
de_result1 <- edgeR_wrapper(cnt[,rownames(grp1)],grp1)
for(i in names(de_result1)){
  write.csv(de_result1[[i]],paste0(de_results_path,i,'.csv'))
}

# PCA 
library(ggfortify)
pca <- prcomp(t(rpkm_data1), scale.=TRUE)


gr <- factor(grp$condition)
pca3d(pca, group=gr, show.ellipses=TRUE,ellipse.ci=0.75, show.plane=FALSE,show.group.labels = T)
autoplot(pca, data = data.frame(grp), colour = 'condition', frame = TRUE, frame.type = 'norm')


# enrichment analysis
library(org.Mm.eg.db)
library(clusterProfiler)
library(pathview)

genes_slt1 <- de_result$normal.BME_50
genes_slt2 <- rownames(genes_slt1)[genes_slt1$FDR < 0.05]
result<- enrichKEGG(genes_slt,organism = 'mmu',keyType = 'uniprot',pAdjustMethod='none',pvalueCutoff = 1.1,qvalueCutoff = 1.1)





















