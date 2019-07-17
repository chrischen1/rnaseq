library(org.Mm.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(pheatmap)

output_path = '~/Dropbox/MU/workspace/will_rnaseq/final/'
cnt_result <- read.csv('~/Dropbox/MU/workspace/liu_rnaseq/counts.csv',as.is = T,row.names = 1)
de_result <- read.csv('~/Dropbox/MU/workspace/liu_rnaseq/de_results.csv',as.is = T,row.names = 1)
meta <- read.csv('~/Dropbox/MU/workspace/liu_rnaseq/meta_info.csv',as.is = T,row.names = 1)

rownames(meta) <- gsub('^D','',rownames(meta))
#heatmap
genes_ens <- de_result[de_result$fdr<0.05,]
cnt_result <- cnt_result[rownames(genes_ens),]
png(paste(output_path,'heatmap.png',sep = ''))
pheatmap(log2(cnt_result+1),show_rownames = F,cluster_cols = F,annotation_col = meta[,2,drop=F],scale = 'column')
dev.off()

eg = bitr(rownames(genes_ens), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

genes_ens1 <- genes_ens[eg$ENSEMBL,]
genes_ens1 <- cbind(genes_ens1,'ENTREZID'=eg$ENTREZID)

# Pathway Enrichment Analysis
x <- enrichPathway(gene=genes_ens1$ENTREZID,organism = 'mouse',pvalueCutoff=0.05, readable=T)

for(i in c(seq(20,153,20),153)){
  png(paste(output_path,'pathway1_',i,'.png',sep = ''),width = 1280+i*15,height = 960+i*10)
  barplot(x, showCategory=i)
  dev.off()
  
  png(paste(output_path,'pathway2_',i,'.png',sep = ''),width = 1280+i*15,height = 960+i*10)
  dotplot(x, showCategory=i)
  dev.off()

}


cnetplot(x, categorySize="pvalue", foldChange=2^(genes_ens1$logFC),showCategory=2)

# GSEA
geneList <- 2^(genes_ens1$logFC)
names(geneList) <- genes_ens1$ENTREZID
geneList <- sort(geneList,decreasing = T)
y <- gsePathway(geneList, nPerm=10000,organism='mouse',pvalueCutoff=0.2,pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(y)

gseaplot(y, geneSetID = res$ID[1])
