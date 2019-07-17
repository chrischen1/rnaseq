final_data_path = '~/Dropbox/MU/workspace/ngs_tools/fast_rnaseq/test_data/'
output_path = '~/Dropbox/MU/workspace/ngs_tools/fast_rnaseq/results/'
  
cnt <- read.csv(paste0(final_data_path,'counts.csv'),row.names = 1,check.names = F)
grp <- read.csv(paste0(final_data_path,'group.csv'),row.names = 1,check.names = F,as.is = T)
rpkm_data <- read.csv(paste0(final_data_path,'rpkm.csv'),row.names = 1,check.names = F)
de_results_path <- paste0(output_path,'/DE/')
enrichment_results_path <- paste0(output_path,'/enrichment/')

dir.create(de_results_path,showWarnings = F)
dir.create(enrichment_results_path,showWarnings = F)


source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')

# Filtering
rpkm_data1 <- rpkm_data[apply(rpkm_data, 1, function(x)(sum(x>0)))!=0,]
cnt1 <- cnt[rownames(rpkm_data1),]

# heatmap
library(pheatmap)
if(length(unique(grp$group))==1){
  grp1 <- grp[,2,drop=F]
}else{
  grp1 <- grp[,1:2]
}
pheatmap(log2(rpkm_data1+0.01),show_rownames = F,scale = 'column', annotation_col = grp1)


# differential expression analysis
de_result <- edgeR_wrapper(cnt1,grp)

# PCA 
library(ggfortify)
library(pca3d)
pca <- prcomp(t(rpkm_data1), scale.=TRUE)
gr <- factor(grp$condition)
pca3d(pca, group=gr, show.ellipses=TRUE,ellipse.ci=0.75, show.plane=FALSE,show.group.labels = T)
autoplot(pca, data = data.frame(grp), colour = 'condition', frame = TRUE, frame.type = 'norm')

# analysis for normal vs normal_treatment

library(org.Mm.eg.db)
eg <- bitr(rownames(cnt), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
de_slt1 <- de_result$normal.normal_treatment
de_slt2 <- de_slt1[eg$ENSEMBL,]
de_slt2 <- cbind(de_slt2,'ENTREZID'=eg$ENTREZID)

# Volcano plot
plot_volcano_pval(logFC = de_slt1$logFC,P_Value = de_slt1$FDR,show_lines = F,left = 0,right = 0)

geneList_all <- de_slt2$logFC
names(geneList_all) <- de_slt2$ENTREZID
geneList <- geneList_all[de_slt1$FDR<0.05]

de_slt2 <- de_slt1[de_slt1$FDR<0.05,]
de_slt2 <- de_slt1[eg$ENSEMBL,]

# Pathway Enrichment Analysis
pathway_res <- enrichPathway(gene= names(geneList),organism = 'mouse',pvalueCutoff=0.05,qvalueCutoff = 1, readable=T)

geneList <- sort(geneList,decreasing = T)

enrichment_barplot(pathway_res, showCategory=10)
enrichment_dotplot(pathway_res, showCategory=10)
enrichment_emapplot(pathway_res, showCategory=30)
enrichment_cnetplot(pathway_res, categorySize="pvalue", foldChange=geneList)

viewPathway(pathway_res@result$Description[1], readable=TRUE, foldChange=geneList,organism = 'mouse')

# GSEA
gsea_res <- gsePathway(geneList, nPerm=10000,organism='mouse',pvalueCutoff=0.8,pAdjustMethod="BH", verbose=FALSE)
emapplot(gsea_res, color="pvalue")
gseaplot(gsea_res, geneSetID = gsea_res@result$ID[1])


# KEGG Pathway Enrichment Analysis
setwd(enrichment_results_path)
kegg_res <- enrichKEGG(names(geneList),organism = 'mmu',keyType = 'ncbi-geneid',pAdjustMethod='none',pvalueCutoff=0.05)
pathview(gene.data = geneList, pathway.id = "mmu05222", species = "mmu",
         limit = list(gene=5, cpd=1),kegg.dir = enrichment_results_path)

# GO enrichment analysis
library(topGO)
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Mm.eg.db", ID="entrez")
go_data <- new("topGOdata",ontology="BP",allGenes= geneList_all,annot=annFUN.GO2genes,
              GO2genes=allGO2genes,geneSel= selection,nodeSize=10)

results.ks <- runTest(go_data, algorithm="classic", statistic="ks")
goEnrichment <- GenTable(go_data, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
goEnrichment$KS <- as.numeric(goEnrichment$KS)

enrichment_goplot(goEnrichment)
par(cex = 0.25)
showSigOfNodes(go_data, score(results.ks), firstSigNodes = 3, useInfo = 'all')

# circus plot
go_res <- enrichGO(gene = names(geneList),OrgDb=org.Mm.eg.db,ont = "CC",pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
go_res2 <- go_res@result
go_res2 <- go_res2[order(go_res2$p.adjust),]
go_res3 <- enrichment_parser(go_res2,geneList_all)

eg2 <- bitr(rownames(go_res3), fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
go_res3 <- go_res3[eg2$ENTREZID,]
rownames(go_res3) <- eg2$SYMBOL # replace ENTREZID with gene SYMBOL
g_circus <- enrichment_circusplot(go_res3,space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)


# saving results
dir.create(de_results_path,showWarnings = F)
for(i in names(de_result)){
  write.csv(de_result[[i]],paste0(de_results_path,i,'.csv'))
}














