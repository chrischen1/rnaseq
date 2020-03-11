source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/fast_pipeline/utils.R')
source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/fast_pipeline/plot_utils.R')
source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/de_rnaseq.R')
library(org.Mm.eg.db)

final_data_path = '~/Dropbox/MU/workspace/will_rnaseq/final/'
output_path = '~/Documents/analysis_results/'

cnt <- read.csv(paste0(final_data_path,'counts.csv'),row.names = 1,check.names = F)
grp <- read.csv(paste0(final_data_path,'group.csv'),row.names = 1,check.names = F,as.is = T)
rpkm_data <- read.csv(paste0(final_data_path,'rpkm.csv'),row.names = 1,check.names = F)
de_results_path <- paste0(output_path,'/DE/')
comparisions_results_path <- paste0(output_path,'/comparisions/')

dir.create(output_path,showWarnings = F)
dir.create(de_results_path,showWarnings = F)
dir.create(comparisions_results_path,showWarnings = F)


# Filtering
rpkm_data1 <- rpkm_data[apply(rpkm_data, 1, function(x)(sum(x>0)))>22,]
cnt1 <- cnt[rownames(rpkm_data1),]

# differential expression analysis
de_result <- edgeR_wrapper(cnt1,grp)

for(i in names(de_result)){
    write.csv(de_result[[i]],paste0(de_results_path,i,'.csv'))
}

# heatmap
if(length(unique(grp$group))==1){
    grp1 <- grp[,2,drop=F]
}else{
    grp1 <- grp[,1:2]
}
tiff(paste0(output_path,'heatmap.tiff'))
pheatmap(log2(rpkm_data1+0.01),show_rownames = F,scale = 'column', annotation_col = grp1)
dev.off()

# PCA 
pca <- prcomp(t(rpkm_data1), scale.=TRUE)
gr <- factor(grp$condition)
pca3d(pca, group=gr, show.ellipses=TRUE,ellipse.ci=0.75, show.plane=FALSE,show.group.labels = T)
snapshotPCA3d(paste0(output_path,'PCA3D.png'))
rgl.close()

tiff(paste0(output_path,'PCA2D.tiff'))
autoplot(pca, data = data.frame(grp), colour = 'condition', frame = TRUE, frame.type = 'norm')
dev.off()

# analysis for each DE pair
eg <- bitr(rownames(cnt1), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

for(i in names(de_result)){
    out_dir_i <- paste0(comparisions_results_path,'/',i,'/')
    dir.create(out_dir_i,showWarnings = F)
    de_slt1 <- de_result[[i]]
    de_slt2 <- de_slt1[eg$ENSEMBL,]
    de_slt2 <- cbind.data.frame(de_slt2,'ENTREZID'=eg$ENTREZID,stringsAsFactors=F)
    
    # Volcano plot
    tiff(paste0(out_dir_i,'volcano_plot.tiff'))
    plot_volcano_pval(logFC = de_slt1$logFC,P_Value = de_slt1$FDR,show_lines = F,left = 0,right = 0)
    dev.off()
    
    geneList_all <- de_slt2$logFC
    names(geneList_all) <- de_slt2$ENTREZID
    
    geneList_dup <- geneList_all[de_slt2$FDR < max(sort(de_slt2$FDR)[100],0.05)]
    geneList <- c()
    for(geneList_name in names(geneList_dup)){
        geneList[geneList_name] <- mean(geneList_dup[names(geneList_dup)==geneList_name])
    }
    geneList <- sort(geneList,decreasing = T)
    
    # Pathway Enrichment Analysis
    pathway_res <- enrichPathway(gene = names(geneList),organism = 'mouse',pvalueCutoff=0.05,qvalueCutoff = 0.05, readable=T)
    write.csv(pathway_res@result,paste0(out_dir_i,'pathwayEnrichment.csv'))
    
    tiff(paste0(out_dir_i,'enrichment_barplot.tiff'),width = 1920,height = 1080)
    enrichment_barplot(pathway_res, showCategory=20,colorBy = 'p.adjust',orderBy = 'p.adjust',color='p.adjust')
    dev.off()
    
    tiff(paste0(out_dir_i,'enrichment_dotplot.tiff'),width = 1920,height = 1080)
    enrichment_dotplot(pathway_res, showCategory=20,colorBy = 'p.adjust',orderBy = 'p.adjust',color='p.adjust')
    dev.off()
    
    tiff(paste0(out_dir_i,'enrichment_emapplot.tiff'),width = 1920,height = 1080)
    enrichment_emapplot(pathway_res, showCategory=30)
    dev.off()
    
    tiff(paste0(out_dir_i,'enrichment_cnetplot.tiff'),width = 1920,height = 1080)
    enrichment_cnetplot(pathway_res, categorySize="pvalue", foldChange=geneList)
    dev.off()
    
    pathway_network_path_i <- paste0(out_dir_i,'/pathway_network/')
    dir.create(pathway_network_path_i,showWarnings = F)
    for (pathway_idx in seq_len(min(10,nrow(pathway_res@result)))) {
        tryCatch(expr = {
            g <- viewPathway(pathway_res@result$Description[pathway_idx], readable=TRUE, foldChange=geneList,organism = 'mouse')
        },error = function(x){print(paste0(pathway_res@result$Description[pathway_idx],' not found.'))})
        tiff(paste0(pathway_network_path_i,pathway_res@result$ID[pathway_idx],'.tiff'),width = 1920,height = 1080)
        plot(g)
        dev.off()
    }
    
    # GSEA
    gsea_res <- gsePathway(geneList, nPerm=10000,organism='mouse',pvalueCutoff=0.8,pAdjustMethod="BH", verbose=FALSE)
    write.csv(gsea_res@result,paste0(out_dir_i,'GSEAEnrichment.csv'))
    # emapplot(gsea_res, color="pvalue")
    gsea_path_i <- paste0(out_dir_i,'/GSEA/')
    dir.create(gsea_path_i,showWarnings = F)
    
    for (pathway_idx in seq_len(min(10,nrow(gsea_res@result)))) {
        g <- gseaplot(gsea_res, geneSetID = gsea_res@result$ID[pathway_idx],title = pathway_res@result$Description[pathway_idx])
        tiff(paste0(gsea_path_i,pathway_res@result$ID[pathway_idx],'.tiff'),width = 1920,height = 1080)
        plot(g)
        dev.off()
    }
    
    
    # KEGG Pathway Enrichment Analysis
    kegg_path_i <- paste0(out_dir_i,'/KEGG/')
    dir.create(kegg_path_i,showWarnings = F)
    setwd(kegg_path_i)
    kegg_res <- enrichKEGG(names(geneList),organism = 'mmu',keyType = 'ncbi-geneid',pAdjustMethod='BH',pvalueCutoff=0.05)
    write.csv(kegg_res@result,paste0(out_dir_i,'KEGGEnrichment.csv'))
    
    for (pathway_idx in seq_len(min(10,nrow(kegg_res@result)))) {
        tiff(paste0(kegg_path_i,pathway_res@result$ID[pathway_idx],'.tiff'),width = 1920,height = 1080)
        tryCatch(expr = {
            pathview(gene.data = geneList, pathway.id = kegg_res@result$ID[pathway_idx], species = "mmu",
                     limit = list(gene=5, cpd=1),kegg.dir = kegg_path_i)
        },error = function(x){print(paste0(kegg_res@result$ID[pathway_idx],' not found.'))})
        dev.off()
    }
    
    # GO enrichment analysis
    allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Mm.eg.db", ID="entrez")
    go_data <- new("topGOdata",ontology="BP",allGenes= geneList_all,annot=annFUN.GO2genes,
                   GO2genes=allGO2genes,geneSel= selection,nodeSize=10)
    
    results.ks <- runTest(go_data, algorithm="classic", statistic="ks")
    goEnrichment <- GenTable(go_data, KS=results.ks, orderBy="KS", topNodes=20)
    goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
    goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
    goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
    goEnrichment$KS <- as.numeric(goEnrichment$KS)
    goEnrichment$KS[is.na(goEnrichment$KS)] <- 1e-30
    goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
    write.csv(goEnrichment,paste0(out_dir_i,'goEnrichment.csv'))
    
    go_path_i <- paste0(out_dir_i,'/GO/')
    dir.create(go_path_i,showWarnings = F)
    g <- enrichment_goplot(goEnrichment)
    
    tiff(paste0(go_path_i,'enrichment_goplot.tiff'),width = 1920,height = 1080)
    plot(g)
    dev.off()
    
    tiff(paste0(go_path_i,'go_showSigOfNodes.tiff'),width = 1920,height = 1080)
    par(cex = 0.25)
    showSigOfNodes(go_data, score(results.ks), firstSigNodes = 3, useInfo = 'all')
    dev.off()
    
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
    tiff(paste0(go_path_i,'circus_goplot.tiff'),width = 1920,height = 1080)
    plot(g_circus)
    dev.off()
}

