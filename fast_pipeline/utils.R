library(clusterProfiler)
library(pathview)
library(ReactomePA)
library(FGNet)
library(DOSE)
library(igraph)
library(ggraph)
library(reshape2)
library(dplyr)
library(GOplot)
library(topGO)
library(ggfortify)
library(pca3d)
library(pheatmap)
library(rgl)

color_palette <- function(colors) colorRampPalette(colors)(n = 299)

sig_palette <- color_palette(c("red", "yellow", "blue"))

heatmap_palette <- color_palette(c("red", "yellow", "green"))

overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}

fc_readable <- function(x, foldChange = NULL) {
    if (is.null(foldChange))
        return(NULL)
    
    if(x@readable) {
        gid <- names(foldChange)
        if (is(x, 'gseaResult')) {
            ii <- gid %in% names(x@geneList)
        } else {
            ii <- gid %in% x@gene
        }
        gid[ii] <- x@gene2Symbol[gid[ii]]
        names(foldChange) <- gid
    }
    return(foldChange)
}

fc_palette <- function(fc) {
    if (all(fc > 0, na.rm=TRUE)) {
        palette <- color_palette(c("blue", "red"))
    } else if (all(fc < 0, na.rm=TRUE)) {
        palette <- color_palette(c("green", "blue"))
    } else {
        palette <- color_palette(c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"))
    }
    return(palette)
}

update_n <- function(x, showCategory) {
    if (!is.numeric(showCategory)) {
        return(showCategory)
    }
    
    ## geneSets <- geneInCategory(x) ## use core gene for gsea result
    n <- showCategory
    if (nrow(x) < n) {
        n <- nrow(x)
    }
    
    return(n)
}

extract_geneSets <- function(x, n) {
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x@result)
    geneSets <- geneSets[y$ID]
    names(geneSets) <- y$Description
    if (is.numeric(n)) {
        return(geneSets[1:n])
    }
    return(geneSets[n]) ## if n is a vector of Description
}


list2graph <- function(inputList) {
    x <- list2df(inputList)
    g <- graph.data.frame(x, directed=FALSE)
    return(g)
}

list2df <- function(inputList) {
    ldf <- lapply(1:length(inputList), function(i) {
        data.frame(categoryID=rep(names(inputList[i]),
                                  length(inputList[[i]])),
                   Gene=inputList[[i]])
    })
    
    do.call('rbind', ldf)
}

enrichment_parser <- function(mat,logfc,idx = 1:10){
    mat <- mat[idx,]
    all_genes <- strsplit(mat$geneID,'/')
    all_genes_uniq <- unique(unlist(all_genes))
    res_mat <- matrix(0,nrow = length(all_genes_uniq),ncol = nrow(mat)+1)
    rownames(res_mat) <- all_genes_uniq
    colnames(res_mat) <- c(mat$Description,'logFC')
    for(i in 1:nrow(mat)){
        res_mat[all_genes[[i]],i] <- 1
    }
    res_mat[,ncol(res_mat)] <- logfc[rownames(res_mat)]
    return(res_mat)
}

selection <- function(allScore){ return(allScore < 0.05)}

go_enrich <- function(gene_all_etz_i,allGO2genes,annFUN.GO2genes,selection,ontology="BP"){
    go_data <- new("topGOdata",ontology=ontology,allGenes= gene_all_etz_i,annot=annFUN.GO2genes,
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
    return(list(goEnrichment,results.ks,go_data))
}