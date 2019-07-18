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



