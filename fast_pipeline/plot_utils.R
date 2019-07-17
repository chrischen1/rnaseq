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

enrichment_dotplot <- function(enrichment_result,x = 'GeneRatio',showCategory=10,color = "p.adjust",size = 'Count',
                               split = NULL,font.size=12,title = "",orderBy = "p.adjust",
                               colorBy = "p.adjust",decreasing=TRUE){
    df <- data.frame(enrichment_result@result)
    idx <- order(df[[orderBy]], decreasing = decreasing)[seq_len(showCategory)]
    df <- df[idx,]
    # df$Description <- factor(df$Description, levels=rev(unique(df$Description[idx])))
    ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
        geom_point() + scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
        ylab(NULL) + ggtitle(title) + theme_dose(font.size) + scale_size(range=c(3, 8))
}

enrichment_barplot <- function(enrichment_result, x="Count", colorBy='p.adjust',orderBy = "p.adjust", 
                               showCategory=8, font.size=12, title="", decreasing = TRUE,color = "p.adjust",
                               ...) {
    if (x == "geneRatio" || x == "GeneRatio") {
        x <- "GeneRatio"
    }else if (x == "count" || x == "Count") {
        x <- "Count"
    }
    
    df <- data.frame(enrichment_result@result)
    idx <- order(df[[orderBy]], decreasing = decreasing)[seq_len(showCategory)]
    df <- df[idx,]    
    if(colorBy %in% colnames(df)) {
        p <- ggplot(df, aes_string(x = "Description", y = x, fill = colorBy)) + 
            theme_dose(font.size) +
            scale_fill_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE))
    } else {
        p <- ggplot(df, aes_string(x = "Description", y = x, fill = "Description")) + 
            theme_dose(font.size) +
            theme(legend.position="none")
    }
    p + geom_bar(stat = "identity") + coord_flip() + ggtitle(title) 
}

enrichment_emapplot <- function(enrichment_result, showCategory = 30, color="p.adjust", layout = "kk", ...) {
    n <- showCategory
    
    geneSets <- setNames(strsplit(geneID(enrichment_result), "/", fixed=TRUE), rownames(enrichment_result@result))
    y <- as.data.frame(enrichment_result@result)
    if (is.numeric(n)) {
        y <- y[1:n,]
    } else {
        y <- y[match(n, y$Description),]
        n <- length(n)
    }
    
    
    if (n == 0) {
        stop("no enriched term found...")
    } else if (n == 1) {
        g <- graph.empty(0, directed=FALSE)
        g <- add_vertices(g, nv = 1)
        V(g)$name <- y$Description
        V(g)$color <- "red"
        return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
    } else {
        id <- y[,1]
        geneSets <- geneSets[id]
        
        n <- nrow(y) #
        w <- matrix(NA, nrow=n, ncol=n)
        colnames(w) <- rownames(w) <- y$Description
        
        for (i in 1:n) {
            for (j in i:n) {
                w[i,j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }
        
        wd <- melt(w)
        wd <- wd[wd[,1] != wd[,2],]
        wd <- wd[!is.na(wd[,3]),]
        g <- graph.data.frame(wd[,-3], directed=FALSE)
        E(g)$width=sqrt(wd[,3] * 5)
        g <- delete.edges(g, E(g)[wd[,3] < 0.2])
        idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))
        
        cnt <- sapply(geneSets[idx], length)
        V(g)$size <- cnt
        
        colVar <- y[idx, color]
        V(g)$color <- colVar
    }
    
    
    p <- ggraph(g, layout=layout)
    
    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
    }
    
    p + geom_node_point(aes_(color=~color, size=~size)) +
        geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
        scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
        ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
        scale_size(range=c(3, 8))
}

enrichment_cnetplot <- function(x,showCategory = 5,foldChange   = NULL,layout = "kk",
                                  colorEdge = FALSE,circular = FALSE,node_label = TRUE,...) {
    
    if (circular) {
        layout <- "linear"
        geom_edge <- geom_edge_arc
    } else {
        geom_edge <- geom_edge_link
    }
    
    geneSets <- extract_geneSets(x, showCategory)
    
    g <- list2graph(geneSets)
    
    foldChange <- fc_readable(x, foldChange)
    
    size <- sapply(geneSets, length)
    V(g)$size <- min(size)/2
    
    n <- length(geneSets)
    V(g)$size[1:n] <- size
    
    if (colorEdge) {
        E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
        edge_layer <- geom_edge(aes_(color = ~category), alpha=.8)
    } else {
        edge_layer <- geom_edge(alpha=.8, colour='darkgrey')
    }
    
    if (!is.null(foldChange)) {
        fc <- foldChange[V(g)$name[(n+1):length(V(g))]]
        V(g)$color <- NA
        V(g)$color[(n+1):length(V(g))] <- fc
        palette <- fc_palette(fc)
        p <- ggraph(g, layout=layout, circular = circular) +
            edge_layer +
            geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
            scale_color_gradientn(name = "fold change", colors=palette, na.value = "#E5C494")
    } else {
        V(g)$color <- "#B3B3B3"
        V(g)$color[1:n] <- "#E5C494"
        p <- ggraph(g, layout=layout, circular=circular) +
            edge_layer +
            geom_node_point(aes_(color=~I(color), size=~size))
    }
    
    p <- p + scale_size(range=c(3, 10), breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
        theme_void()
    
    if (node_label){
        p <- p + geom_node_text(aes_(label=~name), repel=TRUE)
    }
    return(p)
}

enrichment_goplot <- function(goEnrichment,title = ''){
    ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
        stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
        xlab("Biological process") +
        ylab("-log2 Enrichment KS-test pval") +
        ggtitle(title) +
        scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
        theme_bw(base_size=24) +
        theme(
            legend.position='none',
            legend.background=element_rect(),
            plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
            axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
            axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
            axis.title=element_text(size=24, face="bold"),
            legend.key=element_blank(),     #removes the border
            legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
            legend.text=element_text(size=18),  #Text size
            title=element_text(size=18)) +
        guides(colour=guide_legend(override.aes=list(size=2.5))) +
        coord_flip()
}

plot_volcano_pval <- function(logFC,P_Value,plot_title='',left=0.67,right=1.5,hide_black_dots =F,
                              show_lines=T,cex=2,col_cutoff=1,xlim=c(-6,6),ylim=c(0,20)){
    # Red indicates P_Value<0.05 and log2Fold_Change<-1, green is P_Value<0.05 and log2Fold_Change>1)
    # red indicates up regulated, green is downregulated
    
    if(hide_black_dots){
        plot(logFC, -log10(P_Value),col='white',
             cex = cex, main=plot_title, pch=20,xaxt="n", xlab= '', ylab="",cex.lab=1.8, xlim=xlim, ylim=ylim,cex.axis=1.8)
        
    }else{
        plot(logFC,-log10(P_Value),cex=cex,main=plot_title, pch=20, xaxt="n", xlab= '', ylab="",cex.lab=1.8, xlim=xlim, ylim=ylim,cex.axis=1.8)
    }
    axis(1, at = seq(-4,4, by = 0.5), las=2,cex.axis=1.8)
    title(xlab="log"["2"]~"(Fold change)", line=4, cex.lab=1.8)
    title(ylab=("-log"["10"]~"(FDR)"), line=2, cex.lab=1.8)
    #lines
    if(show_lines){
        abline(v=log2(right),lty=3,lwd=5,col='black')
        text(log2(right)*1.2,ylim[2]*0.9, paste("",right), col = "black", adj = c(0, -.1),cex=1.5)
        abline(v=log2(left),lty=3,lwd=5,col='black')
        text(log2(left)-0.9,ylim[2]*0.9, paste("",left), col = "black", adj = c(0, -.1),cex=1.5)
    }
    abline(v=log2(1),lty=3,lwd=5,col='black')
    
    
    up_ind <- P_Value<.05 & logFC > col_cutoff
    down_ind <- P_Value<.05 & logFC < (-col_cutoff)
    
    points(logFC[up_ind],  -log10(P_Value)[up_ind], pch=20, col="red",cex=cex)
    points(logFC[down_ind],  -log10(P_Value)[down_ind], pch=20, col="green",cex=cex)
}

enrichment_circusplot <-function (data, title, space, gene.order, gene.size, gene.space, 
                              nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col, border.size, 
                              process.label, limit) 
{
    y <- id <- xpro <- ypro <- xgen <- ygen <- lx <- ly <- ID <- logFC <- NULL
    Ncol <- dim(data)[2]
    if (missing(title)) 
        title <- ""
    if (missing(space)) 
        space = 0
    if (missing(gene.order)) 
        gene.order <- "none"
    if (missing(gene.size)) 
        gene.size <- 3
    if (missing(gene.space)) 
        gene.space <- 0.2
    if (missing(lfc.col)) 
        lfc.col <- c("brown1", "azure", "cornflowerblue")
    if (missing(lfc.min)) 
        lfc.min <- -3
    if (missing(lfc.max)) 
        lfc.max <- 3
    if (missing(border.size)) 
        border.size <- 0.5
    if (missing(process.label)) 
        process.label <- 11
    if (missing(limit)) 
        limit <- c(0, 0)
    if (gene.order == "logFC") 
        data <- data[order(data[, Ncol], decreasing = T), ]
    if (gene.order == "alphabetical") 
        data <- data[order(rownames(data)), ]
    if (sum(!is.na(match(colnames(data), "logFC"))) > 0) {
        if (nlfc == 1) {
            cdata <- check_chord(data[, 1:(Ncol - 1)], limit)
            lfc <- sapply(rownames(cdata), function(x) data[match(x, 
                                                                  rownames(data)), Ncol])
        }
        else {
            cdata <- check_chord(data[, 1:(Ncol - nlfc)], limit)
            lfc <- sapply(rownames(cdata), function(x) data[, 
                                                            (Ncol - nlfc + 1)])
        }
    }
    else {
        cdata <- check_chord(data, limit)
        lfc <- 0
    }
    if (missing(ribbon.col)) 
        colRib <- grDevices::rainbow(dim(cdata)[2])
    else colRib <- ribbon.col
    nrib <- colSums(cdata)
    ngen <- rowSums(cdata)
    Ncol <- dim(cdata)[2]
    Nrow <- dim(cdata)[1]
    colRibb <- c()
    for (b in 1:length(nrib)) colRibb <- c(colRibb, rep(colRib[b], 
                                                        202 * nrib[b]))
    r1 <- 1
    r2 <- r1 + 0.1
    xmax <- c()
    x <- 0
    for (r in 1:length(nrib)) {
        perc <- nrib[r]/sum(nrib)
        xmax <- c(xmax, (pi * perc) - space)
        if (length(x) <= Ncol - 1) 
            x <- c(x, x[r] + pi * perc)
    }
    xp <- c()
    yp <- c()
    l <- 50
    for (s in 1:Ncol) {
        xh <- seq(x[s], x[s] + xmax[s], length = l)
        xp <- c(xp, r1 * sin(x[s]), r1 * sin(xh), r1 * sin(x[s] + 
                                                               xmax[s]), r2 * sin(x[s] + xmax[s]), r2 * sin(rev(xh)), 
                r2 * sin(x[s]))
        yp <- c(yp, r1 * cos(x[s]), r1 * cos(xh), r1 * cos(x[s] + 
                                                               xmax[s]), r2 * cos(x[s] + xmax[s]), r2 * cos(rev(xh)), 
                r2 * cos(x[s]))
    }
    df_process <- data.frame(x = xp, y = yp, id = rep(c(1:Ncol), 
                                                      each = 4 + 2 * l))
    xp <- c()
    yp <- c()
    logs <- NULL
    x2 <- seq(0 - space, -pi - (-pi/Nrow) - space, length = Nrow)
    xmax2 <- rep(-pi/Nrow + space, length = Nrow)
    for (s in 1:Nrow) {
        xh <- seq(x2[s], x2[s] + xmax2[s], length = l)
        if (nlfc <= 1) {
            xp <- c(xp, (r1 + 0.05) * sin(x2[s]), (r1 + 0.05) * 
                        sin(xh), (r1 + 0.05) * sin(x2[s] + xmax2[s]), 
                    r2 * sin(x2[s] + xmax2[s]), r2 * sin(rev(xh)), 
                    r2 * sin(x2[s]))
            yp <- c(yp, (r1 + 0.05) * cos(x2[s]), (r1 + 0.05) * 
                        cos(xh), (r1 + 0.05) * cos(x2[s] + xmax2[s]), 
                    r2 * cos(x2[s] + xmax2[s]), r2 * cos(rev(xh)), 
                    r2 * cos(x2[s]))
        }
        else {
            tmp <- seq(r1, r2, length = nlfc + 1)
            for (t in 1:nlfc) {
                logs <- c(logs, data[s, (dim(data)[2] + 1 - t)])
                xp <- c(xp, (tmp[t]) * sin(x2[s]), (tmp[t]) * 
                            sin(xh), (tmp[t]) * sin(x2[s] + xmax2[s]), 
                        tmp[t + 1] * sin(x2[s] + xmax2[s]), tmp[t + 
                                                                    1] * sin(rev(xh)), tmp[t + 1] * sin(x2[s]))
                yp <- c(yp, (tmp[t]) * cos(x2[s]), (tmp[t]) * 
                            cos(xh), (tmp[t]) * cos(x2[s] + xmax2[s]), 
                        tmp[t + 1] * cos(x2[s] + xmax2[s]), tmp[t + 
                                                                    1] * cos(rev(xh)), tmp[t + 1] * cos(x2[s]))
            }
        }
    }
    if (lfc[1] != 0) {
        if (nlfc == 1) {
            df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), 
                                                            each = 4 + 2 * l), logFC = rep(lfc, each = 4 + 
                                                                                               2 * l))
        }
        else {
            df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:(nlfc * 
                                                                     Nrow)), each = 4 + 2 * l), logFC = rep(logs, 
                                                                                                            each = 4 + 2 * l))
        }
    }
    else {
        df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), 
                                                        each = 4 + 2 * l))
    }
    aseq <- seq(0, 180, length = length(x2))
    angle <- c()
    for (o in aseq) if ((o + 270) <= 360) 
        angle <- c(angle, o + 270)
    else angle <- c(angle, o - 90)
    df_texg <- data.frame(xgen = (r1 + gene.space) * sin(x2 + 
                                                             xmax2/2), ygen = (r1 + gene.space) * cos(x2 + xmax2/2), 
                          labels = rownames(cdata), angle = angle)
    df_texp <- data.frame(xpro = (r1 + 0.15) * sin(x + xmax/2), 
                          ypro = (r1 + 0.15) * cos(x + xmax/2), labels = colnames(cdata), 
                          stringsAsFactors = FALSE)
    cols <- rep(colRib, each = 4 + 2 * l)
    x.end <- c()
    y.end <- c()
    processID <- c()
    for (gs in 1:length(x2)) {
        val <- seq(x2[gs], x2[gs] + xmax2[gs], length = ngen[gs] + 
                       1)
        pros <- which((cdata[gs, ] != 0) == T)
        for (v in 1:(length(val) - 1)) {
            x.end <- c(x.end, sin(val[v]), sin(val[v + 1]))
            y.end <- c(y.end, cos(val[v]), cos(val[v + 1]))
            processID <- c(processID, rep(pros[v], 2))
        }
    }
    df_bezier <- data.frame(x.end = x.end, y.end = y.end, processID = processID)
    df_bezier <- df_bezier[order(df_bezier$processID, -df_bezier$y.end), 
                           ]
    x.start <- c()
    y.start <- c()
    for (rs in 1:length(x)) {
        val <- seq(x[rs], x[rs] + xmax[rs], length = nrib[rs] + 
                       1)
        for (v in 1:(length(val) - 1)) {
            x.start <- c(x.start, sin(val[v]), sin(val[v + 1]))
            y.start <- c(y.start, cos(val[v]), cos(val[v + 1]))
        }
    }
    df_bezier$x.start <- x.start
    df_bezier$y.start <- y.start
    df_path <- bezier(df_bezier, colRib)
    if (length(df_genes$logFC) != 0) {
        tmp <- sapply(df_genes$logFC, function(x) ifelse(x > 
                                                             lfc.max, lfc.max, x))
        logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, 
                                                lfc.min, x))
        df_genes$logFC <- logFC
    }
    theme_blank <- theme(axis.line = element_blank(), axis.text.x = element_blank(),
                         axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
                         axis.title.y = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
    
    g <- ggplot() + geom_polygon(data = df_process, aes(x, y, 
                                                        group = id), fill = "gray70", inherit.aes = F, color = "black") + 
        geom_polygon(data = df_process, aes(x, y, group = id), 
                     fill = cols, inherit.aes = F, alpha = 0.6, color = "black") + 
        geom_point(aes(x = xpro, y = ypro, size = factor(labels, 
                                                         levels = labels), shape = NA), data = df_texp) + 
        guides(size = guide_legend("Catgory", ncol = 4, byrow = T, 
                                   override.aes = list(shape = 22, fill = unique(cols), 
                                                       size = 8))) + theme(legend.text = element_text(size = process.label)) + 
        geom_text(aes(xgen, ygen, label = labels, angle = angle), 
                  data = df_texg, size = gene.size) + geom_polygon(aes(x = lx, 
                                                                       y = ly, group = ID), data = df_path, fill = colRibb, 
                                                                   color = "black", size = border.size, inherit.aes = F) + 
        labs(title = title) + theme_blank
    if (nlfc >= 1) {
        g + geom_polygon(data = df_genes, aes(x, y, group = id, 
                                              fill = logFC), inherit.aes = F, color = "black") + 
            scale_fill_gradient2("logFC", space = "Lab", low = lfc.col[3], 
                                 mid = lfc.col[2], high = lfc.col[1], guide = guide_colorbar(title.position = "top", 
                                                                                             title.hjust = 0.5), breaks = c(min(df_genes$logFC), 
                                                                                                                            max(df_genes$logFC)), labels = c(round(min(df_genes$logFC),2), 
                                                                                                                                                             round(max(df_genes$logFC),2))) + theme(legend.position = "bottom", 
                                                                                                                                                                                                    legend.background = element_rect(fill = "transparent"), 
                                                                                                                                                                                                    legend.box = "horizontal", legend.direction = "horizontal")
    }
    else {
        g + geom_polygon(data = df_genes, aes(x, y, group = id), 
                         fill = "gray50", inherit.aes = F, color = "black") + 
            theme(legend.position = "bottom", legend.background = element_rect(fill = "transparent"), 
                  legend.box = "horizontal", legend.direction = "horizontal")
    }
}

bezier <- function(data, process.col){
    x <- c()
    y <- c()
    Id <- c()
    sequ <- seq(0, 1, by = 0.01)
    N <- dim(data)[1]
    sN <- seq(1, N, by = 2)
    if (process.col[1] == '') col_rain <- grDevices::rainbow(N) else col_rain <- process.col
    for (n in sN){
        xval <- c(); xval2 <- c(); yval <- c(); yval2 <- c()
        for (t in sequ){
            xva <- (1 - t) * (1 - t) * data$x.start[n] + t * t * data$x.end[n]
            xval <- c(xval, xva)
            xva2 <- (1 - t) * (1 - t) * data$x.start[n + 1] + t * t * data$x.end[n + 1]
            xval2 <- c(xval2, xva2)
            yva <- (1 - t) * (1 - t) * data$y.start[n] + t * t * data$y.end[n]  
            yval <- c(yval, yva)
            yva2 <- (1 - t) * (1 - t) * data$y.start[n + 1] + t * t * data$y.end[n + 1]
            yval2 <- c(yval2, yva2)			
        }
        x <- c(x, xval, rev(xval2))
        y <- c(y, yval, rev(yval2))
        Id <- c(Id, rep(n, 2 * length(sequ)))
    }
    df <- data.frame(lx = x, ly = y, ID = Id)
    return(df)
}

check_chord <- function(mat, limit){
    
    if(all(colSums(mat) >= limit[2]) & all(rowSums(mat) >= limit[1])) return(mat)
    
    tmp <- mat[(rowSums(mat) >= limit[1]),]
    mat <- tmp[,(colSums(tmp) >= limit[2])]
    
    mat <- check_chord(mat, limit)
    return(mat)
}

rle_plot <- function(x,grp,main=''){
    library(ggplot2)
    library(reshape2)
    x <- log2(x+1)
    m <- x-apply(x,1,median)
    m1 <- melt(m)
    m1 <- cbind.data.frame(m1,'condition'=grp$condition[m1$Var2],stringsAsFactors=F)
    g <- ggplot(m1, aes(x=Var2, y=value, fill=condition)) +geom_boxplot()+ xlab(main)
    return(g)
}

ecdf_plot <- function(x){
    p_seq <- seq(0,1,0.2)
    cdf <- c()
    for(i in p_seq){
        cdf <- c(cdf,mean(x<=i))
    }
    lo <- loess(cdf~p_seq)
    plot(p_seq,cdf,xlim=c(0,1),ylim=c(0,1))
    lines(predict(lo), col='red', lwd=2)
}