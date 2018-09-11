
modify_ylab <- function(y_labs){
  for(i in 1:length(y_labs)){
    y <- as.numeric(y_labs[i])
    if(length(grep('.',y_labs[i],fixed = T))>0){
      y_labs[i] <- as.character(round(y))
    }
    if(y < 10){
      y_labs[i] <- paste(y_labs[i],'\t    ',sep = '')
    }else if(y < 100){
      y_labs[i] <- paste(y_labs[i],'\t   ',sep = '')
    }else if(y < 1000){
      y_labs[i] <- paste(y_labs[i],'\t  ',sep = '')
    }
  }
  return(y_labs)
}

normality_check = function(treatment_info) {
  control_name = treatment_info$control_name
  treatment_name = treatment_info$treatment_name
  cis_genes = treatment_info$cis_genes
  split2cis = treatment_info$split2cis
  control_rowMeans = treatment_info$control_counts
  treatment_rowMeans = treatment_info$treatment_counts
  GeneID = treatment_info$GeneID
  outputfile = treatment_info$outputfile
  left_line = treatment_info$left_line
  right_line = treatment_info$right_line
  data = get_all_summary_matrix(GeneID,control_name,treatment_name,control_rowMeans,treatment_rowMeans)
  ratio_summary_matrix = data$ratio_summary_matrix
  filtered_ratio_summary_matrix = get_ratio_summary_matrix(ratio_summary_matrix) 
  ratio = (filtered_ratio_summary_matrix[,2])/(filtered_ratio_summary_matrix[,1])  
  if(split2cis!=1){
    res1 <- get_normal_pvalue(ratio,outputfile)
    return(res1)
  }else{
    cis_ratio_summary = filtered_ratio_summary_matrix[rownames(filtered_ratio_summary_matrix) %in% as.character(unlist(cis_genes)),]
    cis_ratio = (cis_ratio_summary[,2])/(cis_ratio_summary[,1])
    res1 <- get_normal_pvalue(cis_ratio,paste(outputfile,'cis',sep = '_'))
    trans_ratio_summary = filtered_ratio_summary_matrix[!(rownames(filtered_ratio_summary_matrix) %in% as.character(unlist(cis_genes))),]
    trans_ratio = (trans_ratio_summary[,2])/(trans_ratio_summary[,1])  
    res1 <- rbind(res1,get_normal_pvalue(trans_ratio,paste(outputfile,'trans',sep = '_')))
    return(res1)
  }
}

plot_distribution_split <- function(cis_ratio,trans_ratio,left_line,right_line,title_name=''){
  df <- data.frame('ratio'=c(cis_ratio,trans_ratio),'type'=c(rep('cis',length(cis_ratio)),rep('trans',length(trans_ratio))))
  xlabel = c(sprintf("%.1f", round(seq(0,5.85,0.5),2)),'>6.0')
  g1 <- ggplot(df, aes(x=ratio, color=type,fill=type)) + geom_histogram(binwidth=0.05,alpha=0.3, position="identity")+
    theme(legend.position="bottom")+geom_vline(xintercept = c(left_line,1,right_line),color = "black", size=0.5)+
    labs(title=title_name, x="Ratio", y ="Frequency")+
    theme(plot.title = element_text(color="#666666", face="bold", size=35)) +
    theme(legend.text = element_text( size=22),legend.title = element_text( size=22), legend.key.size = unit(1,"cm")) +
    theme(axis.title = element_text(color="#666666", face="bold", size=32),axis.text=element_text(size=35))+
    theme(panel.grid.major = element_line(colour = "black", linetype = "dotted"),panel.grid.minor.y = element_blank())
  x_ranges = ggplot_build(g1)$layout$panel_ranges[[1]]$x.range
  y_labs = ggplot_build(g1)$layout$panel_ranges[[1]]$y.labels
  y_labs = modify_ylab(y_labs)
  g1 <- g1 + annotate("text", x = c(left_line-0.14,1.1,right_line+0.14), y = ggplot_build(g1)$layout$panel_ranges[[1]]$y.range[2]*0.98, label = sprintf("%.2f",round(c(left_line,1,right_line),2)),size=10)+
    scale_x_continuous("Ratio",limits = c(-0.3,6.3),breaks = seq(0,6,0.5),  labels=xlabel,minor_breaks = seq(-0.3,6,0.1)) +
    scale_y_continuous("Frequency",breaks = as.numeric(y_labs),labels = y_labs)
  
  return(g1)
}

plot_distribution = function(ratio , title_name ,left_line,right_line,max_ratio = 6){
  ratio[ratio>max_ratio]=max_ratio
  library(ggplot2)
  df = data.frame(ratio=ratio)
  summary(ratio)
  xlabel = c(sprintf("%.1f", round(seq(0,5.85,0.5),2)),'>6.0')
  g1<-ggplot(df, aes(x=ratio)) + geom_histogram(binwidth=0.05, color="black", fill="cornflowerblue") +
    geom_freqpoly(binwidth=0.05, size=1,col="gray24")+ theme_bw()  +
    # scale_fill_gradient("Frequency",low = "green", high = "red") +
    labs(title=title_name, x="Ratio", y ="Frequency")+
    theme(plot.title = element_text(color="#666666", face="bold", size=35)) +
    theme(legend.text = element_text( size=22),legend.title = element_text( size=22), legend.key.size = unit(1,"cm")) +
    theme(axis.title = element_text(color="#666666", face="bold", size=32),axis.text=element_text(size=35))+
    theme(panel.grid.major = element_line(colour = "black", linetype = "dotted"),panel.grid.minor.y = element_blank())+
    geom_vline(xintercept = c(left_line,1,right_line),color = "black", size=1.5)+theme(plot.title = element_text(hjust = 0.5))
  return(g1)
}
  
plot_distribution2 = function(ratio , title_name ,left_line,right_line,mid_shift=1.1){
  ratio[ratio>6]=6
  library(ggplot2)
  df = data.frame(ratio=ratio)
  summary(ratio)
  xlabel = c(sprintf("%.1f", round(seq(0,5.85,0.5),2)),'>6.0')
  g1<-ggplot(df, aes(x=ratio)) + geom_histogram(binwidth=0.05, color="black", fill="cornflowerblue") +
    geom_freqpoly(binwidth=0.05, size=1,col="gray24")+ theme_bw()  +
    # scale_fill_gradient("Frequency",low = "green", high = "red") +
    labs(title=title_name, x="Ratio", y ="Frequency")+
    theme(plot.title = element_text(color="#666666", face="bold", size=35)) +
    theme(legend.text = element_text( size=22),legend.title = element_text( size=22), legend.key.size = unit(1,"cm")) +
    theme(axis.title = element_text(color="#666666", face="bold", size=32),axis.text=element_text(size=35))+
    theme(panel.grid.major = element_line(colour = "black", linetype = "dotted"),panel.grid.minor.y = element_blank())+
    geom_vline(xintercept = c(left_line,1,right_line),color = "black", size=1.5)+theme(plot.title = element_text(hjust = 0.5))
  x_ranges = ggplot_build(g1)$layout$panel_ranges[[1]]$x.range
  y_labs = ggplot_build(g1)$layout$panel_ranges[[1]]$y.labels
  y_labs = modify_ylab(y_labs)
  g1 = g1 + annotate("text", x = c(left_line-0.14,mid_shift,right_line+0.14), y = ggplot_build(g1)$layout$panel_ranges[[1]]$y.range[2]*0.98, label = sprintf("%.2f",round(c(left_line,1,right_line),2)),size=10)+
    scale_x_continuous("Ratio",limits = c(-0.3,6.3),breaks = seq(0,6,0.5),  labels=xlabel,minor_breaks = seq(-0.3,6,0.1)) +
    scale_y_continuous("Frequency",breaks = as.numeric(y_labs),labels = y_labs)
  return(g1)
}

# Normalality check
get_normal_pvalue <- function(x,test_desc){
  test_desc <- gsub('.+/','',test_desc)
  test_desc <- gsub('.jpeg','',test_desc)
  test_desc <- gsub('.jpg','',test_desc)
  
  if(length(x) <7)return(c(test_desc,rep('no_meaningful_pval',3)))
  x[x>6] <- 6
  x[x<1/6] <- 1/6
  log2x <- log2(x)
  library(nortest)
  pval1 <- ad.test(log2x)
  pval2 <- cvm.test(log2x)
  pval3 <- lillie.test(log2x)

  res1 <- c(test_desc,pval1$p.value,pval2$p.value,pval3$p.value,mean(x),median(x),sd(x))
  names(res1) <- c('Description',pval1$method,pval2$method,pval3$method,'mean','median','sd')
  return(res1)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,layout.pos.col = matchidx$col))
    }
  }
}

plot_distribution_pairs <- function(r,cis_genes,title_name = '',left_line = 0.67,right_line = 1.5,max_ratio = 6){
  g <- multiplot(plotlist = list(plot_distribution(r[cis_genes],title_name = 'Cis genes',left_line = left_line,right_line = right_line,max_ratio = max_ratio),
                                 plot_distribution(r[!names(r)%in%cis_genes],title_name = 'Trans genes',left_line = left_line,right_line = right_line,max_ratio = max_ratio)))
  return(g)
}

plot_volcano <- function(Fold_Change,Expre,P_Value,plot_title='',left=0.67,right=1.5,hide_black_dots =F,show_lines=T,cex=2,xlim=6,ylim=20){
  # Red indicates P_Value<0.05 and log2Fold_Change<-1, green is P_Value<0.05 and log2Fold_Change>1)
  # red indicates up regulated, green is downregulated
  
  if(hide_black_dots){
    plot(log2(Fold_Change), log2(Expre),col='white',
         cex = cex, main=plot_title, pch=20,xaxt="n", xlab= '', ylab="",cex.lab=1.8, xlim=c(-xlim,xlim), ylim=c(0,xlim),cex.axis=1.8)
    
  }else{
    plot(log2(Fold_Change), log2(Expre),cex=cex,main=plot_title, pch=20, xaxt="n", xlab= '', ylab="",cex.lab=1.8, xlim=c(-xlim,xlim), ylim=c(0,20),cex.axis=1.8)
  }
  axis(1, at = seq(-4,4, by = 0.5), las=2,cex.axis=1.8)
  title(xlab="log"["2"]~"(Ratio)", line=4, cex.lab=1.8)
  title(ylab=("log"["2"]~"(S Expression)"), line=2, cex.lab=1.8)
  #lines
  if(show_lines){
    abline(v=log2(1.5),lty=3,lwd=5,col='black')
    text(log2(right),0, paste("",right), col = "black", adj = c(0, -.1),cex=1.5)
    abline(v=log2(2/3),lty=3,lwd=5,col='black')
    text(log2(left),0, paste("",left), col = "black", adj = c(0, -.1),cex=1.5)
  }
  abline(v=log2(1),lty=3,lwd=5,col='black')
  

  up_ind <- P_Value<.05 & log2(Fold_Change) > 1
  down_ind <- P_Value<.05 & log2(Fold_Change) < -1
  
  points(log2(Fold_Change)[up_ind], log2(Expre)[up_ind], pch=20, col="red",cex=cex)
  points(log2(Fold_Change)[down_ind], log2(Expre)[down_ind], pch=20, col="green",cex=cex)
}

plot_volcano_pairs <- function(Fold_Change,Expre,P_Value,cis_genes,plot_title='',left=0.67,right=1.5,hide_black_dots =F,show_lines=T,cex=1,xlim=6,ylim=20){
  trans_genes <- names(Fold_Change)[!names(Fold_Change)%in%cis_genes]
  plot_volcano(Fold_Change = Fold_Change[cis_genes],plot_title='cis',Expre = sum_expre[cis_genes],P_Value = P_Value[cis_genes])
  plot_volcano(Fold_Change = Fold_Change[trans_genes],plot_title='trans',Expre = sum_expre[trans_genes],P_Value = P_Value[trans_genes])
}
  
plot_volcano_pval <- function(Fold_Change,P_Value,plot_title='',left=0.67,right=1.5,hide_black_dots =F,show_lines=T,cex=2,col_cutoff=1,xlim=6,ylim=20){
  # Red indicates P_Value<0.05 and log2Fold_Change<-1, green is P_Value<0.05 and log2Fold_Change>1)
  # red indicates up regulated, green is downregulated
  
  if(hide_black_dots){
    plot(log2(Fold_Change), -log10(P_Value),col='white',
         cex = cex, main=plot_title, pch=20,xaxt="n", xlab= '', ylab="",cex.lab=1.8, xlim=c(-xlim,xlim), ylim=c(0,ylim),cex.axis=1.8)
    
  }else{
    plot(log2(Fold_Change),-log10(P_Value),cex=cex,main=plot_title, pch=20, xaxt="n", xlab= '', ylab="",cex.lab=1.8, xlim=c(-xlim,xlim), ylim=c(0,ylim),cex.axis=1.8)
  }
  axis(1, at = seq(-4,4, by = 0.5), las=2,cex.axis=1.8)
  title(xlab="log"["2"]~"(Ratio)", line=4, cex.lab=1.8)
  title(ylab=("-log"["10"]~"(p-value)"), line=2, cex.lab=1.8)
  #lines
  if(show_lines){
    abline(v=log2(1.5),lty=3,lwd=5,col='black')
    text(log2(right),0, paste("",right), col = "black", adj = c(0, -.1),cex=1.5)
    abline(v=log2(2/3),lty=3,lwd=5,col='black')
    text(log2(left),0, paste("",left), col = "black", adj = c(0, -.1),cex=1.5)
  }
  abline(v=log2(1),lty=3,lwd=5,col='black')
  
  
  up_ind <- P_Value<.05 & log2(Fold_Change) > col_cutoff
  down_ind <- P_Value<.05 & log2(Fold_Change) < (-col_cutoff)
  
  points(log2(Fold_Change)[up_ind],  -log10(P_Value)[up_ind], pch=20, col="red",cex=cex)
  points(log2(Fold_Change)[down_ind],  -log10(P_Value)[down_ind], pch=20, col="green",cex=cex)
}

plot_cis_trans <- function(x,gene_list,chr,image_name){
  x1 <- x[x$gene_id %in% gene_list,]
  cis_ratio <- (x1$treatment/x1$control)[x1$chr_id==as.character(chr)]
  trans_ratio <- (x1$treatment/x1$control)[x1$chr_id!=as.character(chr)]
  cis_ratio[cis_ratio>6] <- 6
  trans_ratio[trans_ratio>6] <- 6
  g_cis <- plot_distribution(cis_ratio,title_name = '',left_line = 0.67,right_line = 1.5)
  g_trans <- plot_distribution(trans_ratio,title_name = '',left_line = 0.67,right_line = 1.5)
  jpeg(image_name, width = 40, height = 28, units = 'in', res = 300)
  multiplot(g_cis,g_trans)
  dev.off()
}

plot_ploid <- function(x,gene_list,image_name){
  x1 <- x[x$gene_id %in% gene_list,]
  ratio1 <- (x1$treatment/x1$control)
  ratio1[ratio1>6] <- 6
  g <- plot_distribution(ratio1,title_name = '',left_line = 0.67,right_line = 1.5)
  jpeg(image_name, width = 40, height = 28, units = 'in', res = 300)
  plot(g)
  dev.off()
}

get_m_gene_list <- function(type='gbM'){
  m <- read.csv('~/Arabidopsis_project_finalization/data/FigureSet9/status_mat.csv',as.is = T,row.names = 1)
  m2 <- m[,c(4:8,3,2)]
  m3 <- apply(m2,2,function(x)(x==m$diploid& x == type))
  gene_list <- list()
  for(i in 1:7){
    gene_list[[i]] <- rownames(m)[m3[,i]]
  }
  return(gene_list)
}


get_group <- function(x){
  grp = rep('trt',length(x))
  grp[grep('cont',tolower(x))] <- 'ctr'
  return(as.factor(grp))
}

get_norm_counts <- function(genes_data,spikes_data){
  library(RUVSeq)
  
  genes_data_f <- genes_data[apply(genes_data, 1, function(x) length(x[x>=10])>=3 ),]
  spikes_data_f <- spikes_data[apply(spikes_data, 1, function(x) length(x[x>=10])>=3 ),]
  
  all_data_f <- rbind(genes_data_f,spikes_data_f[,colnames(genes_data_f)])
  grp <- get_group(colnames(all_data_f))
  set = newSeqExpressionSet(as.matrix(all_data_f),phenoData = data.frame(grp,row.names = colnames(all_data_f)))
  
  set = betweenLaneNormalization(set,which = 'upper')
  set1 = RUVg(set, cIdx = rownames(spikes_data_f) , k = 1)
  counts_normalized <- normCounts(set1)
  counts_normalized_gene = counts_normalized[rownames(genes_data_f),]
  return(counts_normalized_gene)
}

plot_chromesome <- function(ratio,chr){
  library(ggplot2)
  seq <- 1:length(ratio)
  m <- cbind.data.frame(seq,ratio,chr)
  colnames(m) <- c('seq','ratio','chr')
  m$chr <- factor(chr,levels(factor(chr))[c(1,3:10,2)])
  g1 <- ggplot(m, aes(x=seq, y=ratio,color=chr)) +
    geom_point(size = 0.8)+ theme_bw()+ facet_grid(~chr,scales="free_x")+ylim(c(-4,4))+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    labs(x="", y="Ratio (log2)", title="") +  theme(legend.position="none") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  return(g1)
}

plot_vol_edgeR <- function(counts,grp,plot_title='',hide_black_dots =F,cex=1,show_lines=T,pval_plot=F,gene_keep=NULL){
  de1 <- data.frame(edgeR_wrapper(counts[,rownames(grp)],grp)[[1]])
  if(!is.null(gene_keep)){
    de1 <- de1[gene_keep,]
    counts <- counts[gene_keep,]
  }
  if(!pval_plot){
    plot_volcano(Fold_Change = 2^(de1$logFC),P_Value = de1$FDR,show_lines = show_lines,hide_black_dots = hide_black_dots,cex=cex,Expre = apply(counts, 1, mean),plot_title=plot_title)
  }else{
    plot_volcano_pval(Fold_Change = 2^(de1$logFC),P_Value = de1$FDR,show_lines = show_lines,hide_black_dots = hide_black_dots,cex=cex,plot_title=plot_title)
    
  }
}
       
plot_pca <- function(i,j,df_pca,col,pch=1){
  if(i==j){
    plot(df_pca$x[,i], df_pca$x[,j],col='white',xlab = '',ylab = '',xaxt='n',yaxt='n')
    legend('center',legend = paste('PC',i),cex = 3,bty = 'n')
  }else{
    per_sdv <- round(df_pca$sdev/sum(df_pca$sdev),4)*100
    plot(df_pca$x[,j], df_pca$x[,i],col = col,pch = as.numeric(pch),xlab = paste('PC',j,' ',per_sdv[j],'%',sep = ''),ylab = paste('PC',i,' ',per_sdv[i],'%',sep = ''))

  }
}

plot_pca_2pc <- function(exp,grp){
  library(ggplot2)
  if(sum(colnames(exp)%in%rownames(grp))!=nrow(grp)){
    warning('row names of group table is different with colnames of expression matrix')
    return(0)
  }
  exp <- exp[,rownames(grp)]
  pca_res  <- prcomp(t(exp))
  per_sdv <- round(pca_res$sdev/sum(pca_res$sdev),4)*100
  df_pca <- data.frame('PC1'=pca_res$x[,1],'PC2'=pca_res$x[,2],'Condition'=grp$condition)
  g <- ggplot(data = df_pca,aes(x=PC1, y=PC2,color=Condition))+geom_point(size=3)+labs(x = paste('PC1 ',per_sdv[1],'%',sep = ''),y = paste('PC2 ',per_sdv[2],'%',sep = ''))
  plot(g)
}

plot_pca_3pc <- function(exp,grp){
  grp <- grp[colnames(exp),]
  df_pca <- prcomp(t(exp))
  par(mfrow=c(3,3))
  for(i in 1:3){
    for(j in 1:3){
      plot_pca(i,j,df_pca,col = factor(grp[,2]),pch=factor(grp[,1]))
    }
  }
}

# param: df: a dataframe with 3 cloumns, V1 is name of enriched catagory, V2 is -log10 p adjusted and V3 is count of genes in the pathway
barplot_enrichment <- function(df){
  library(ggplot2)
  g <- ggplot(data=df, aes(x=V1,y=V2)) +geom_bar(stat="identity",fill="cornflowerblue")+ xlab("") + ylab("-Log10 adjusted pvalue")+theme(axis.text=element_text(size=22,face="bold"),axis.title=element_text(size=24,face="bold"))+
    geom_text(aes(label=V3), hjust=1.0,vjust=0, color="white",position = position_dodge(1), size=9.5) +coord_flip()
  return(g)
}


enrichment_circus <-function (data, title, space, gene.order, gene.size, gene.space, 
                              nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col, border.size, 
                              process.label, limit) 
{
  library(GOplot)
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

enrichment_parser <- function(mat,logfc,idx = 1:10){
  mat <- mat[idx,]
  all_genes <- strsplit(mat$SYMBOL,'/')
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



