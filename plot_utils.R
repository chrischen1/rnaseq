
#plot_distribution = function(ratio , title_name , image_name,left_line,right_line){
#  ratio[ratio>6]=6
#  library(ggplot2)
#  df = data.frame(ratio=ratio)
#  summary(ratio)
#  xlabel = c(sprintf("%.1f", round(seq(0,5.85,0.5),2)),'>6.0')
#  g1<-ggplot(df, aes(x=ratio)) + geom_histogram(binwidth=0.05, color="black", fill="cornflowerblue") +
#    geom_freqpoly(binwidth=0.05, size=1,col="gray24")+ theme_bw()  +
#    # scale_fill_gradient("Frequency",low = "green", high = "red") +
#    labs(title=NULL, x="Ratio", y ="Frequency")+
#    theme(plot.title = element_text(color="#666666", face="bold", size=35)) +
#    theme(legend.text = element_text( size=22),legend.title = element_text( size=22), legend.key.size = unit(1,"cm")) +
#    theme(axis.title = element_text(color="#666666", face="bold", size=32),axis.text=element_text(size=31))+
#    theme(panel.grid.major = element_line(colour = "black", linetype = "dotted"),panel.grid.minor.y = element_blank())+
#    geom_vline(xintercept = c(left_line,1,right_line),color = "black", size=1.5)+theme(plot.title = element_text(hjust = 0.5))
#    x_ranges = ggplot_build(g1)$layout$panel_ranges[[1]]$x.range
#    g1 = g1 + annotate("text", x = c(left_line-0.1,1.1,right_line+0.1), y = ggplot_build(g1)$layout$panel_ranges[[1]]$y.range[2]*0.94, label = sprintf("%.2f", round(c(left_line,1,right_line),2)) ,size=12)+
#      scale_x_continuous("Ratio",breaks = seq(0,6,0.5),  labels=xlabel,minor_breaks = seq(round(x_ranges[1],1),round(x_ranges[2],1),0.1))
#  #jpeg(image_name, width = 2080 , height = 700)
#  #plot(g1)
#  #dev.off()
#  return(g1)
#}
#

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




plot_distribution = function(ratio , title_name ,left_line,right_line,mid_shift=1.1){
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
  x <- log2(x)
  library(nortest)
  pval1 <- ad.test(x)
  pval2 <- cvm.test(x)
  pval3 <- lillie.test(x)

  res1 <- c(test_desc,pval1$p.value,pval2$p.value,pval3$p.value,mean(2^x),median(2^x))
  names(res1) <- c('Description',pval1$method,pval2$method,pval3$method,'mean','median')
  return(res1)
}



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


plot_volcano <- function(Fold_Change,Expre,P_Value,plot_title='',left=0.67,right=1.5,hide_black_dots =F,show_lines=T,cex=1,xlim=6,ylim=20){
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
  title(ylab=("log"["2"]~"(Avg. Expression)"), line=2, cex.lab=1.8)
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

plot_volcano_pval <- function(Fold_Change,P_Value,plot_title='',left=0.67,right=1.5,hide_black_dots =F,show_lines=T,cex=1,col_cutoff=1,xlim=6,ylim=20){
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

plot_vol_edgeR <- function(counts,grp,plot_title='',hide_black_dots =F,cex=1,show_lines=T,pval_plot=F){
  de1 <- edgeR_wrapper(counts[,rownames(grp)],grp)
  if(!pval_plot){
    plot_volcano(Fold_Change = 2^(de1[,3]),P_Value = de1[,1],show_lines = show_lines,hide_black_dots = hide_black_dots,cex=cex,Expre = apply(counts, 1, mean),plot_title=plot_title)
  }else{
    plot_volcano_pval(Fold_Change = 2^(de1[,3]),P_Value = de1[,1],show_lines = show_lines,hide_black_dots = hide_black_dots,cex=cex,plot_title=plot_title)
    
  }
}
                                     
