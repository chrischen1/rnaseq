get_all_summary_matrix = function(GeneId,control_name,treatment_name,control_rowMean, treatment_rowMean){
  #get gene id
  gene_id = GeneId
  head(gene_id)
  
  #get control average
  control_avg = control_rowMean
  
  #get treatment average
  treatment_avg = treatment_rowMean
  #combine ratio into frame
  ratio_summary = data.frame(controlavg=control_avg,treatmentavg=treatment_avg)
  row.names(ratio_summary) <- gene_id
  head(ratio_summary)
  ratio_summary_matrix <- data.matrix(ratio_summary)
  head(ratio_summary_matrix)
  
  data = list("ratio_summary_matrix" = ratio_summary_matrix)
  return(data)
  
}


get_ratio_summary_matrix = function(ratio_summary_matrix ){
  Chr_select = rowSums(ratio_summary_matrix) != 0
  Chr_ratio_summary_matrix =ratio_summary_matrix[Chr_select,]
  Chr_ratio_summary_matrix[Chr_ratio_summary_matrix==0]=0.0000000001
  return(Chr_ratio_summary_matrix)
}


#
#
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

# junk functions (used as black box)
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

get_result_by_treatment_name = function(treatment_info) {
  
  ##########################  set information 
  control_name = treatment_info$control_name
  treatment_name = treatment_info$treatment_name
  cis_genes = treatment_info$cis_genes
  split2cis = treatment_info$split2cis
  
  
  #import control data
  control_rowMeans = treatment_info$control_counts
  head(control_rowMeans)
  
  #import treatment data
  treatment_rowMeans = treatment_info$treatment_counts
  head(treatment_rowMeans)
  
  #import geneid
  GeneID = treatment_info$GeneID
  head(GeneID)
  
  outputfile = treatment_info$outputfile
  
  left_line = treatment_info$left_line
  right_line = treatment_info$right_line
  
  
  # check if gene number is the same
  try(if(length(control_rowMeans) != length(treatment_rowMeans)) stop("The gene id number is not equal"))
  
  data = get_all_summary_matrix(GeneID,control_name,treatment_name,control_rowMeans,treatment_rowMeans)
  ratio_summary_matrix = data$ratio_summary_matrix
  
  filtered_ratio_summary_matrix = get_ratio_summary_matrix(ratio_summary_matrix) 
  
  head(filtered_ratio_summary_matrix)
  
  ratio = (filtered_ratio_summary_matrix[,2])/(filtered_ratio_summary_matrix[,1])  
  if(split2cis!=1)
  {
      # title_name = paste('Histogram for ratio of ',treatment_name,'/',control_name,' for all genes',sep="")
      title_name = NULL
      ggplot1 <- plot_distribution(ratio,title_name ,left_line,right_line)
      image_name = outputfile
      jpeg(image_name, width = 32, height = 16, units = 'in', res = 300)
      plot(ggplot1)
      dev.off()      
  }else{
  
    ## set cis genes 
    cis_ratio_summary = filtered_ratio_summary_matrix[rownames(filtered_ratio_summary_matrix) %in% as.character(unlist(cis_genes)),]
    cis_ratio = (cis_ratio_summary[,2])/(cis_ratio_summary[,1])  
    
    # title_name = paste('Histogram for ratio of ',treatment_name,'/',control_name,' for all cis genes',sep="")
    title_name = NULL
    ggplot2 <- plot_distribution(cis_ratio,title_name,left_line,right_line)
    
    
    ## set trans genes 
    trans_ratio_summary = filtered_ratio_summary_matrix[!(rownames(filtered_ratio_summary_matrix) %in% as.character(unlist(cis_genes))),]
    trans_ratio = (trans_ratio_summary[,2])/(trans_ratio_summary[,1])  
    
    # title_name = paste('Histogram for ratio of ',treatment_name,'/',control_name,' for all trans genes',sep="")
    title_name = NULL
    ggplot3 <- plot_distribution(trans_ratio,title_name,left_line,right_line)
    
    image_name = outputfile
    jpeg(image_name, width = 40, height = 28, units = 'in', res = 300)
    multiplot(ggplot2,ggplot3)
    dev.off()
  }
}

# for most figs
plot_setup<- function(data_dir,results_dir,normal_check=F){
  files = list.files(data_dir)
  normal_check_res <- NULL
  for (file in files){
    l = line_position(file)
    if(length(grep('Trisomy',file))>0){
      cis_id <- gsub('.*Trisomy_(\\d+).+','\\1',file)
    }else{
      cis_id = 0
    }
    prefix = gsub('.txt','',file)
    output_figure = paste(results_dir,prefix,'.jpeg',sep = '')
    datafile=paste(data_dir,file,sep = '')
    
    dataset = read.table(datafile,header=T,sep="\t",as.is = T)
    chr = dataset[,1]
    geneid = dataset[,2]
    control_exp = dataset[,4]
    treatment_exp = dataset[,3]
    cis_genes= geneid[chr==cis_id]
    split2cis=0
    if(cis_id!=0)
    {
      split2cis=1
    }
    treatment_info = list("control_name" = "control", "control_counts" = control_exp, "treatment_name" = "treatment", "treatment_counts" = treatment_exp,"GeneID" = geneid,"cis_genes" = cis_genes,"split2cis"=split2cis,"outputfile"=output_figure,"left_line" = l$leftline,"right_line" = l$rightline)
    if(normal_check){
      normal_check_res <- rbind(normal_check_res,normality_check(treatment_info))
    }else{
      get_result_by_treatment_name(treatment_info)
    }
  }
  if(!is.null(normal_check_res)){
    return(normal_check_res)
  }
}


# for TF
plot_tf_trisomy <- function(tf_ratio,target_table,image_name,chr,line_pos,normality_check=F){
  target_ratio_cis <- (target_table$treatment/target_table$control)[target_table$chr_id==chr]
  target_ratio_trans <- (target_table$treatment/target_table$control)[target_table$chr_id!=chr]
  if(normality_check){
    return(rbind(get_normal_pvalue(tf_ratio,paste(image_name,'tf_ratio',sep = '_')),
                 get_normal_pvalue(target_ratio_cis,paste(image_name,'target_ratio_cis',sep = '_')),
                 get_normal_pvalue(target_ratio_trans,paste(image_name,'target_ratio_trans',sep = '_'))))
  }else{
    ggplot1 <- plot_distribution(tf_ratio,NULL,line_pos$leftline,line_pos$rightline)
    ggplot2 <- plot_distribution(target_ratio_cis,NULL,line_pos$leftline,line_pos$rightline)
    ggplot3 <- plot_distribution(target_ratio_trans,NULL,line_pos$leftline,line_pos$rightline)
    jpeg(image_name, width = 40, height = 28, units = 'in', res = 300)
    multiplot(ggplot1,ggplot2,ggplot3)
    dev.off()
  }
}

plot_tf_ploid <- function(tf_table,target_table,image_name,line_pos,normality_check=F){
  tf_ratio <- tf_table$treatment/tf_table$control
  target_ratio <- target_table$treatment/target_table$control
  if(normality_check){
    return(rbind(get_normal_pvalue(tf_ratio,paste(image_name,'tf_ratio',sep = '_')),
                 get_normal_pvalue(target_ratio,paste(image_name,'target_ratio',sep = '_'))))
  }else{
    ggplot1 <- plot_distribution(tf_ratio,NULL,line_pos$leftline,line_pos$rightline)
    ggplot2 <- plot_distribution(target_ratio,NULL,line_pos$leftline,line_pos$rightline)
    jpeg(image_name, width = 40, height = 28, units = 'in', res = 300)
    multiplot(ggplot1,ggplot2)
    dev.off()
  }
}

plot_tf <- function(data_dir,results_dir,normality_check=F){
  normality_check_res <- NULL
  #Trisomy
  for(i in 1:5){
    tf_table_name <- paste(data_dir,'Transcription_Factor_Trisomy_',i,'.txt',sep = '')
    target_cis_table_name <- paste(data_dir,'Trisomy_',i,'_Cis_Target.txt',sep = '')
    target_trans_table_name <- paste(data_dir,'Trisomy_',i,'_Trans_Target.txt',sep = '')
    tf_table <- read.table(tf_table_name,header=T,sep="\t",as.is = T)
    cis_table <- read.table(target_cis_table_name,header=T,sep="\t",as.is = T)
    trans_table <- read.table(target_trans_table_name,header=T,sep="\t",as.is = T)
    l <- line_position('Trisomy')
    tf_cis_ratio <- (tf_table$treatment/tf_table$control)[tf_table$chr_id==i]
    tf_trans_ratio <- (tf_table$treatment/tf_table$control)[tf_table$chr_id!=i]
    #cis_tf
    image_name1 = paste(results_dir,'Trisomy_',i,'_Cis_Transcription_Factors.jpg',sep = '')
    #trans_tf
    image_name2 = paste(results_dir,'Trisomy_',i,'_Trans_Transcription_Factors.jpg',sep = '')
    if(normality_check){
      normality_check_res <- rbind(normality_check_res,plot_tf_trisomy(tf_ratio = tf_cis_ratio,target_table = cis_table,image_name = image_name1,chr=i,line_pos=l,normality_check = T))
      normality_check_res <- rbind(normality_check_res,plot_tf_trisomy(tf_ratio = tf_trans_ratio,target_table = trans_table,image_name = image_name2,chr=i,line_pos=l,normality_check = T))
      
    }else{
      plot_tf_trisomy(tf_ratio = tf_cis_ratio,target_table = cis_table,image_name = image_name,chr=i,line_pos=l)
      plot_tf_trisomy(tf_ratio = tf_trans_ratio,target_table = trans_table,image_name = image_name,chr=i,line_pos=l)
    }
  }
  
  #Triploid
  tf_table1 <- read.table(paste(data_dir,'Transcription_Factor_Triploid.txt',sep=''),header=T,sep="\t",as.is = T)
  target_table1 <- read.table(paste(data_dir,'Target_Gene_Triploid.txt',sep=''),header=T,sep="\t",as.is = T)
  image_name1 = paste(results_dir,'Triploid_Transcription_Factors.jpg',sep = '')
  
  #Tetraploid
  tf_table2 <- read.table(paste(data_dir,'Transcription_Factor_Tetraploid.txt',sep=''),header=T,sep="\t",as.is = T)
  target_table2 <- read.table(paste(data_dir,'Target_Gene_Tetraploid.txt',sep=''),header=T,sep="\t",as.is = T)
  image_name2 = paste(results_dir,'Tetraploid_Transcription_Factors.jpg',sep = '')
  
  
  if(normality_check){
    normality_check_res <- rbind(normality_check_res,plot_tf_ploid(tf_table1,target_table1,image_name1,line_pos=l,normality_check = T))
    normality_check_res <- rbind(normality_check_res,plot_tf_ploid(tf_table2,target_table2,image_name2,line_pos=line_position('Tetraploid'),normality_check = T))
  }else{
    plot_tf_ploid(tf_table,target_table,image_name,line_pos=l)
    plot_tf_ploid(tf_table,target_table,image_name,line_pos=line_position('Tetraploid'))
  }
  if(!is.null(normality_check_res))return(normality_check_res)
}

# rearranged functions

normality_check_from_list <- function(m,gene_list,desc){
  
  m2 <- m[gene_list,]
  return(get_normal_pvalue(m2$treatment/m2$control,test_desc = desc))
}

normality_check_setup <- function(gene_list,gene_list_desc,file_path='~/Arabidopsis_project_finalization/data/FigureSet1/'){
  gene_list2 <- gene_list
  if(!is.list(gene_list)){
    gene_list <- list()
    for(i in 1:7){
      gene_list[[i]] <- gene_list2
    }
  }
  
  pval_mat <- NULL
  for(i in 1:5){
    ratio_data <- read.delim(paste(file_path,'Trisomy_',i,'.txt',sep = ''),as.is = T)
    ratio_data_cis <- ratio_data[ratio_data$gene_id %in% gene_list[[i]] & ratio_data$chr_id==i,]
    ratio_data_trans <- ratio_data[ratio_data$gene_id %in% gene_list[[i]] & ratio_data$chr_id!=i,]
    cis_res <- get_normal_pvalue(ratio_data_cis$treatment/ratio_data_cis$control,test_desc = paste(gene_list_desc,'_cis_Trisomy_',i,sep = ''))
    trans_res <- get_normal_pvalue(ratio_data_trans$treatment/ratio_data_trans$control,test_desc = paste(gene_list_desc,'_trans_Trisomy_',i,sep = ''))
    pval_mat <- rbind(pval_mat,cis_res,trans_res)
  }
  m1 <- read.delim(paste(file_path,'Triploid.txt',sep = ''),as.is = T)
  m2 <- read.delim(paste(file_path,'Tetraploid.txt',sep = ''),as.is = T)
  m1 <- m1[m1$gene_id %in% gene_list[[6]],]
  m2 <- m2[m2$gene_id %in% gene_list[[7]],]
  
  m1_res <- get_normal_pvalue(m1$treatment/m1$control,test_desc = paste(gene_list_desc,'_Triploid',sep = ''))
  m2_res <- get_normal_pvalue(m2$treatment/m2$control,test_desc = paste(gene_list_desc,'_Tetraploid',sep = ''))
  pval_mat <- rbind(pval_mat,m1_res,m2_res)
  return(pval_mat)
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

line_position = function(file_name){
  file_name <- tolower(file_name)
  if(length(grep('trisomy',file_name))!=0 | length(grep('triploid',file_name))!=0| length(grep('duplicate_all_cis_tf_trans_target',file_name))!=0| length(grep('single_all_cis_tf_trans_target',file_name))!=0){
    l <- list('leftline'=0.67,'rightline'=1.5)
  }else{
    l <- list('leftline'=0.5,'rightline'=2)
  }
  return(l)
}

diff_ratio_test <- function(dataset,id,xlab_title,plot_title,figurename,gene_list=''){
  if(gene_list!=''){
    dataset <- dataset[dataset$Gene.id %in% gene_list,]
  }
  # Make volcano plot
  Fold_Change = dataset$Avg_treatment
  Expre = dataset$Sum_expression
  P_Value = dataset$p.value
  
  
  cis_select = grepl(id,dataset$Gene.id)
  dataset_cis=dataset[cis_select,]
  Fold_Change = dataset_cis$Avg_treatment
  Expre = dataset_cis$Sum_expression
  P_Value = dataset_cis$p.value
  jpeg(figurename, width = 16, height = 21, units = 'in', res = 300)
  par(mfrow=c(2,1))
  plot(log2(Fold_Change), log2(dataset_cis$Sum_expression), pch=20, xaxt="n", xlab= '', ylab="",cex.lab=1.8, xlim=c(-3.5,3.5), ylim=c(0,20),main=plot_title,cex.axis=1.8)
  axis(1, at = seq(-4,4, by = 0.5), las=2,cex.axis=1.8)
  title(xlab=xlab_title, line=4, cex.lab=1.8)
  title(ylab=("log"["2"]~"(Avg. Expression)"), line=2, cex.lab=1.8)
  #head(log2(Fold_Change))
  #head(log2(Fold_Change))
  #lines
  abline(v=log2(1.5),lty=3,lwd=5,col='black')
  text(0.7,0, "Ratio = 1.5", col = "black", adj = c(0, -.1),cex=1.5)
  abline(v=log2(2/3),lty=3,lwd=5,col='black')
  text(-1.8,0, "Ratio = 0.67", col = "black", adj = c(0, -.1),cex=1.5)
  abline(v=log2(1),lty=3,lwd=5,col='black')
  
  # Red indicates P_Value<0.05 and log2Fold_Change<-1.3, green is P_Value<0.05 and log2Fold_Change>1.3)
  # red indicates up regulated, green is downregulated
  # adjust the numbers, colors, or variables as you deem fit
  
  #with(subset(res, P_Value<.05 & Fold_Change>1.3), points(log2(Fold_Change), -log10(P_Value), pch=20, col="red"))
  with(subset(dataset_cis, P_Value<.05 & Fold_Change>1), points(log2(Avg_treatment), log2(Sum_expression), pch=20, col="red"))
  with(subset(dataset_cis, P_Value<.05 & Fold_Change<1), points(log2(Avg_treatment), log2(Sum_expression), pch=20, col="green"))
  
  
  #################  trans
  dataset_trans=dataset[!cis_select,]
  Fold_Change = dataset_trans$Avg_treatment
  Expre = dataset_trans$Sum_expression
  P_Value = dataset_trans$p.value
  with(dataset_trans, plot(log2(Fold_Change), log2(Sum_expression), pch=20, xaxt="n", xlab='', ylab="",cex.lab=1.8,cex.axis=1.8, xlim=c(-3.5,3.5), ylim=c(0,20),main=plot_title))
  axis(1, at = seq(-4,4, by = 0.5), las=2,cex.axis=1.8)
  title(xlab=xlab_title, line=4, cex.lab=1.8)
  title(ylab=("log"["2"]~"(Avg. Expression)"), line=2, cex.lab=1.8)
  #head(log2(Fold_Change))
  #head(log2(Fold_Change))
  #lines
  abline(v=log2(1.5),lty=3,lwd=5,col='black')
  text(0.7,0, "Ratio = 1.5", col = "black", adj = c(0, -.1),cex=1.5)
  abline(v=log2(2/3),lty=3,lwd=5,col='black')
  text(-1.8,0, "Ratio = 0.67", col = "black", adj = c(0, -.1),cex=1.5)
  abline(v=log2(1),lty=3,lwd=5,col='black')
  
  # Red indicates P_Value<0.05 and log2Fold_Change<-1.3, green is P_Value<0.05 and log2Fold_Change>1.3)
  # red indicates up regulated, green is downregulated
  # adjust the numbers, colors, or variables as you deem fit
  
  #with(subset(res, P_Value<.05 & Fold_Change>1.3), points(log2(Fold_Change), -log10(P_Value), pch=20, col="red"))
  with(subset(dataset_trans, P_Value<.05 & Fold_Change>1), points(log2(Avg_treatment), log2(Sum_expression), pch=20, col="red"))
  
  with(subset(dataset_trans, P_Value<.05 & Fold_Change<1), points(log2(Avg_treatment), log2(Sum_expression), pch=20, col="green"))
  
  
  dev.off()
}

diff_ratio_test_polid <- function(dataset,xlab_title,plot_title,figurename,gene_list=''){
  if(gene_list!=''){
    dataset <- dataset[dataset$Gene.id %in% gene_list,]
  }
  # Make volcano plot
  Fold_Change = dataset$Avg_treatment
  Expre = dataset$Sum_expression
  P_Value = dataset$p.value
  
  
  jpeg(figurename, width = 16, height = 11, units = 'in', res = 300)
  par(mfrow=c(1,1))
  
  
  with(dataset, plot(log2(Fold_Change), log2(Sum_expression), pch=20, xaxt="n", xlab='', ylab="",cex.lab=1.8,cex.axis=1.8, xlim=c(-3.5,3.5), ylim=c(0,20)))
  axis(1, at = seq(-4,4, by = 0.5), las=2,cex.axis=1.8)
  title(xlab=xlab_title, line=4, cex.lab=1.8)
  title(ylab=("log"["2"]~"(Avg. Expression)"), line=2, cex.lab=1.8)
  #head(log2(Fold_Change))
  #head(log2(Fold_Change))
  #lines
  if(length(grep('Tetraploid',figurename))>0){
    abline(v=1,lty=3,lwd=5,col='black')
    text(0.7,0, "Ratio = 2.0", col = "black", adj = c(0, -.1),cex=1.5)
    abline(v=-1,lty=3,lwd=5,col='black')
    text(-1.8,0, "Ratio = 0.5", col = "black", adj = c(0, -.1),cex=1.5)
    abline(v=0,lty=3,lwd=5,col='black')
  }else{
    abline(v=log2(1.5),lty=3,lwd=5,col='black')
    text(0.7,0, "Ratio = 1.5", col = "black", adj = c(0, -.1),cex=1.5)
    abline(v=log2(2/3),lty=3,lwd=5,col='black')
    text(-1.8,0, "Ratio = 0.67", col = "black", adj = c(0, -.1),cex=1.5)
    abline(v=log2(1),lty=3,lwd=5,col='black')
  }

  
  # Red indicates P_Value<0.05 and log2Fold_Change<-1.3, green is P_Value<0.05 and log2Fold_Change>1.3)
  # red indicates up regulated, green is downregulated
  # adjust the numbers, colors, or variables as you deem fit
  
  #with(subset(res, P_Value<.05 & Fold_Change>1.3), points(log2(Fold_Change), -log10(P_Value), pch=20, col="red"))
  with(subset(dataset, P_Value<.05 & Fold_Change>1), points(log2(Avg_treatment), log2(Sum_expression), pch=20, col="red"))
  
  with(subset(dataset, P_Value<.05 & Fold_Change<1), points(log2(Avg_treatment), log2(Sum_expression), pch=20, col="green"))
  
  
  dev.off()
}

plot_vol <- function(out_path,gene_list=NULL,gene_list_name='',logratio_path = '~/Arabidopsis_project_finalization/data/FigureSet2/Trisomy_for_t_test/'){
  if(is.null(gene_list)){
    gene_list1 <- list()
    for(i in 1:7){
      gene_list1[[i]] <- read.table(paste(logratio_path,"Trisomy1_ratio_for_t_test_exp_logratio.diff",sep=''), header=TRUE,sep="\t",as.is = T)$Gene.id
    }
  }else if(!is.list(gene_list)){
    gene_list1 <- list()
    for(i in 1:7){
      gene_list1[[i]] <- gene_list
    }
  }else{
    gene_list1 <- gene_list
  }
  
  ## trisomy1
  dataset <- read.table(paste(logratio_path,"Trisomy1_ratio_for_t_test_exp_logratio.diff",sep=''), header=TRUE,sep="\t",as.is = T)
  xlab_title = "log"["2"]~"(Trisomy1 Ratio)"
  plot_title = NULL
  figurename = paste(out_path,gene_list_name,'Diff_ratio_trisomy1_diploid_log_cis_trans.jpeg',sep = '')
  diff_ratio_test(dataset[dataset$Gene.id %in% gene_list1[[1]],],'AT1G',xlab_title,plot_title,figurename)
  
  
  ## trisomy2
  dataset <- read.table(paste(logratio_path,"Trisomy2_ratio_for_t_test_exp_logratio.diff",sep=''), header=TRUE,sep="\t",as.is = T)
  xlab_title = "log"["2"]~"(Trisomy2 Ratio)"
  plot_title = NULL
  figurename = paste(out_path,gene_list_name,'Diff_ratio_trisomy2_diploid_log_cis_trans.jpeg',sep = '')
  diff_ratio_test(dataset[dataset$Gene.id %in% gene_list1[[2]],],'AT2G',xlab_title,plot_title,figurename)
  
  ## trisomy3
  dataset <- read.table(paste(logratio_path,"Trisomy3_ratio_for_t_test_exp_logratio.diff",sep=''), header=TRUE,sep="\t",as.is = T)
  xlab_title = "log"["2"]~"(Trisomy3 Ratio)"
  plot_title = NULL
  figurename = paste(out_path,gene_list_name,'Diff_ratio_trisomy3_diploid_log_cis_trans.jpeg',sep = '')
  diff_ratio_test(dataset[dataset$Gene.id %in% gene_list1[[3]],],'AT3G',xlab_title,plot_title,figurename)
  
  ## trisomy4
  dataset <- read.table(paste(logratio_path,"Trisomy4_ratio_for_t_test_exp_logratio.diff",sep=''), header=TRUE,sep="\t",as.is = T)
  xlab_title = "log"["2"]~"(Trisomy4 Ratio)"
  plot_title = NULL
  figurename = paste(out_path,gene_list_name, 'Diff_ratio_trisomy4_diploid_log_cis_trans.jpeg',sep = '')
  diff_ratio_test(dataset[dataset$Gene.id %in% gene_list1[[4]],],'AT4G',xlab_title,plot_title,figurename)
  
  ## trisomy5
  dataset <- read.table(paste(logratio_path,"Trisomy5_ratio_for_t_test_exp_logratio.diff",sep=''), header=TRUE,sep="\t",as.is = T)
  xlab_title = "log"["2"]~"(Trisomy5 Ratio)"
  plot_title = NULL
  figurename = paste(out_path,gene_list_name,'Diff_ratio_trisomy5_diploid_log_cis_trans.jpeg',sep = '')
  diff_ratio_test(dataset[dataset$Gene.id %in% gene_list1[[5]],],'AT5G',xlab_title,plot_title,figurename)
  
  ## triploid
  dataset <- read.table(paste(logratio_path,"Triploid_ratio_for_t_test_exp_logratio.diff",sep=''), header=TRUE,sep="\t",as.is = T)
  xlab_title = "log"["2"]~"(Triploid Ratio)"
  plot_title = NULL
  figurename = paste(out_path,gene_list_name, 'Diff_ratio_Triploid_diploid_log.jpeg',sep = '')
  diff_ratio_test_polid(dataset[dataset$Gene.id %in% gene_list1[[6]],],xlab_title,plot_title,figurename)
  
  ## tetraploid
  dataset <- read.table(paste(logratio_path,"Tetraploid_ratio_for_t_test_exp_logratio.diff",sep=''), header=TRUE,sep="\t",as.is = T)
  xlab_title = "log"["2"]~"(Tetraploid Ratio)"
  plot_title =NULL
  figurename = paste(out_path,gene_list_name, 'Diff_ratio_Tetraploid_diploid_log.jpeg',sep = '')
  diff_ratio_test_polid(dataset[dataset$Gene.id %in% gene_list1[[7]],],xlab_title,plot_title,figurename)
}

plot_vol_single <- function(dataset,plot_title,gene_list=NULL,left=0.67,right=1.5){
  if(!is.null(gene_list)){
    dataset <- dataset[dataset$Gene.id %in% gene_list,]
  }
  # Make volcano plot
  Fold_Change = dataset$Avg_treatment
  Expre = dataset$Sum_expression
  P_Value = dataset$p.value
  
  plot(log2(Fold_Change), log2(Expre), pch=20, xaxt="n", xlab= '', ylab="",cex.lab=1.8, xlim=c(-3.5,3.5), ylim=c(0,20),main=plot_title,cex.axis=1.8)
  axis(1, at = seq(-4,4, by = 0.5), las=2,cex.axis=1.8)
  title(xlab="log"["2"]~"(Ratio)", line=4, cex.lab=1.8)
  title(ylab=("log"["2"]~"(Avg. Expression)"), line=2, cex.lab=1.8)
  #lines
  abline(v=log2(1.5),lty=3,lwd=5,col='black')
  text(log2(right),0, paste("Ratio = ",right), col = "black", adj = c(0, -.1),cex=1.5)
  abline(v=log2(2/3),lty=3,lwd=5,col='black')
  text(log2(left),0, paste("Ratio = ",left), col = "black", adj = c(0, -.1),cex=1.5)
  abline(v=log2(1),lty=3,lwd=5,col='black')
  
  # Red indicates P_Value<0.05 and log2Fold_Change<-1.3, green is P_Value<0.05 and log2Fold_Change>1.3)
  # red indicates up regulated, green is downregulated
  # adjust the numbers, colors, or variables as you deem fit
  
  #with(subset(res, P_Value<.05 & Fold_Change>1.3), points(log2(Fold_Change), -log10(P_Value), pch=20, col="red"))
  with(subset(dataset, P_Value<.05 & Fold_Change>1), points(log2(Avg_treatment), log2(Sum_expression), pch=20, col="red"))
  with(subset(dataset, P_Value<.05 & Fold_Change<1), points(log2(Avg_treatment), log2(Sum_expression), pch=20, col="green"))
}

plot_vol_tf <- function(data_dir,results_dir){
  dataset <- read.table("Trisomy1_ratio_for_t_test_exp_logratio.diff", header=TRUE,sep="\t")
  #Trisomy
  for(i in 1:5){
    tf_table_name <- paste(data_dir,'Transcription_Factor_Trisomy_',i,'.txt',sep = '')
    target_cis_table_name <- paste(data_dir,'Trisomy_',i,'_Cis_Target.txt',sep = '')
    target_trans_table_name <- paste(data_dir,'Trisomy_',i,'_Trans_Target.txt',sep = '')
    tf_table <- read.table(tf_table_name,header=T,sep="\t",as.is = T)
    cis_table <- read.table(target_cis_table_name,header=T,sep="\t",as.is = T)
    trans_table <- read.table(target_trans_table_name,header=T,sep="\t",as.is = T)
    #cis_tf
    gene_list <- tf_table$gene_id[tf_table$chr_id==i]
    image_name = paste(results_dir,'vol_Trisomy_',i,'_Cis_Transcription_Factors.jpg',sep = '')
    target_genes_cis <- cis_table$gene_id[cis_table$chr_id==i]
    target_genes_trans <- cis_table$gene_id[cis_table$chr_id!=i]
    jpeg(image_name, width = 16, height = 21, units = 'in', res = 300)
    par(mfrow=c(3,1))
    plot_vol_single(dataset,'',gene_list)
    plot_vol_single(dataset,'',target_genes_cis)
    plot_vol_single(dataset,'',target_genes_trans)
    dev.off()
    
    #trans_tf
    gene_list <- tf_table$gene_id[tf_table$chr_id!=i]
    image_name = paste(results_dir,'vol_Trisomy_',i,'_Trans_Transcription_Factors.jpg',sep = '')
    target_genes_cis <- trans_table$gene_id[trans_table$chr_id==i]
    target_genes_trans <- trans_table$gene_id[trans_table$chr_id!=i]
    jpeg(image_name, width = 16, height = 21, units = 'in', res = 300)
    par(mfrow=c(3,1))
    plot_vol_single(dataset,'',gene_list)
    plot_vol_single(dataset,'',target_genes_cis)
    plot_vol_single(dataset,'',target_genes_trans)
    dev.off()
    
  }
  
  #Triploid
  tf_table <- read.table(paste(data_dir,'Transcription_Factor_Triploid.txt',sep=''),header=T,sep="\t",as.is = T)
  target_table <- read.table(paste(data_dir,'Target_Gene_Triploid.txt',sep=''),header=T,sep="\t",as.is = T)
  tf_genes <- tf_table$gene_id
  target_genes <- target_table$gene_id
  image_name = paste(results_dir,'vol_Triploid_Transcription_Factors.jpg',sep = '')
  jpeg(image_name, width = 16, height = 21, units = 'in', res = 300)
  par(mfrow=c(2,1))
  plot_vol_single(dataset,'',tf_genes)
  plot_vol_single(dataset,'',target_genes)
  dev.off()
  
  #Tetraploid
  tf_table <- read.table(paste(data_dir,'Transcription_Factor_Tetraploid.txt',sep=''),header=T,sep="\t",as.is = T)
  target_table <- read.table(paste(data_dir,'Target_Gene_Tetraploid.txt',sep=''),header=T,sep="\t",as.is = T)
  tf_genes <- tf_table$gene_id
  target_genes <- target_table$gene_id
  image_name = paste(results_dir,'vol_Tetraploid_Transcription_Factors.jpg',sep = '')
  jpeg(image_name, width = 16, height = 21, units = 'in', res = 300)
  par(mfrow=c(2,1))
  plot_vol_single(dataset,'',tf_genes,left = 0.5,right = 2)
  plot_vol_single(dataset,'',target_genes,left = 0.5,right = 2)
  dev.off()
}

plot_fig4 <- function(control_file,treatment_file,figurename,left_line=0.67,right_line=1.5,control_rep_num = 3,treatment_rep_num = 2,normality_check=F){
  control = read.table(control_file,header=T,sep="\t",as.is = T)
  treatment = read.table(treatment_file,header=T,sep="\t",as.is = T)
  gene_id = control$Gene.ID
  control_avg = control$Average
  treatment_avg = treatment$Average
  ratio_summary = data.frame(controlavg=control_avg,treatmentavg=treatment_avg)
  row.names(ratio_summary) <- gene_id
  ratio_summary_matrix <- data.matrix(ratio_summary)
  control_rep = control[,seq(2,2+control_rep_num-1)]
  treatment_rep = treatment[,seq(2,2+treatment_rep_num-1)]
  rep_summary = data.frame(controlrep=control_rep,treatmentrep=treatment_rep)
  row.names(rep_summary) <- gene_id
  rep_summary_matrix <- data.matrix(rep_summary)
  Chr_id = gsub('.+/(.+).txt','\\1',treatment_file)
  search = paste('AT',gsub('\\D','',Chr_id),'G',sep = '')
  Chr_select = grepl(search,row.names(rep_summary_matrix))
  none_Chr_select = !Chr_select
  Chr_rep_summary_matrix =rep_summary_matrix[Chr_select,] 
  none_Chr_rep_summary_matrix =rep_summary_matrix[none_Chr_select,] 
  Chr_ratio_summary_matrix =ratio_summary_matrix[Chr_select,]
  none_Chr_ratio_summary_matrix =ratio_summary_matrix[none_Chr_select,]
  Chr_select2 = rowSums(Chr_ratio_summary_matrix) != 0
  none_Chr_select2 = rowSums(none_Chr_ratio_summary_matrix) != 0
  Chr_rep_summary_matrix =Chr_rep_summary_matrix[Chr_select2,] 
  none_Chr_rep_summary_matrix =none_Chr_rep_summary_matrix[none_Chr_select2,] 
  Chr_ratio_summary_matrix =Chr_ratio_summary_matrix[Chr_select2,]
  none_Chr_ratio_summary_matrix =none_Chr_ratio_summary_matrix[none_Chr_select2,]
  Chr_ratio_summary_matrix[Chr_ratio_summary_matrix==0]=0.0000000001
  none_Chr_ratio_summary_matrix[none_Chr_ratio_summary_matrix==0]=0.0000000001
  ratio = (Chr_ratio_summary_matrix[,2])/(Chr_ratio_summary_matrix[,1])
  none_ratio = (none_Chr_ratio_summary_matrix[,2])/(none_Chr_ratio_summary_matrix[,1])
  low_expressed = none_Chr_ratio_summary_matrix[none_Chr_ratio_summary_matrix[,1]<1,]
  ratio1 = (low_expressed[,2])/(low_expressed[,1])
  ratio1[ratio1>6]=6 # set large number to same value
  # find moderate expressed genes in deploid, RPKM > 1 and < 100
  moderate_expressed = none_Chr_ratio_summary_matrix[none_Chr_ratio_summary_matrix[,1]>1 & none_Chr_ratio_summary_matrix[,1]<100,]
  ratio2 = (moderate_expressed[,2])/(moderate_expressed[,1])
  ratio2[ratio2>6]=6 # set large number to same value
  # find high expressed genes in deploid, RPKM>100
  high_expressed = none_Chr_ratio_summary_matrix[none_Chr_ratio_summary_matrix[,1]>100,]
  ratio3 = (high_expressed[,2])/(high_expressed[,1])
  ratio3[ratio3>6]=6 # set large number to same value
  if(normality_check){
    p1 <- get_normal_pvalue(ratio1,paste(figurename,'low',sep = '_'))
    p2 <- get_normal_pvalue(ratio2,paste(figurename,'mod',sep = '_'))
    p3 <- get_normal_pvalue(ratio3,paste(figurename,'high',sep = '_'))
    return(rbind(p1,p2,p3))
  }else{
    # plot ratio distribution
    g1 <- plot_distribution(ratio1,NULL,left_line,right_line,mid_shift = 1.15)
    g2 <- plot_distribution(ratio2,NULL,left_line,right_line,mid_shift = 1.15)
    g3 <- plot_distribution(ratio,NULL,left_line,right_line,mid_shift = 1.12)
    jpeg(figurename, width = 22, height = 20, units = 'in', res = 300)
    multiplot(g1,g2,g3,cols=1)
    dev.off()
  }
}

#plot_fig16
plot_fig16 <- function(data_dir,results_dir,file_ids,chr_ids,normality_check=F){
  # for most figs
  files = list.files(data_dir)
  files = files[file_ids]
  names(chr_ids) <- files
  normality_check_res <- NULL
  for (i in files){
    file_name = i
    chr_id = chr_ids[i]
    prefix = gsub('.txt','',file_name)
    output_figure = paste(results_dir,prefix,'.jpeg',sep = '')
    datafile=paste(data_dir,file_name,sep = '')
    dataset = read.table(datafile,header=T,sep="\t",as.is = T)
    chr = dataset[,1]
    geneid = dataset[,2]
    control_exp = dataset[,4]
    treatment_exp = dataset[,3]
    cis_genes= geneid[chr==LETTERS[chr_id]]
    treatment_info = list("control_name" = "control", "control_counts" = control_exp, "treatment_name" = "treatment", "treatment_counts" = treatment_exp,"GeneID" = geneid,"cis_genes" = cis_genes,"split2cis"=1,"outputfile"=output_figure,"left_line" = 0.5,"right_line" = 2)
    if(normality_check){
      normality_check_res <- rbind(normality_check_res,normality_check(treatment_info))
    }else{
      get_result_by_treatment_name(treatment_info)
    }
  }
  if(!is.null(normality_check_res)){
    return(normality_check_res)
  }
}

#fig10
get_df <- function(x,type='gbM'){
  percent_all <- apply(x,2,function(x)mean(x==type))
  df1 <- cbind.data.frame(names(percent_all),percent_all,stringsAsFactors=F)
  colnames(df1) <- c('condition','percentage')
  df1$percentage <- round(df1$percentage,3)
  df1 <- df1[c("diploid","triploid","tetraploid","trisomy1","trisomy2","trisomy3","trisomy4","trisomy5"),]
  df1 <- within(df1,condition <- factor(condition,levels=df1$condition))
  return(df1)
}

get_df2 <- function(x){
  di_m_m <- mean(x[,1]=='gbM' & x[,2]=='gbM')
  di_m_um <- mean(x[,1]=='gbM' & x[,2]=='UM')
  di_um_m <- mean(x[,1]=='UM' & x[,2]=='gbM')
  di_um_um <- mean(x[,1]=='UM' & x[,2]=='UM')
  
  percent_all <- c(di_m_m,di_m_um,di_um_m,di_um_um)
  df1 <- cbind.data.frame(c('Diploid_gbm_condition_gbm','Diploid_gbm_condition_um','Diploid_um_condition_gbm','Diploid_um_condition_um'),
                          percent_all,stringsAsFactors=F)
  colnames(df1) <- c('condition','percentage')
  df1$percentage <- round(df1$percentage,3)
  df1 <- within(df1,condition <- factor(condition,levels=df1$condition))
  return(df1)
}

get_df3 <- function(x,type='gbM'){
  df1 <- NULL
  x <- x[,grep('trisomy',colnames(x))]
  for(i in colnames(x)){
    chr_id <- gsub('trisomy','',i)
    cis_ind <- grep(paste('^AT',chr_id,sep = ''),rownames(x))
    trans_ind <- grep(paste('^AT',chr_id,sep = ''),rownames(x),invert = T)
    df1 <- rbind(df1,c(i,'cis',mean(x[cis_ind,i]==type)),c(i,'trans',mean(x[trans_ind,i]==type)))
  }
  df1[,3] <- round(as.numeric(df1[,3]),3)
  colnames(df1) <- c('sample','type','percentage')
  df1 <- as.data.frame(df1,stringsAsFactors=F)
  df1$sample <- as.factor(df1$sample)
  df1$type <- as.factor(df1$type)
  df1$percentage <- as.numeric(df1$percentage)
  return(df1)
}

df2barplot <- function(df){
  library(ggplot2)
  
  g <- ggplot(data=df, aes(x=condition, y=percentage)) +geom_bar(stat="identity", fill="steelblue",width = 0.5)+
    geom_text(aes(label=percentage), vjust=1.6, color="black", size=11)+labs(x = "")+
    theme(axis.text=element_text(size=35),axis.title=element_text(size=35,face="bold"))
  return(g)
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

plot_fig9 <- function(m){
  p_um <- ggplot(data=m, aes(x=sample, y=percentage, fill=type)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+labs(x = "")+
    geom_text(aes(label=percentage), vjust=1.6, color="black",position = position_dodge(0.9), size=9)+
    theme(axis.text=element_text(size=35),axis.title=element_text(size=35,face="bold"),
          legend.title=element_text(size=32),legend.text=element_text(size=32))+ scale_fill_brewer(palette="Blues")
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
