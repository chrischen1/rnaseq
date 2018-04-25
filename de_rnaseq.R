## de_rnaseq.R - functions for the differential expression analysis
##
## LSP RNAseq bcbio pipeline 
## by Artem Sokolov, Chris Chen, et al.

## Retrieves count file and group information file from command line arguments, 
## Returns a named list of values which is used by the main() function in run_de.R
get_args <- function(){
  library(optparse)
  ## Define available options
  option_list = list(
    make_option(c("-c", "--count"), type="character", default=NULL,
                help="Path to .count file from bcbio output, which is a ensumbl ID by sample ID matrix", metavar="character"),
    make_option(c("-a", "--annotation"), type="character", default=NULL, 
                help="Path to group information file, which is a dataframe with 3 columns: group, condition and control
                \n\tgroup: contains information which treatment samples will be compared against control cases in each group
                \n\tcondition: indicates type of treatment, replicates have same condition
                \n\tcontrol: TRUE for controls and FALSE for treatments
                \n\torder of samples in annotation must be the same as samples in count table", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="Path to save differential analysis results", metavar="character"),
    make_option(c("-p", "--pairwise"), type="logical", default=TRUE, 
                help="If the P-values and FDR are given pairwise or as ANOVA-like test for any differences", metavar="TRUE/FALSE"),
    make_option(c("-s", "--symbol"), type="logical", default=FALSE, 
                help="If gene symbols will be added to the output", metavar="TRUE/FALSE")
    )
  
  ## Parse the arguments
  opt_parser = OptionParser(option_list=option_list)
  argv = parse_args(opt_parser)
  
  ## Basic verification
  if (is.null(argv$count) || is.null(argv$annotation) || is.null(argv$output)){
    print_help(opt_parser)
    stop("Count table, annotation and output path must be provided.\n 
         usage: Rscript run_de.R -c path/to/rnaseq.count -a path/to/group_info.txt -o path/to/output", call.=FALSE)
  }
  
  return( argv )
  }

#' transform TPM to RPKM
#'
#' @param combined output file end with .combined from bcbio.
#' @param tx2gene output file which maps ensumble ID to gene from bcbio.
#' @param spikes a vector of string defining the name of spikes.
#' @return p by n matrix for p genes across n samples
tpm2rpkm <- function(combined,tx2gene,spikes = NULL){
  library(reshape2)
  library(dplyr)
  colnames(combined) <- tolower(colnames(combined))
  gene_mapping <- cbind('transcript'= c(tx2gene$V1,spikes$GenBank),'gene' = c(tx2gene$V2,spikes$ERCC_ID))
  genes <- gene_mapping[,2]
  names(genes) <- gene_mapping[,1]
  lib_size <- data.frame('numreads'=combined$numreads,'sample'=combined$sample)
  x <- lib_size %>% group_by(sample) %>% summarise_each(funs(sum))
  scale_factor <- x$numreads/1000000
  names(scale_factor) <- x$sample
  
  combined$RPM <- combined$numreads/scale_factor[combined$sample]
  combined$RPKM <- combined$RPM/(combined$effectivelength/1000)
  combined$gene <- genes[combined$name]
  
  rpkm_combined <- data.frame('sample'=combined$sample,'gene'=combined$gene,'RPKM'=combined$RPKM)
  rpkm_combined_gene <- rpkm_combined %>% group_by(sample,gene)%>% summarise_each(funs(sum))
  
  rpkm_raw <- acast(rpkm_combined_gene,gene~sample)
  return(rpkm_raw[-nrow(rpkm_raw),])
}

#' transform TPM to RPKM
#'
#' @param combined output file end with .combined from bcbio.
#' @param tx2gene output file which maps ensumble ID to gene from bcbio.
#' @param spikes a vector of string defining the name of spikes.
#' @return p by n matrix for p genes across n samples
sf2tpm <- function(combined,tx2gene,spikes = NULL){
  library(reshape2)
  library(dplyr)
  colnames(combined) <- tolower(colnames(combined))
  gene_mapping <- cbind('transcript'= c(tx2gene$V1,spikes$GenBank),'gene' = c(tx2gene$V2,spikes$ERCC_ID))
  genes <- gene_mapping[,2]
  names(genes) <- gene_mapping[,1]
  combined$gene <- genes[combined$Name]
  combined2 <- combined[!is.na(combined[,'gene']),]
  tpm_combined <- data.frame('sample'=combined2$sample,'gene'=combined2$gene,'tpm_raw'=combined2$tpm)
  tpm_combined_gene <- tpm_combined %>% group_by(sample,gene)%>% summarise_each(funs(sum))
  
  tpm_raw <- acast(tpm_combined_gene,gene~sample)
  return(tpm_raw)
}


#' get length of genes given a list of gene names
#'
#' @param gene_id_type types of gene ID, see ensembl bioMart package for possible types
#' @param dataset name of ensembl dataset
#' @return a named list of gene lengths, duplicate gene length records will be averaged
get_gene_length <- function(gene_names,gene_id_type='ensembl_gene_id',dataset="hsapiens_gene_ensembl"){
  library(biomaRt)
  library(reshape2)
  library(dplyr)
  
  ensembl <- useEnsembl(biomart="ensembl", dataset=dataset)
  target_gene_raw <- getBM(attributes=c(gene_id_type,'transcript_length'),filters = gene_id_type, values = gene_names, mart = ensembl)
  target_gene <- target_gene_raw %>% group_by(ensembl_gene_id)%>% summarise_all(funs(mean))
  l <- target_gene$transcript_length
  names(l) <- target_gene$ensembl_gene_id
  return(l)
}

#' transform counts to TPM
#'
#' @param m p by n counts matrix for p genes across n samples
#' @param gene_length a named list of gene lengths
#' @return p by n TPM matrix for p genes across n samples
counts2tpm <- function(m,gene_length){
  m <- m[names(gene_length),]
  rpk <- m/gene_length
  scaling <- colSums(rpk)/1000000
  for(i in 1:ncol(rpk)){
    rpk[,i] <- rpk[,i]/scaling[i]
  }
  rownames(rpk) <- names(gene_length)
  return(rpk)
}


#' get hgnc_symbol from ensembl_gene_id  
#'
#' @param ens vector of ensembl_gene_ids.
#' @return a dataframe with 2 columns: ensembl_gene_id and hgnc_symbol
ens2symbol <- function(ens){
  library(biomaRt)
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  target_gene <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),filters = 'ensembl_gene_id', values = ens, mart = ensembl)
  return(target_gene)
}

#' get gene_id1 from gene_id2
#'
#' @param ids vector of gene_ids.
#' @param gene_id1 format of original gene id, must be valid filters name in Ensembl
#' @param gene_id2 format of destination gene id, must be valid attributes name in Ensembl
#' @return a dataframe with 2 columns: gene_id1 and gene_id2
gene_id_mapping <- function(ids,gene_id1='ensembl_gene_id',gene_id2='hgnc_symbol',dataset="hsapiens_gene_ensembl"){
  library(biomaRt)
  ensembl <- useEnsembl(biomart="ensembl", dataset=dataset)
  target_gene <- getBM(attributes=c(gene_id1,gene_id2),filters = gene_id1, values = ids, mart = ensembl)
  return(target_gene)
}

#' conversion a gene by sample matrix from gene_id1 from gene_id2 by names in t[,1] to names in t[,2]
#'
#' @param m a gene by sample matrix
#' @param t conversion table
#' @return a gene by sample matrix with new gene_id
gene_matrix_conversion_by_table <- function(m,t){
  library(reshape2)
  library(dplyr)
  genes <- t[,2]
  names(genes) <- t[,1]
  m <- cbind.data.frame(m,'name'=genes[rownames(m)],stringsAsFactors=F)
  m <- m[!is.na(m$name),]
  m2 <- m %>% group_by(name)%>% summarise_each(funs(sum))
  m3 <- as.matrix(m2[,-1])
  rownames(m3) <- m2$name
  return(m3)
}


#' conversion a gene by sample matrix from gene_id1 from gene_id2
#'
#' @param m a gene by sample matrix
#' @param gene_id1 format of original gene id, must be valid filters name in Ensembl
#' @param gene_id2 format of destination gene id, must be valid attributes name in Ensembl
#' @return a gene by sample matrix with new gene_id
gene_matrix_conversion <- function(m,gene_id1='ensembl_gene_id',gene_id2='hgnc_symbol',dataset="hsapiens_gene_ensembl"){
  
  target_genes <- gene_id_mapping(rownames(m),gene_id1,gene_id2,dataset)
  target_genes <- target_genes[target_genes[,1]!=''&target_genes[,2]!='',]
  gene_matrix_conversion_by_table(m,target_genes)
}


## Resolves a filename by downloading the file if it's a synapse ID
## Returns a filename that can be directly used for loading by, e.g., read.delim
resolve.filename <- function( fn, syn.local   = "~/data/")
{
  library(synapseClient)
  if( substr( fn, 0, 3 ) == "syn" )
  {
    dir.create(syn.local,showWarnings = F)
    s <- synGet( fn, downloadLocation = syn.local )
    return( s@filePath )
  }
  return( fn )
}

#' generate .csv file used by bcbio
#' 
#' @param sample_path path of .fastq files
#' @return a csv file contain basic sample meta info required for bcbio
get_sample_csv <- function(sample_path){
  x=grep('\\.fastq',list.files(sample_path),value = T)
  y=gsub('\\.fastq','',x)
  z=cbind('samplename'=y,'description'=y)
  write.csv(z,paste(sample_path,'samples.csv',sep = ''),row.names = F,quote = F)
}

#' wrapper for getting fold change, pvalue and FDR, by per cell line per time point
#' 
#' @param cnt p by n matrix for p genes across n samples
#' @param grp_table dataframe with 3 columns: group, condition and control
#'  group: contains information which treatment samples will be compared against control cases in each group
#'  condition: indicates type of treatment, replicates have same condition, do NOT use numbers only.
#'  control: TRUE for controls and FALSE for treatments
#'  order of well in samples annotation must be the same as the columns in count table
#' @param combine_fdr T for combine FDR and p-values with group and F for compute pairwisely
#' @param w n by p matrix for n samples and p factors for batch effect correction from RUVSeq
#' @param CommonDisp and TagwiseDisp used internally for passing overal dispersion to comparisons without replicates
#' @return list of 3 if combine_fdr = F: pmat,fdr_mat and logFC: all are p by m matrix for p genes across m types of treatments
#'         p by m+4 matrix for p genes across m types of treatments and p-value, LR,logCPM and FDR
edgeR_wrapper <- function(cnt,grp_table,combine_fdr = F,w = NULL,CommonDisp = NULL,TagwiseDisp = NULL){
  library(edgeR)
  if(sum(rownames(grp_table) %in% colnames(cnt)) < nrow(grp_table)){
    warning('rownames in group table not compatible with colnames in count table')
    return(0)
  }
  cnt <- cnt[,rownames(grp_table)]
  design <- model.matrix(~condition,data = grp_table)
  # add RUV batch effect correction when w exists
  if(!is.null(w))  design <- cbind(design,w)
  y <- DGEList(counts=cnt, group=grp_table$condition)
  # Calculate overall dispersions when called first time
  if(is.null(CommonDisp)){
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    CommonDisp <- y$common.dispersion
    TagwiseDisp <- y$tagwise.dispersion
  }
  if(length(grp_table$condition)==unique(length(grp_table$condition))){
    # When both control and treatment lacking replicates, use overall dispersion instead
    y$common.dispersion <- CommonDisp
    y$tagwise.dispersion <- TagwiseDisp
  }else{
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
  }
  # only anova-like FDR/Pvalues is required
  if(combine_fdr){
    y <- calcNormFactors(y)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef=2:(ncol(design)))
    lrt_tab <- topTags(lrt,n = Inf)$table[rownames(cnt),]
    colnames(lrt_tab) <- gsub('logFC.condition','',colnames(lrt_tab))
    return(lrt_tab)
  }
  # pairwise FDR/Pvalues is required
  p_mat <- fdr_mat <- logFC <- NULL
  col_names <- c()
  for(i in unique(grp_table$group)){
    grp_table_i <- grp_table[grp_table$group==i,]
    ctr_row <- rownames(grp_table_i)[grp_table_i$control==T]
    for (j in unique(grp_table_i$condition[grp_table_i$control!=T])){
      j_row <- rownames(grp_table_i)[grp_table_i$condition==j]
      grp_new <- rbind(grp_table[c(ctr_row,j_row),])
      cnt_new <- cnt[,rownames(grp_new)]
      result_new <- edgeR_wrapper(cnt_new,grp_new,combine_fdr = T,CommonDisp = CommonDisp,TagwiseDisp = TagwiseDisp)
      if(is.null(p_mat)){
        p_mat <- result_new$PValue
        fdr_mat <- result_new$FDR
        logFC <- result_new$logFC
        if(length(unique(grp_table$condition))==2){
          return_mat <- cbind('pval'=p_mat,'fdr'=fdr_mat,'logFC'=logFC)
          rownames(return_mat) <- rownames(cnt_new)
          return(return_mat)
        }
      }else{
        p_mat <- cbind(p_mat,result_new$PValue)
        fdr_mat <- cbind(fdr_mat,result_new$FDR)
        logFC <- cbind(logFC,result_new$logFC)
      }
      col_names <- c(col_names,j)
    }
  }
  colnames(p_mat) <- colnames(fdr_mat) <- colnames(logFC) <- col_names
  rownames(p_mat) <- rownames(fdr_mat) <- rownames(logFC) <- rownames(cnt)
  return(list('pmat'=p_mat,'fdr_mat'=fdr_mat,'logFC'=logFC))
}

# examine if there is task running on slurm
slurm_running <- function(job_name){
  sacct_out <- system('sacct |grep ING',intern = T)
  return(length(grep(job_name,sacct_out))>0)
}

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

kegg_pathway_enrichment <- function(gene_names,organism='mmu',keyType = 'uniprot',pAdjustMethod='none',pvalueCutoff = 0.05){
  library(clusterProfiler)
  result<- enrichKEGG(gene_names,organism = organism,keyType = keyType,pAdjustMethod=pAdjustMethod,pvalueCutoff = pvalueCutoff)
  pathway_slt <- result@result
  
}

plot_pathway_kegg <- function(gene_names,log_fc,pathway.id,out.suffix='',organism='mmu',gene.idtype='entrez'){
  library(pathview)
  names(log_fc) <- gene_names
  pv.out <- pathview(gene.data =log_fc, pathway.id = pathway.id,out.suffix=out.suffix,species = organism, kegg.native = F,gene.idtype = gene.idtype,pdf.size=c(12,12))

}



