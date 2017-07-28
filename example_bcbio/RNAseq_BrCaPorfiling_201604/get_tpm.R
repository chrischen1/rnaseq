source('../de_rnaseq.R')
combined_id = 'syn10084069'
tx2gene_id  = 'syn10084072'
meta_id     = 'syn10084520'
syn.local   = "~/data/"

## Log in to Synapse and retrieve data
synapseLogin(rememberMe = T)
combined <- read.delim(resolve.filename(combined_id,syn.local = syn.local),sep = '\t',as.is = T)
tx2gene  <- read.csv(resolve.filename(tx2gene_id,syn.local = syn.local),header = F,as.is = T)
sample_info <- read.delim(resolve.filename(meta_id,syn.local = syn.local),as.is = T)

# Get TPM results
tpm_raw <- sf2tpm(combined,tx2gene)

# Mapping samples and merge replicates by averaging
tpm <- log2(tpm_raw+1)
samples <- sample_info$cellines
rep_tag <- sample_info$rep_tag
names(samples) <- names(rep_tag) <- sample_info$sample_id
colnames(tpm) <- gsub('-bio-replic','',colnames(tpm),fixed = T)
colnames(tpm) <- gsub('-tech-replic','',colnames(tpm),fixed = T)
tpm_merge <- NULL
samples_merge <- c()
for(i in unique(samples)){
  samples_merge <- c(samples_merge,i)
  col_slct <- colnames(tpm) %in% names(samples)[samples==i]
  if(sum(col_slct)>1){
    new_col <- apply(tpm[,col_slct],1,mean)
  }else{
    new_col <- tpm[,col_slct]
  }
  if(is.null(tpm_merge)){
    tpm_merge <- new_col
  }else{
    tpm_merge <- cbind(tpm_merge,new_col)
  }
}
colnames(tpm_merge) <- samples_merge
tpm_hgnc <- gene_matrix_conversion(tpm_merge,gene_id1='ensembl_gene_id',gene_id2='hgnc_symbol')

tpm_merge2 <- rbind('id'=colnames(tpm_merge),tpm_merge)
tpm_hgnc2 <- rbind('id'=colnames(tpm_hgnc),tpm_hgnc)
write.table(tpm_merge2,'RNAseq-tpm-ens.tsv',sep = '\t',quote = F,col.names = F)
write.table(tpm_hgnc2,'RNAseq-tpm.tsv',sep = '\t',quote = F,col.names = F)
