source('../../de_rnaseq.R')
synapseLogin(rememberMe = T)
syn.local   = "~/data/"

#get combined.sf if not exist
# change working directory to  201604_de/work/
sf_files <- list.files(path="./",pattern='*\\.sf',recursive=T)
sf_info <- NULL
for(i in sf_files){
  si <- read.delim(i,as.is = T)
  si <- cbind(si,'sample'=gsub('salmon/(.+)/quant/.+','\\1',i),'id'=si$Name)
  sf_info <- rbind(sf_info,si)
}
sf <- sf_info
colnames(sf) <- c('name','length','effectiveLength','tpm','numreads','sample','id')
# change this according to project name:
write.table(sf,'../final/runDate_201604_de/combined.sf',sep='\t')

# change working directory to  201604_de/final/runDate_201604_de/
sf_file = 'combined.sf'
tx_file = 'tx2gene.csv'

sf <- read.delim(sf_file,as.is = T)
tx <- read.csv(tx_file,as.is = T,header=FALSE)
rpkm <- tpm2rpkm(sf,tx)

sample_info <- read.delim(resolve.filename('syn10084520',syn.local = syn.local),as.is = T)
samples <- sample_info$cellines
rep_tag <- sample_info$rep_tag
names(samples) <- names(rep_tag) <- sample_info$sample_id
colnames(rpkm) <- gsub('-bio-replic','',colnames(rpkm),fixed = T)
colnames(rpkm) <- gsub('-tech-replic','',colnames(rpkm),fixed = T)
colnames(rpkm) <- samples[colnames(rpkm)]

genes <- ens2symbol(rownames(rpkm))
genes <- genes[genes$hgnc_symbol!='' & !(genes$ensembl_gene_id %in% genes$ensembl_gene_id[duplicated(genes$ensembl_gene_id)]) & 
                 !(genes$hgnc_symbol %in% genes$hgnc_symbol[duplicated(genes$hgnc_symbol)]),]
gene_name <- genes$hgnc_symbol
names(gene_name) <- genes$ensembl_gene_id

#merge replicates
rpkm_merge_rep <- NULL
sample_names <- c()
for(i in unique(colnames(rpkm))){
  sample_names <- c(sample_names,i)
  if(sum(colnames(rpkm)==i)==1){
    rpkm_merge_rep <- rbind(rpkm_merge_rep,rpkm[,colnames(rpkm)==i])
  }else{
    rpkm_merge_rep <- rbind(rpkm_merge_rep,apply(rpkm[,colnames(rpkm)==i],1,mean))
  }
}
rownames(rpkm_merge_rep) <- sample_names
rpkm_merge_rep <- t(rpkm_merge_rep)
rpkm_merge_rep_symbol <- rpkm_merge_rep[rownames(rpkm_merge_rep)%in%names(gene_name),]
rownames(rpkm_merge_rep_symbol) <- gene_name[rownames(rpkm_merge_rep_symbol)]

rpkm_merge_rep_final <- cbind('id'=rownames(rpkm_merge_rep),log2(rpkm_merge_rep+1))
rpkm_merge_rep_symbol_final <- cbind('id'=rownames(rpkm_merge_rep_symbol),log2(rpkm_merge_rep_symbol+1))

write.table(rpkm_merge_rep_final,'./rpkm_basal_merge_rep_ens.tsv',sep = '\t',row.names = F,quote = F)
write.table(rpkm_merge_rep_symbol_final,'./RNAseq-rpkm.tsv',sep = '\t',row.names = F,quote = F)
