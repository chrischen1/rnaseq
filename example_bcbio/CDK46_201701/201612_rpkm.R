source('../../de_rnaseq.R')

#get combined.sf if not exist
# change working directory to /work
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
write.table(sf,'../final/runDate_201612_de/combined.sf',sep='\t')

# change working directory to  201612_de/final/runDate_201612_de/
sf_file = 'combined.sf'
tx_file = 'tx2gene.csv'

sf <- read.delim(sf_file,as.is = T)
tx <- read.csv(tx_file,as.is = T,header=FALSE)
sample_info <- read.delim(resolve.filename('syn10155362',syn.local = syn.local),as.is = T)

samples <- paste(sample_info$Cell.line,sample_info$Treatment,sample_info$Time.point)
names(samples) <- sample_info$Sample
samples <- gsub(' N/A N/A','',samples)
rpkm <- tpm2rpkm(sf,tx)
colnames(rpkm) <- samples[gsub('.+S(\\d+)_.+','\\1',colnames(rpkm))]

genes <- ens2symbol(rownames(rpkm))
genes <- genes[genes$hgnc_symbol!='',]
rpkm2 <- rpkm[rownames(rpkm)%in%genes$ensembl_gene_id,]
gene_v <- genes$hgnc_symbol
names(gene_v) <- genes$ensembl_gene_id
rpkm2 <- cbind.data.frame(rpkm2,'symbol'=gene_v[rownames(rpkm2)],stringsAsFactors=F)
rpkm3 <- rpkm2 %>% group_by(symbol) %>% summarise_each(funs(sum))
rpkm_symbol <- rpkm3[,-1]
rownames(rpkm_symbol) <- rpkm3$symbol

write.table(log2(rpkm+1),'./201612_rpkm_ens.tsv',sep = '\t')
write.table(log2(rpkm_symbol+1),'./201612_rpkm_symbol.tsv',sep = '\t')
