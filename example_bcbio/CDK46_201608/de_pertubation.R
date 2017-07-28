source('../de_rnaseq.R')
count_file = 'syn10164324'
mata_file  = 'syn10155361'
tx2gene_file = 'syn10084072'
out_folder = './'
syn.local   = "~/data/"

# data cleaning
counts <- read.delim(resolve.filename(count_file,syn.local = syn.local),as.is = T,row.names = 1)

# get counts from salmon, run on orchestra
tx2gene <- read.csv(resolve.filename(tx2gene_file,syn.local = syn.local),as.is = T,header = F)

# read meta data and change it to starndard format
sample_annotations <- read.delim(resolve.filename(mata_file,syn.local = syn.local), stringsAsFactors=FALSE)

# remove sample with bad quality(HS578T control 24h)
sample_annotations <- sample_annotations[-33,]

counts <- counts[,-8]
counts <- counts[apply(counts,1,min)>4,]

# we want to use control in any time for both 6h and 24h group, so we calculate them separately in different group table
#24h
samples_24 <- sample_annotations[sample_annotations$Time==24 | sample_annotations$Conc==0,]
samples_24$Time <-24
grp_table_24 <- data.frame(cbind('group' = paste(samples_24$CellLine,samples_24$Time,sep = '_'),
                                 'condition' = paste(samples_24$CellLine,samples_24$DrugName,samples_24$Conc,samples_24$Time,sep = '_'),
                                 'control' = samples_24$DrugName=='-'),stringsAsFactors = F)
rownames(grp_table_24) <- samples_24$well

cnt_24 <- counts[,samples_24$well]

raw_result_24 <- edgeR_wrapper(cnt_24,grp_table_24)
logfc24 <- raw_result_24$logFC

#same for 6h
samples_6 <- sample_annotations[sample_annotations$Time==6 | sample_annotations$Conc==0,]
samples_6$Time <-6
grp_table_6 <- data.frame(cbind('group' = paste(samples_6$CellLine,samples_6$Time,sep = '_'),
                                'condition' = paste(samples_6$CellLine,samples_6$DrugName,samples_6$Conc,samples_6$Time,sep = '_'),
                                'control' = samples_6$DrugName=='-'),stringsAsFactors = F)
rownames(grp_table_6) <- samples_6$well

cnt_6 <- counts[,samples_6$well]
raw_result_6 <- edgeR_wrapper(cnt_6,grp_table_6)

logFC_raw <- cbind(raw_result_6$logFC,raw_result_24$logFC[rownames(raw_result_6$logFC),])
pval_raw <- cbind(raw_result_6$pmat,raw_result_24$pmat[rownames(raw_result_6$logFC),])
fdr_raw <- cbind(raw_result_6$fdr_mat,raw_result_24$fdr_mat[rownames(raw_result_6$logFC),])

cnt_sybl <- ens2symbol(rownames(logFC_raw))
cnt_sybl <- cnt_sybl[cnt_sybl$ensembl_gene_id!=''&cnt_sybl$hgnc_symbol!='',]
genes <- cnt_sybl$hgnc_symbol
names(genes) <- cnt_sybl$ensembl_gene_id

logFC <- logFC_raw[rownames(logFC_raw)%in%names(genes),]
rownames(logFC) <- genes[rownames(logFC)]

pval <- pval_raw[rownames(pval_raw)%in%names(genes),]
rownames(pval) <- genes[rownames(pval)]

fdr <- fdr_raw[rownames(fdr_raw)%in%names(genes),]
rownames(fdr) <- genes[rownames(fdr)]

logFC_ens <- cbind('SYMBOL'=genes[rownames(logFC_raw)],logFC_raw)
pval_ens <- cbind('SYMBOL'=genes[rownames(pval_raw)],pval_raw)
fdr_ens <- cbind('SYMBOL'=genes[rownames(fdr_raw)],fdr_raw)
rownames(logFC_ens) <- rownames(pval_ens) <- rownames(fdr_ens) <- rownames(logFC_raw)

write.csv(logFC_ens,paste(out_folder,'logFC_ens.csv',sep = ''))
write.csv(pval_ens,paste(out_folder,'pval_ens.csv',sep = ''))
write.csv(fdr_ens,paste(out_folder,'fdr_ens.csv',sep = ''))
