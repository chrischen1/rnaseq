source('../../de_rnaseq.R')
synapseLogin(rememberMe = T)
syn.local   = "~/data/"

# for MCF7
count_file = "syn10155359"
mata_file  = "syn10155362"

# data cleaning
counts <- read.delim(resolve.filename(count_file,syn.local = syn.local),as.is = T,row.names = 1)
colnames(counts) <- gsub('.+_S(\\d+)_R.+','\\1',colnames(counts))

# read meta data and change it to starndard format
sample_annotations <- read.delim(resolve.filename(mata_file,syn.local = syn.local), stringsAsFactors=FALSE)
sample_annotations <- sample_annotations[sample_annotations$Treatment != 'N/A',]
counts <- counts[,colnames(counts)%in%sample_annotations$Sample]
counts <- counts[apply(counts,1,min)>4,]

trt <- strsplit(sample_annotations$Treatment,split = ' ')
samples <- data.frame(cbind('well'=sample_annotations$Sample,'CellLine'=sample_annotations$Cell.line,'DrugName'=unlist(lapply(trt,function(x)x[2])),
                            'Conc'=unlist(lapply(trt,function(x)x[1])),'Time'=gsub(' hours','',sample_annotations$Time.point)),stringsAsFactors = F)
samples$DrugName[is.na(samples$DrugName)] <- '-'
samples$Conc[samples$Conc=='ctrl'] <- '0.0'
samples$Conc <- as.numeric(samples$Conc)
samples$ctrl <- samples$Conc==0
samples$Time[samples$ctrl] <- 0

# we want to use control in any time for both 6h and 24h group, so we calculate them separately in different group table
#24h
samples_24 <- samples[samples$Time==24 | samples$Time==0,]
samples_24$Time[samples_24$Time==0] <-24
grp_table_24 <- data.frame(cbind('group' = paste(samples_24$CellLine,samples_24$Time,sep = '_'),
                                 'condition' = paste(samples_24$CellLine,samples_24$DrugName,samples_24$Conc,samples_24$Time,sep = '_'),
                                 'control' = samples_24$ctrl),stringsAsFactors = F)
rownames(grp_table_24) <- samples_24$well

cnt_24 <- counts[,samples_24$well]
raw_result_24 <- edgeR_wrapper(cnt_24,grp_table_24,combine_fdr = F)

#same for 6h
samples_6 <- samples[samples$Time==6 | samples$Time==0,]
samples_6$Time[samples_6$Time==0] <-6
grp_table_6 <- data.frame(cbind('group' = paste(samples_6$CellLine,samples_6$Time,sep = '_'),
                                'condition' = paste(samples_6$CellLine,samples_6$DrugName,samples_6$Conc,samples_6$Time,sep = '_'),
                                'control' = samples_6$ctrl),stringsAsFactors = F)
rownames(grp_table_6) <- samples_6$well

cnt_6 <- counts[,samples_6$well]
raw_result_6 <- edgeR_wrapper(cnt_6,grp_table_6,combine_fdr = F)

cnt_sybl <- ens2symbol(unique(c(rownames(raw_result_6)),rownames(raw_result_24)))
cnt_sybl <- cnt_sybl[cnt_sybl$ensembl_gene_id!=''&cnt_sybl$hgnc_symbol!='',]
genes <- cnt_sybl$hgnc_symbol
names(genes) <- cnt_sybl$ensembl_gene_id

logFC_raw <- cbind(raw_result_6$logFC,raw_result_24$logFC)
pval_raw <- cbind(raw_result_6$pmat,raw_result_24$pmat)
fdr_raw <- cbind(raw_result_6$fdr_mat,raw_result_24$fdr_mat)

logFC_ens <- cbind('SYMBOL'=genes[rownames(logFC_raw)],logFC_raw)
pval_ens <- cbind('SYMBOL'=genes[rownames(pval_raw)],pval_raw)
fdr_ens <- cbind('SYMBOL'=genes[rownames(fdr_raw)],fdr_raw)
rownames(logFC_ens) <- rownames(pval_ens) <- rownames(fdr_ens) <- rownames(logFC_raw)
logFC_ens_mcf7 <- logFC_ens
pval_ens_mcf7 <- pval_ens
fdr_ens_mcf7 <- fdr_ens


# for other cell lines in CDK46_response_201608

count_file = "syn10164324"
mata_file  = "syn10155361"

# data cleaning
counts <- read.delim(resolve.filename(count_file,syn.local = syn.local),as.is = T,row.names = 1)

# read meta data and change it to starndard format
sample_annotations <- read.delim(resolve.filename(mata_file,syn.local = syn.local), stringsAsFactors=FALSE)

# remove sample with bad quality(HS578T control 24h)
counts <- counts[,-8]
counts <- counts[apply(counts,1,min)>4,]
sample_annotations <- sample_annotations[-33,]

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

genes_intersect <- intersect(rownames(fdr_raw),rownames(logFC_ens_mcf7))

logFC_ens_final <- cbind(logFC_ens_mcf7[genes_intersect,],logFC_raw[genes_intersect,])
pval_ens_final <- cbind(pval_ens_mcf7[genes_intersect,],pval_raw[genes_intersect,])
fdr_ens_final <- cbind(fdr_ens_mcf7[genes_intersect,],fdr_raw[genes_intersect,])

write.csv(logFC_ens_final,'LogFC_perturbation.csv')
write.csv(pval_ens_final,'pvalue_perturbation.csv')
write.csv(fdr_ens_final,'FDR_perturbation.csv')
