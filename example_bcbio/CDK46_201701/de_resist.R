source('../../de_rnaseq.R')
synapseLogin(rememberMe = T)
syn.local   = "~/data/"

count_file = "syn10155359"
mata_file  = "syn10155362"

# data cleaning
counts <- read.delim(resolve.filename(count_file,syn.local = syn.local),as.is = T,row.names = 1)
colnames(counts) <- gsub('.+_S(\\d+)_R.+','\\1',colnames(counts))

# read meta data and change it to starndard format
sample_annotations <- read.delim(resolve.filename(mata_file,syn.local = syn.local), stringsAsFactors=FALSE)
sample_annotations <- sample_annotations[sample_annotations$Treatment == 'N/A',]
sample_annotations <- sample_annotations[sample_annotations$Sample!=30,] #remove low quality sample
counts <- counts[,colnames(counts)%in%sample_annotations$Sample]
counts <- counts[apply(counts,1,min)>4,]

trt <- strsplit(sample_annotations$Treatment,split = ' ')
samples <- data.frame(cbind('well'=sample_annotations$Sample,'CellLine'=gsub('(.+)\\-PR| .+','\\1',sample_annotations$Cell.line),'DrugName'='-',
                 'Conc'=sample_annotations$Cell.line,'Time'=0),stringsAsFactors = F)

grp_table <- data.frame(cbind('group' = gsub('(.+)\\-PR| .+','\\1',sample_annotations$Cell.line),
                                 'condition' = sample_annotations$Cell.line,'control' = F),stringsAsFactors = F)
grp_table$control[grep('parent',sample_annotations$Cell.line)] <- T
rownames(grp_table) <- samples$well
grp_table$condition <- gsub(' on| off','',grp_table$condition)

# remove MCF7 as lack of replicates
grp_table <- grp_table[grp_table$group!='MCF7',]
cnt <- counts[,rownames(grp_table)]
raw_result <- edgeR_wrapper(cnt,grp_table)

logFC_raw <- raw_result$logFC
pval_raw <- raw_result$pmat
fdr_raw <- raw_result$fdr_mat

cnt_sybl <- ens2symbol(rownames(logFC_raw))
cnt_sybl <- cnt_sybl[cnt_sybl$ensembl_gene_id!=''&cnt_sybl$hgnc_symbol!='',]
genes <- cnt_sybl$hgnc_symbol
names(genes) <- cnt_sybl$ensembl_gene_id

logFC_ens <- cbind('SYMBOL'=genes[rownames(logFC_raw)],logFC_raw)
pval_ens <- cbind('SYMBOL'=genes[rownames(pval_raw)],pval_raw)
fdr_ens <- cbind('SYMBOL'=genes[rownames(fdr_raw)],fdr_raw)
rownames(logFC_ens) <- rownames(pval_ens) <- rownames(fdr_ens) <- rownames(logFC_raw)

write.csv(logFC_ens,'LogFC_resist.csv')
write.csv(pval_ens,'pvalue_resist.csv')
write.csv(fdr_ens,'FDR_resist.csv')

