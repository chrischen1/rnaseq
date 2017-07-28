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
cnt <- counts[,rownames(grp_table)]

x <- cbind(cnt[,2]/cnt[,1],cnt[,3]/cnt[,1],cnt[,5]/cnt[,4],cnt[,7]/cnt[,6],cnt[,8]/cnt[,6],cnt[,10]/cnt[,9],cnt[,11]/cnt[,9])
y <- log(x,2)
rownames(y) <- rownames(cnt)
genes <- ens2symbol(ens = rownames(y))
aaa=genes$hgnc_symbol
names(aaa) <- genes$ensembl_gene_id
y <- cbind(y,'SYMBOL'=aaa[rownames(y)])
y <- y[,c(8,1:7)]
colnames(y) <- c('SYMBOL','Hs578T_on','Hs578T_off','MCF7_on','T47D_on','T47D_off','HCC1806_on','HCC1806_off')
write.table(y,'exact_LogFC_resist.csv')
