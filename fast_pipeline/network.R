source('https://raw.githubusercontent.com/chrischen1/rnaseq/master/fast_pipeline/plot_utils.R')

infile = '~/Dropbox/MU/workspace/mike20200308/list_03022020.csv'
outpath = '~/Dropbox/MU/workspace/mike20200308/'

library(minet)
library(GNET2)
library(bio3d)
library(dplyr)
library(reshape2)
library(igraph)

###############
x <- read.csv(infile)
exp1 <- x[,11:22]
exp2 <- log2(exp1+1)
sample_x <- gsub('\\.+I.+','',colnames(exp2))
sample_x <- gsub('\\.','_',sample_x)
condition_x <- gsub('\\d+','',sample_x)
colnames(exp2) <- sample_x
protein_names <- gsub('\\|.+','',as.character(x$Accession))
rownames(exp2) <- protein_names
protein_symbol <- gsub('_MOUSE','',gsub('.+\\|','',as.character(x$Accession)))
names(protein_symbol) <- protein_names
##########################

# fa_out <- ''
# for (i in 1:nrow(x)) {
#   fa_out <- paste0(fa_out,'>',gsub('.+\\|','',as.character(x$Accession[i])),'\n')
#   fa_out <- paste0(fa_out,uniprot(protein_names[i])$sequence,'\n')
# 
# }
# writeLines(fa_out,paste0(outpath,'protein_list.fa'))

mi <- build.mim(t(exp2))

rownames(mi) <- colnames(mi) <- rownames(exp2)
n1 <- aracne(mi)
n2 <- clr(mi)
n3 <- minet(mi)
n4 <- mrnet(mi)

exp_var <- apply(exp2, 1, var)
gnet_res1 <- gnet(exp2,names(exp_var)[exp_var>median(exp_var)],init_method = 'kmeans',max_iter = 20)
gnet_res2 <- gnet(exp2,names(exp_var)[exp_var<=median(exp_var)],init_method = 'kmeans',max_iter = 20)

el1 <- extract_edges(gnet_res1)
el2 <- extract_edges(gnet_res2)

el2a <- el2
el2a$regulator <- el2$target
el2a$target <- el2$regulator
el_all <- rbind.data.frame(el1,el2a,stringsAsFactors = F)
el_all2 <- el_all %>% group_by_(.dots = c('regulator','target')) %>% summarise_all(list('score' = sum))
el_all2$score <- el_all2$score/max(el_all2$score)
el_all3 <- data.frame(el_all2)

n1l <- melt(n1)
n1l <- n1l[n1l$value>0,]
n1l$value <- n1l$value/max(n1l$value)

n2l <- melt(n2)
n2l <- n2l[n2l$value>0,]
n2l$value <- n2l$value/max(n2l$value)

n3l <- melt(n3)
n3l <- n3l[n3l$value>0,]
n3l$value <- n3l$value/max(n3l$value)

n4l <- melt(n4)
n4l <- n4l[n4l$value>0,]
n4l$value <- n4l$value/max(n4l$value)

colnames(el_all3) <- colnames(n1l) <- colnames(n2l) <- colnames(n3l) <- colnames(n4l) <- c('x','y','score')
all_interactions <- unique(rbind.data.frame(el_all3[,1:2],n1l[,1:2],n2l[,1:2],n3l[,1:2],n4l[,1:2]))

g_all <- graph_from_edgelist(as.matrix(all_interactions),directed = F)
el_all4 <- unique(get.edgelist(g_all))


sum_scores <- function(el_all,el_input){
  scores <- rep(0,nrow(el_all))
  for (i in 1:nrow(el_all)) {
    idx <- (el_input$x == el_all[i,1] & el_input$y == el_all4[i,2]) | (el_input$y == el_all[i,1] & el_input$x == el_all[i,2])
    scores[i] <- sum(el_input$score[idx])
  }
  return(scores)
}
s1 <- sum_scores(el_all4,el_all3)
s2 <- sum_scores(el_all4,n1l)
s3 <- sum_scores(el_all4,n2l)
s4 <- sum_scores(el_all4,n3l)
s5 <- sum_scores(el_all4,n4l)

el_all_final <- cbind.data.frame(el_all4,'score'=s1+s2+s3+s4+s5)
el_all_final$score <- el_all_final$score/max(el_all_final$score)
el_all_final2 <- unique(el_all_final)



tiff(paste0(outpath,'pred_interactions.tiff'),width = 2560,height = 1440)
g1 <- plot_network(el_all_final,n = 200,protein_symbol = protein_symbol)
dev.off()

el_all_final2[,1] <- protein_symbol[el_all_final2[,1]]
el_all_final2[,2] <- protein_symbol[el_all_final2[,2]]

write.table(el_all_final2,paste0(outpath,'pred_interactions.txt'),row.names = F,col.names = F,quote = F,sep = '\t')

