# merge fastq files of sample sample from 2 runs to get one new fastq file

x=list.files('/groups/sorger/Marc/RNAseq_data/CDK46_response_201608/160824_CM3415-1_NS500233_fastq')
p1='/groups/sorger/Marc/RNAseq_data/CDK46_response_201608/160824_CM3415-1_NS500233_fastq/'
p2='/groups/sorger/Marc/RNAseq_data/CDK46_response_201608/160825_CM3415-2_NS500144_fastq/'
for(i in x){
  i2 <- gsub('20160824','20160825',i)
  i2 <- gsub('CM3415-1','CM3415-2',i2)
  m <- gsub('20160824_(.+)_CM.+','\\1',i)
  system(paste("bsub -q short -W 3:00 'cat ",p1,i,' ',p2,i2,' > /groups/sorger/cchris/fastq_merged_201608/',m,".fastq.gz'",sep = ''))
}