# environment set
raw_data_path = '/storage/htc/birchlerlab/CENH3_ChIP-seq/data/ChIP_data/'
out_dir = '/storage/htc/birchlerlab/CENH3_ChIP-seq/data/'

filter_data_path = paste0(out_dir,'ChIP_data_trimmed_filtered/')
qc_result_path = paste0(out_dir,'ChIP_data_qc/')
sam_data_path = paste0(out_dir,'sams/')
align_log_path = paste0(out_dir,'align_log/')
bam_data_path = paste0(out_dir,'bams/')
bam_merge_path = paste0(out_dir,'bams_merge/')
peakcalling_results = paste0(out_dir,'peakcalling/')

log_path = paste0(out_dir,'../logs/')

ref_genome_path = '/storage/htc/birchlerlab/CENH3_ChIP-seq/data/B73/index_bowtie2/B73_Bcentremere'
dir.create(qc_result_path)
dir.create(sam_data_path)
dir.create(align_log_path)
dir.create(bam_data_path)
dir.create(bam_merge_path)
dir.create(filter_data_path)
dir.create(peakcalling_results)
dir.create(log_path)

all_samples <- unique(gsub('_R.+','',list.files(raw_data_path)))
setwd(log_path)

#1.Filter low quality reads
for(i in all_samples){
  infile1 <- paste0(raw_data_path,i,'_R1_001.fastq.gz')
  infile2 <- paste0(raw_data_path,i,'_R2_001.fastq.gz')
  out_files <- paste0(filter_data_path,i,
                      c('_R1_paired.fastq.gz','_R1_unpaired.fastq.gz','_R2_paired.fastq.gz','_R2_unpaired.fastq.gz')
                      ,collapse = ' ')
  system(paste('sbatch /storage/htc/birchlerlab/CENH3_ChIP-seq/scripts/bash_sub.sbatch java -jar /cluster/software/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE',infile1,infile2,
               out_files,'LEADING:30 TRAILING:29 SLIDINGWINDOW:4:30 MINLEN:29 AVGQUAL:28'))
}

#2. QC with fastqc
for(i in list.files(filter_data_path,full.names = T)){
  system(paste0('sbatch /storage/htc/birchlerlab/CENH3_ChIP-seq/scripts/bash_sub.sbatch fastqc ',
                i,' -o ',qc_result_path))
}

#3. Map reads to genome (B73 + B cenromere sequences, should bind together) and merge bam files for replicates
for (i in all_samples) {
  infile1 <- paste0(filter_data_path,i,'_R1_paired.fastq.gz')
  infile2 <- paste0(filter_data_path,i,'_R2_paired.fastq.gz')
  infile3 <- paste0(filter_data_path,i,c('_R1_unpaired.fastq.gz','_R2_unpaired.fastq.gz'),collapse = ',')
  outfile <- paste0(sam_data_path,i,'.sam')
  metric_file <- paste0(align_log_path,i,'.log')
  system(paste('sbatch /storage/htc/birchlerlab/CENH3_ChIP-seq/scripts/bash_sub.sbatch bowtie2 -n 3 -x ',ref_genome_path,
               '-1',infile1,'-2',infile2,'-U',infile3,'-S',outfile,'--met-file',metric_file))
}
#4. sam to bam
for (i in all_samples) {
  infile <- paste0(sam_data_path,i,'.sam')
  outfile <- paste0(bam_data_path,i,'.bam')
  system(paste('sbatch /storage/htc/birchlerlab/CENH3_ChIP-seq/scripts/bash_sub.sbatch samtools view -bSF4',
               infile,'-o',outfile))
}
#5. merge bam
all_samples1 <- unique(gsub('_.+','',all_samples))
for (i in all_samples1) {
  infile <- paste0(list.files(bam_data_path,pattern = paste0(i,'.+'),full.names = T),collapse = ' ')
  outfile <- paste0(bam_merge_path,i,'.bam')
  system(paste('sbatch /storage/htc/birchlerlab/CENH3_ChIP-seq/scripts/bash_sub.sbatch samtools merge',
               outfile,infile))
}

#6. MACS2 call peaks 
for(i in list.files(bam_merge_path)){
  experiment_name <- gsub('.bam','',i)
  system(paste0('sbatch /storage/htc/birchlerlab/CENH3_ChIP-seq/scripts/bash_sub.sbatch macs2 callpeak -t ',
                bam_merge_path,i,' -n ',experiment_name,' --outdir ',peakcalling_results,experiment_name ,' -f BAM -g 1.2e8 -B -q 0.01 --nomodel'))
}










