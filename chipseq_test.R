system('module load bowtie2/bowtie2-2.2.9')
system('module load trimmomatic/trimmomatic-0.36')
system('module load fastqc')

raw_data_path = '/storage/htc/birchlerlab/CENH3_ChIP-seq/data/ChIP_data/'
trim_data_path = '/storage/htc/birchlerlab/CENH3_ChIP-seq/data/ChIP_data_trimmed/'
qc_result_path = '/storage/htc/birchlerlab/CENH3_ChIP-seq/data/ChIP_data_qc/'
sam_data_path = '/storage/htc/birchlerlab/CENH3_ChIP-seq/data/sams/'
bam_data_path = '/storage/htc/birchlerlab/CENH3_ChIP-seq/data/bams/'
bam_merge_path = '/storage/htc/birchlerlab/CENH3_ChIP-seq/data/bams_merge/'

log_path = '/storage/htc/birchlerlab/CENH3_ChIP-seq/logs/'

ref_genome_path = '/storage/htc/bdm/ref_genomes/B73/index_bowtie2'
dir.create(trim_data_path)
dir.create(qc_result_path)
dir.create(bam_data_path)
dir.create(bam_merge_path)
dir.create(log_path)

all_samples <- unique(gsub('_R.+','',list.files(raw_data_path)))
setwd(log_path)

#1.Trimming adapter and Filter low quality reads(Hua not run)
for(i in all_samples){
  infile1 <- paste0(raw_data_path,i,'_R1_001.fastq.gz')
  infile2 <- paste0(raw_data_path,i,'_R2_001.fastq.gz')
  out_files <- paste0(trim_data_path,i,
                      c('_R1_paired.fastq.gz','_R1_unpaired.fastq.gz','_R2_paired.fastq.gz','_R2_unpaired.fastq.gz')
                      ,collapse = ' ')
  system(paste('sbatch /storage/htc/birchlerlab/CENH3_ChIP-seq/scripts/bash_sub.sbatch java -jar /cluster/software/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE',infile1,infile2,
               out_files,'LEADING:30 TRAILING:29 SLIDINGWINDOW:4:30 MINLEN:29 AVGQUAL:28'))
}
#2. QC with fastqc
for(i in list.files(trim_data_path)){
  system(paste0('sbatch /storage/htc/birchlerlab/CENH3_ChIP-seq/scripts/bash_sub.sbatch fastqc ',
                trim_data_path,i,' -o ',qc_result_path))
}

#3. Map reads to genome (B73 + B cenromere sequences, should bind together) and merge bam files for replicates
for (i in all_samples) {
  infile1 <- paste0(trim_data_path,i,'_R1_paired.fastq.gz')
  infile2 <- paste0(trim_data_path,i,'_R2_paired.fastq.gz')
  infile3 <- paste0(trim_data_path,i,c('_R1_unpaired.fastq.gz','_R2_unpaired.fastq.gz'),collapse = ',')
  outfile <- paste0(sam_data_path,i,'.sam')
  system(paste('sbatch /storage/htc/birchlerlab/CENH3_ChIP-seq/scripts/bash_sub.sbatch bowtie2 -x ',ref_genome_path,
               '-1',infile1,'-2',infile2,'-U',infile3,'-S',outfile))
}
# sam to bam
for (i in all_samples) {
  infile <- paste0(sam_data_path,i,'.sam')
  outfile <- paste0(bam_data_path,i,'.bam')
  print(paste('sbatch /storage/htc/birchlerlab/CENH3_ChIP-seq/scripts/bash_sub.sbatch samtools view -bSF4',
               infile,'>',outfile))
}
#merge bam
all_samples1 <- unique(gsub('_.+','',all_samples))
for (i in all_samples1) {
  infile <- paste0(list.files(bam_data_path,pattern = paste0(i,'.+'),full.names = T),collapse = ' ')
  outfile <- paste0(bam_merge_path,i,'.bam')
  system(paste('sbatch /storage/htc/birchlerlab/CENH3_ChIP-seq/scripts/bash_sub.sbatch samtools merge',
               outfile,infile))
}

#4. MACS2 call peaks 
