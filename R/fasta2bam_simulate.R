# Function to simulate BAMs from a fasta sequence input (i.e. has 1 to N fasta sequences in the the fasta file)

fasta2bam_simulate <- function(clone_fasta_filename,art_bin,haploid_cov,read_length,fragment_size,fragment_size_sd,tmp_dir,fq1,fq2,bwa,refseq,bam,ncores,samtools,logfile,skip_art=FALSE,skip_bwa=FALSE){
  clone_fastq_filename=gsub(".fa","_",clone_fasta_filename)
  if (!skip_art){
    system(paste0(art_bin," -i ",clone_fasta_filename," -p -ss HS25 -f ",haploid_cov," -na -l ",read_length," -m ",fragment_size," -s ",fragment_size_sd," -o ",clone_fastq_filename), wait = TRUE)
  }
  if (!skip_bwa){
    system(paste("mkdir -p",tmp_dir))
    print(paste(fq1,fq2,bam))
    system(paste(bwa,"mem",refseq,fq1,fq2,"-t",ncores,"|",samtools,"sort - -O bam -o",bam,"-T",tmp_dir,"-@",ncores,">",logfile,"2>&1"), wait=TRUE)
    system(paste(samtools,"index",bam), wait = TRUE)
  }
}
