# Create Normal diploid BAM for a chromosome (mainly for copy number paired-analysis)

generate_normal_bam <- function(chrom,fasta_dir,phase_dir,art_bin,haploid_coverage,read_length,fragment_size,fragment_size_sd,tmp_dir,fq1,fq2,bwa,refseq,bam,ncores,samtools,logfile,skip_art=FALSE,skip_bwa=FALSE){
  
  print(paste0("Initiating simulation of chromosome",chrom," - diploid"))
  
  prepare_chrom_fasta(chrom=chrom,
                      fasta_dir = fasta_dir,
                      phase_dir = phase_dir)
  
  paternal_fasta_filename=paste0("chr",chrom,".paternal.fa")
  maternal_fasta_filename=paste0("chr",chrom,".maternal.fa")
  normal_fasta_filename=paste0("chr",chrom,".normal.fa")
  system(paste("cat",maternal_fasta_filename,paternal_fasta_filename,">",normal_fasta_filename), wait = TRUE)
  
  fasta2bam_simulate(clone_fasta_filename = normal_fasta_filename,
                     art_bin = ART,
                     haploid_cov = haploid_coverage,
                     read_length = read_length,
                     fragment_size = fragment_size,
                     fragment_size_sd = fragment_size_sd,
                     tmp_dir = tmp_dir,
                     fq1 = fq1,
                     fq2 = fq2,
                     bwa = bwa,
                     refseq = refseq,
                     bam = bam,
                     ncores = ncores,
                     samtools = samtools,
                     logfile = logfile,
                     skip_art = FALSE,
                     skip_bwa = FALSE)
}
