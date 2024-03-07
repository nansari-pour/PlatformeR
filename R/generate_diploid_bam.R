#' Create Normal diploid BAM for a chromosome
#'
#' The function creates a diploid BAM for a given chromosome mainly for copy number paired-analysis 
#' as germline or to generate a non-CNA chromosome for the tumour 
#' @param chrom Name of the chromosome to be simulated 
#' @param fasta_dir Full path to the chromosome fasta directory
#' @param phase_dir Full path to the phase reference file directory
#' @param art_bin Full path to the ART bin tool for Illumina read simulation
#' @param haploid_coverage Sequencing coverage depth to be simulated per chromosome copy (numeric)
#' @param read_length Length of the Illumina reads to be simulated (integer, usually 150)
#' @param fragment_size Mean size of the sequencing library fragments to be simulated (integer)
#' @param fragment_size_sd Standard deviation of the size of the sequencing library fragments to be simulated (integer)
#' @param tmp_dir Name of the temporary directory for samtools run (character string, e.g. "TMP")
#' @param bwa Full path to the BWA tool for read alignment
#' @param refseq Full path to the fasta reference file to be used by BWA 
#' @param bam Name of the BAM file to be generated (string)
#' @param ncores Number of cores to be used in running BWA and samtools
#' @param samtools Full path to the samtools bin
#' @param logfile Name of the logfile for the BWA run (string)
#' @param skip_art set to TRUE if ART run is already complete - files are expected in the working directory (default = FALSE)
#' @param skip_bwa set to TRUE if BWA run is already complete - files are expected in the working directory (default = FALSE)
#' @author naser.ansari-pour
#' @export

generate_diploid_bam <- function(chrom,fasta_dir,phase_dir,art_bin,haploid_coverage,read_length,fragment_size,fragment_size_sd,tmp_dir,bwa,refseq,bam,ncores,samtools,logfile,skip_art=FALSE,skip_bwa=FALSE){
  
  print(paste0("Initiating simulation of chromosome ",chrom," - diploid"))
  
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
                     fq1 = paste0(gsub(".fa","_",normal_fasta_filename),"1.fq"),
                     fq2 = paste0(gsub(".fa","_",normal_fasta_filename),"2.fq"),
                     bwa = bwa,
                     refseq = refseq,
                     bam = bam,
                     ncores = ncores,
                     samtools = samtools,
                     logfile = logfile,
                     skip_art = skip_art,
                     skip_bwa = skip_bwa)
}
