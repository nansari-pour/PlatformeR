#' @param out_dir Full path to the out directory where genomewide files will be written
#' @param fasta_dir Full path to the chromosome fasta directory
#' @param phase_dir Full path to the phase reference file directory
#' @param art_bin Full path to the ART bin tool for Illumina read simulation
#' @param bwa Full path to the BWA tool for read alignment
#' @param samtools Full path to the samtools bin
#' @param tmp_dir Name of the temporary directory for samtools run (character string, e.g. "TMP")
#' @param ncores Number of cores to be used in running BWA and samtools
#' @param coverage Sequencing coverage depth to be simulated per chromosome copy for the tumour (numeric)
#' @param read_length Length of the Illumina reads to be simulated (integer, usually 150)
#' @param fragment_size Mean size of the sequencing library fragments to be simulated (integer)
#' @param fragment_size_sd Standard deviation of the size of the sequencing library fragments to be simulated (integer)
#' @param simulated_purity The tumour purity to be simulated (numeric, ranging (0,1])
#' @param loh_haplotype The haplotype onto which LOH should be simulated (string, either 'maternal' or 'paternal')
#' @param gain_haplotype The haplotype onto which GAIN should be simulated (string, either 'maternal' or 'paternal')
#' @param generate_normal set to FALSE if a normal germline diploid BAM for the same chromosome is not required e.g. for tumour-only simulations (default=TRUE)
#' @param normal_coverage Sequencing coverage depth to be simulated per chromosome copy for the diploid germline/control sample (numeric)
#' @param generate_diploid_normal set to FALSE if a normal germline diploid BAM is not required for a non-CNA chromosome in the tumour e.g. for tumour-only simulations (default=TRUE)
#' @author naser.ansari-pour
#' @export


PlatformeR <- function(cna_simulate_file,chrom,out_dir,fasta_dir,phase_dir,art_bin,coverage,read_length,fragment_size,fragment_size_sd,tmp_dir,bwa,ncores,samtools,loh_haplotype,gain_haplotype,simulated_purity,generate_normal,normal_coverage,generate_diploid_normal=TRUE){
  # Read in CNA simulation input for chromosomes
  cna_simulate=read.table(cna_simulate_file, header = F, stringsAsFactors = F)
  names(cna_simulate)=c("chr","startpos","endpos","cna","CCF")
  print(cna_simulate)
  # Standardise Chr notation
  cna_simulate$chr=gsub("chr","",cna_simulate$chr)
  cna_chroms=unique(cna_simulate$chr)

  if (chrom %in% cna_chroms){
    print(paste("Chromosomes to be simulated with CNAs =",chrom))
    chrom_dir=paste0(out_dir,"chr",chrom,"/")
    dir.create(chrom_dir,recursive = TRUE,showWarnings = FALSE)
    setwd(chrom_dir)
    print(getwd())

    # run CNA simulation per chromosome

    run_chrom_cna(cna_simulate = cna_simulate,
                  chrom = chrom,
                  fasta_dir = fasta_dir,
                  phase_dir = phase_dir,
                  art_bin = art_bin,
                  haploid_coverage = coverage/2,
                  read_length = read_length,
                  fragment_size = fragment_size,
                  fragment_size_sd = fragment_size_sd,
                  tmp_dir=tmp_dir,
                  bwa = bwa,
                  ncores = ncores,
                  samtools = samtools,
                  loh_haplotype = loh_haplotype,
                  gain_haplotype = gain_haplotype,
                  simulated_purity = simulated_purity,
                  generate_normal = generate_normal,
                  haploid_coverage_normal = normal_coverage/2)

  } else {

    print(paste("Starting Chromosome ",chrom,"- no CNA"))
    chrom_dir=paste0(out_dir,"chr",chrom,"/")
    dir.create(chrom_dir,recursive = TRUE,showWarnings = FALSE)
    setwd(chrom_dir)
    print(getwd())

    generate_diploid_bam(chrom = chrom,
                         fasta_dir = fasta_dir,
                         phase_dir = phase_dir,
                         art_bin = art_bin,
                         haploid_cov = coverage/2,
                         read_length = read_length,
                         fragment_size = fragment_size,
                         fragment_size_sd = fragment_size_sd,
                         tmp_dir = tmp_dir,
                         bwa = bwa,
                         refseq = paste0(fasta_dir,"chr",chrom,".fa"),
                         bam = ifelse(simulated_purity==1,paste0("chr",chrom,".clone.bam"),paste0("chr",chrom,".clone_purity_",simulated_purity,".bam")),
                         ncores = ncores,
                         samtools = samtools,
                         logfile = "BWA_clone.log")

    if (generate_diploid_normal){

      generate_diploid_bam(chrom = chrom,
                           fasta_dir = fasta_dir,
                           phase_dir = phase_dir,
                           art_bin = art_bin,
                           haploid_cov = normal_coverage/2,
                           read_length = read_length,
                           fragment_size = fragment_size,
                           fragment_size_sd = fragment_size_sd,
                           tmp_dir = tmp_dir,
                           bwa = bwa,
                           refseq = paste0(fasta_dir,"chr",chrom,".fa"),
                           bam = paste0("chr",chrom,".normal.bam"),
                           ncores = ncores,
                           samtools = samtools,
                           logfile = "BWA_normal.log")
    }

  }
}