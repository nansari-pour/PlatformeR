#' PlatformeR main function for generating per chromosome BAMs (with CNA or diploid)
#'
#' The function creates a BAM for a chromosome with >=1 CNA simulated in it
#' and creates a diploid germline normal for the same chromosome if required (i.e.paired analysis)
#' @param cna_simulate_file Full path to the file with CNA simulations which is a table with five columns (chr, startpos, endpos, cna, CCF) with no headers
#' @param chrom Name of the chromosome to be simulated
#' @param out_dir Full path to the out directory where genomewide files will be written
#' @param fasta_dir Full path to the chromosome fasta directory
#' @param phase_dir Full path to the phase reference file directory
#' @param art_bin Full path to the ART bin tool for Illumina read simulation
#' @param bwa Full path to the BWA tool for read alignment
#' @param samtools Full path to the samtools bin
#' @param read_length Length of the Illumina reads to be simulated (integer, usually 150)
#' @param fragment_size Mean size of the sequencing library fragments to be simulated (integer)
#' @param fragment_size_sd Standard deviation of the size of the sequencing library fragments to be simulated (integer)
#' @param tmp_dir Name of the temporary directory for samtools run (character string, e.g. "TMP")
#' @param ncores Number of cores to be used in running BWA and samtools
#' @param coverage Sequencing coverage depth to be simulated per chromosome copy for the tumour (numeric)
#' @param normal_coverage Sequencing coverage depth to be simulated per chromosome copy for the diploid germline/control sample (numeric)
#' @param simulated_purity The tumour purity to be simulated (numeric, ranging (0,1])
#' @param chromX_haplotype The haplotype onto which CNA should be simulated for chromosome X (string, default 'paternal')
#' @param loh_haplotype The haplotype onto which LOH should be simulated (string, either 'maternal' or 'paternal')
#' @param gain_haplotype The haplotype onto which GAIN should be simulated (string, either 'maternal' or 'paternal')
#' @param generate_normal set to FALSE if a normal germline diploid BAM for the same chromosome is not required e.g. for tumour-only simulations (default=TRUE)
#' @param generate_diploid_normal set to FALSE if a normal germline diploid BAM is not required for a non-CNA chromosome in the tumour e.g. for tumour-only simulations (default=TRUE)
#' @author naser.ansari-pour
#' @export


PlatformeR <- function(cna_simulate_file,chrom,out_dir,fasta_dir,phase_dir,art_bin,bwa,samtools,read_length,fragment_size,fragment_size_sd,tmp_dir,ncores,coverage,normal_coverage,simulated_purity,chromX_haplotype="paternal",loh_haplotype,gain_haplotype,generate_normal=TRUE,generate_diploid_normal=TRUE){
  # Read in CNA simulation input for chromosomes
  cna_simulate=read.table(cna_simulate_file, header = F, stringsAsFactors = F)
  names(cna_simulate)=c("chr","startpos","endpos","cna","CCF","TCN")
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

    if (chrom=="X"){
      
      run_chromX_cna(cna_simulate = cna_simulate,
                     fasta_dir = fasta_dir,
                     phase_dir = phase_dir,
                     art_bin = ART,
                     haploid_coverage = coverage/2,
                     read_length = read_length,
                     fragment_size = fragment_size,
                     fragment_size_sd = fragment_size_sd,
                     tmp_dir=tmp_dir,
                     bwa = BWA,
                     ncores = ncores,
                     samtools = SAM,
                     cna_haplotype = chromX_haplotype,
                     generate_normal = generate_normal,
                     haploid_coverage_normal = normal_coverage/2,
                     simulated_purity = simulated_purity)
    } else {

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
                  generate_normal = generate_normal,
                  haploid_coverage_normal = normal_coverage/2,
                  simulated_purity = simulated_purity)

  } 
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
                         haploid_cov = ifelse(chrom=="X",coverage/4,coverage/2)
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
                           haploid_cov = ifelse(chrom=="X",normal_coverage/4,normal_coverage/2)
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
