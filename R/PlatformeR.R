# PlatformeR main function

PlatformeR <- function(cna_simulate_file,chrom,out_dir,fasta_dir,phase_dir,art_bin,bwa,samtools,tmp_dir,ncores,coverage,normal_coverage,simulated_purity,loh_haplotype,gain_haplotype,generate_diploid_normal=TRUE){
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
                  art_bin = ART,
                  haploid_coverage = coverage/2,
                  read_length = 150,
                  fragment_size = 400,
                  fragment_size_sd = 10,
                  tmp_dir=tmp_dir,
                  bwa = BWA,
                  ncores = ncores,
                  samtools = SAM,
                  loh_haplotype = loh_haplotype,
                  gain_haplotype = gain_haplotype,
                  generate_normal = TRUE,
                  haploid_coverage_normal = normal_coverage/2)
    
  } else {
    
    print(paste("Starting Chromosome ",chrom,"- normal only"))
    chrom_dir=paste0(out_dir,"chr",chrom,"/")
    dir.create(chrom_dir,recursive = TRUE,showWarnings = FALSE)
    setwd(chrom_dir)
    print(getwd())
    
    generate_diploid_bam(chrom = chrom,
                         fasta_dir = fasta_dir,
                         phase_dir = phase_dir,
                         art_bin = art_bin,
                         haploid_cov = coverage/2,
                         read_length = 150,
                         fragment_size = 400,
                         fragment_size_sd = 10,
                         tmp_dir = tmp_dir,
                         bwa = BWA,
                         refseq = paste0(fasta_dir,"chr",chrom,".fa"),
                         bam = ifelse(simulated_purity==1,paste0("chr",chrom,".clone.bam"),paste0("chr",chrom,".clone_purity_",simulated_purity,".bam")),
                         ncores = ncores,
                         samtools = SAM,
                         logfile = "BWA_clone.log")
    
    if (generate_diploid_normal){
      
      generate_diploid_bam(chrom = chrom,
                           fasta_dir = fasta_dir,
                           phase_dir = phase_dir,
                           art_bin = art_bin,
                           haploid_cov = normal_coverage/2,
                           read_length = 150,
                           fragment_size = 400,
                           fragment_size_sd = 10,
                           tmp_dir = tmp_dir,
                           bwa = BWA,
                           refseq = paste0(fasta_dir,"chr",chrom,".fa"),
                           bam = paste0("chr",chrom,".normal.bam"),
                           ncores = ncores,
                           samtools = SAM,
                           logfile = "BWA_normal.log")
    }
    
  }
}
