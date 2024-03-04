#' Function to simulate normal region BAM
#' 
#' This function simulates a normal bam for a given region of interest on a chromosome
#' @param chrom Name of the chromosome to be simulated (chrom number without 'chr' e.g. 1,2,...,22)
#' @param region The region number to be simulated among the diploid_regions (integer)
#' @param diploid_regions A dataframe with two columns containing startpos and endpos of each diploid region
#' @param art_bin Full path to the ART bin tool for Illumina read simulation
#' @param haploid_cov Sequencing coverage depth to be simulated per chromosome copy (numeric)
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

simulate_normal_region <- function(chrom,region,diploid_regions,art_bin,haploid_cov,read_length,fragment_size,fragment_size_sd,tmp_dir,bwa,refseq,bam,ncores,samtools,logfile,skip_art=FALSE,skip_bwa=FALSE){
region_maternal_fasta_filename=paste0("chr",chrom,".maternal_",region,".fa")
region_maternal_fasta=substr(chrom_fasta_alt,diploid_regions$startpos[region],diploid_regions$endpos[region])
writeLines(paste0(">chr",chrom,"-maternal","-",region),region_maternal_fasta_filename)
cat(region_maternal_fasta, sep = "\n", file = region_maternal_fasta_filename, append = TRUE)
print(paste0(">chr",chrom,"-maternal","-",region," length = ",nchar(region_maternal_fasta)))

region_paternal_fasta_filename=paste0("chr",chrom,".paternal_",region,".fa")
region_paternal_fasta=substr(chrom_fasta_ref,diploid_regions$startpos[region],diploid_regions$endpos[region])
writeLines(paste0(">chr",chrom,"-paternal","-",region),region_paternal_fasta_filename)
cat(region_paternal_fasta, sep = "\n", file = region_paternal_fasta_filename, append = TRUE)
print(paste0(">chr",chrom,"-paternal","-",region," length = ",nchar(region_paternal_fasta)))

region_clone_fasta_filename=paste0("chr",chrom,".normal.clone_region",region,".fa")

system(paste("cat",region_maternal_fasta_filename,region_paternal_fasta_filename,">",region_clone_fasta_filename), wait = TRUE)
system(paste("rm",region_maternal_fasta_filename,region_paternal_fasta_filename))

fasta2bam_simulate(clone_fasta_filename = region_clone_fasta_filename,
                   art_bin = art_bin,
                   haploid_cov = haploid_cov,
                   read_length = read_length,
                   fragment_size = fragment_size,
                   fragment_size_sd = fragment_size_sd,
                   tmp_dir = tmp_dir,
                   fq1 = paste0(gsub(".fa","_",region_clone_fasta_filename),"1.fq"),
                   fq2 = paste0(gsub(".fa","_",region_clone_fasta_filename),"2.fq"),
                   bwa = bwa,
                   refseq = refseq,
                   bam = bam,
                   ncores = ncores,
                   samtools = samtools,
                   logfile = logfile,
                   skip_art = skip_art,
                   skip_bwa = skip_bwa)
}

#' Function to simulate LOH region BAM
#' 
#' This function simulates an LOH bam for a given region of interest on a chromosome
#' @param chrom Name of the chromosome to be simulated (chrom number without 'chr' e.g. 1,2,...,22)
#' @param loh_haplotype The haplotype onto which LOH should be simulated (string, either 'maternal' or 'paternal')
#' @param loh_region The region number to be simulated among the loh_simulate regions (integer)
#' @param loh_simulate A dataframe with two columns containing startpos and endpos of each loh region
#' @param art_bin Full path to the ART bin tool for Illumina read simulation
#' @param haploid_cov Sequencing coverage depth to be simulated per chromosome copy (numeric)
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

simulate_loh_region <- function(chrom,loh_haplotype,loh_region,loh_simulate,art_bin,haploid_cov,read_length,fragment_size,fragment_size_sd,tmp_dir,bwa,refseq,bam,ncores,samtools,logfile,skip_art=FALSE,skip_bwa=FALSE){
  
  if (loh_haplotype=="maternal"){
    # LOH_region_maternal_fasta is NULL - because of LOH, so only generate LOH_region paternal fasta
    LOH_region_clone_fasta_filename=paste0("chr",chrom,".LOH_",loh_region,".fa")
    LOH_region_clone_fasta=substr(chrom_fasta_ref,loh_simulate$startpos[loh_region],loh_simulate$endpos[loh_region])
    writeLines(paste0(">chr",chrom,"-paternal","-LOH-",loh_region),LOH_region_clone_fasta_filename)
    cat(LOH_region_clone_fasta, sep = "\n", file = LOH_region_clone_fasta_filename, append = TRUE)
    print(paste0(">chr",chrom,"-paternal","-LOH-",loh_region," length = ",nchar(LOH_region_clone_fasta)))
  } else if (loh_haplotype=="paternal"){
    # LOH_region_paternal_fasta is NULL - because of LOH, so only generate LOH_region maternal fasta 
    LOH_region_clone_fasta_filename=paste0("chr",chrom,".LOH_",loh_region,".fa")
    LOH_region_clone_fasta=substr(chrom_fasta_alt,loh_simulate$startpos[loh_region],loh_simulate$endpos[loh_region])
    writeLines(paste0(">chr",chrom,"-maternal","-LOH-",loh_region),LOH_region_clone_fasta_filename)
    cat(LOH_region_clone_fasta, sep = "\n", file = LOH_region_clone_fasta_filename, append = TRUE)
    print(paste0(">chr",chrom,"-maternal","-LOH-",loh_region," length = ",nchar(LOH_region_clone_fasta)))
  } else {
    stop("Wrong LOH haplotype designation")
  }
  
  fasta2bam_simulate(clone_fasta_filename = LOH_region_clone_fasta_filename,
                     art_bin = art_bin,
                     haploid_cov = haploid_cov,
                     read_length = read_length,
                     fragment_size = fragment_size,
                     fragment_size_sd = fragment_size_sd,
                     tmp_dir = tmp_dir,
                     fq1 = paste0(gsub(".fa","_",LOH_region_clone_fasta_filename),"1.fq"),
                     fq2 = paste0(gsub(".fa","_",LOH_region_clone_fasta_filename),"2.fq"),
                     bwa = bwa,
                     refseq = refseq,
                     bam = bam,
                     ncores = ncores,
                     samtools = samtools,
                     logfile = logfile,
                     skip_art = skip_art,
                     skip_bwa = skip_bwa)
  
}


#' Function to simulate GAIN region BAM
#' 
#' This function simulates a GAIN bam for a given region of interest on a chromosome
#' @param chrom Name of the chromosome to be simulated (chrom number without 'chr' e.g. 1,2,...,22)
#' @param gain_haplotype The haplotype onto which GAIN should be simulated (string, either 'maternal' or 'paternal')
#' @param gain_region The region number to be simulated among the gain_simulate regions (integer)
#' @param gain_simulate A dataframe with two columns containing startpos and endpos of each gain region
#' @param art_bin Full path to the ART bin tool for Illumina read simulation
#' @param haploid_cov Sequencing coverage depth to be simulated per chromosome copy (numeric)
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

simulate_gain_region <- function(chrom,gain_haplotype,gain_region,gain_simulate,art_bin,haploid_cov,read_length,fragment_size,fragment_size_sd,tmp_dir,bwa,refseq,bam,ncores,samtools,logfile,skip_art=FALSE,skip_bwa=FALSE){
  
  GAIN_region_clone_fasta_filename=paste0("chr",chrom,".GAIN_",gain_region,".fa")
  
  GAIN_region_maternal_fasta_filename=paste0("chr",chrom,".maternal_GAIN_",gain_region,".fa")
  GAIN_region_maternal_fasta=substr(chrom_fasta_alt,gain_simulate$startpos[gain_region],gain_simulate$endpos[gain_region])
  writeLines(paste0(">chr",chrom,"-maternal","-GAIN-",gain_region),GAIN_region_maternal_fasta_filename)
  cat(GAIN_region_maternal_fasta, sep = "\n", file = GAIN_region_maternal_fasta_filename, append = TRUE)
  print(paste0(">chr",chrom,"-maternal","-GAIN-",gain_region," length = ",nchar(GAIN_region_maternal_fasta)))
  
  GAIN_region_paternal_fasta_filename=paste0("chr",chrom,".paternal_GAIN_",gain_region,".fa")
  GAIN_region_paternal_fasta=substr(chrom_fasta_ref,gain_simulate$startpos[gain_region],gain_simulate$endpos[gain_region])
  writeLines(paste0(">chr",chrom,"-paternal","-GAIN-",gain_region),GAIN_region_paternal_fasta_filename)
  cat(GAIN_region_paternal_fasta, sep = "\n", file = GAIN_region_paternal_fasta_filename, append = TRUE)
  print(paste0(">chr",chrom,"-paternal","-GAIN-",gain_region," length = ",nchar(GAIN_region_paternal_fasta)))
  
  if (gain_haplotype=="paternal"){
    GAINED_region_paternal_fasta_filename=gsub("_GAIN_","_GAINED_",GAIN_region_paternal_fasta_filename)
    system(paste("cp",GAIN_region_paternal_fasta_filename,GAINED_region_paternal_fasta_filename),wait=TRUE)
    system(paste("sed -i 's/-GAIN-/-GAINED-/g'",GAINED_region_paternal_fasta_filename),wait=TRUE)
    system(paste("cat",GAIN_region_paternal_fasta_filename,GAIN_region_maternal_fasta_filename,GAINED_region_paternal_fasta_filename,">",GAIN_region_clone_fasta_filename), wait = TRUE)
    system(paste("rm",GAINED_region_paternal_fasta_filename),wait=TRUE)
  } else if (gain_haplotype=="maternal"){
    GAINED_region_maternal_fasta_filename=gsub("_GAIN_","_GAINED_",GAIN_region_maternal_fasta_filename)
    system(paste("cp",GAIN_region_maternal_fasta_filename,GAINED_region_maternal_fasta_filename),wait=TRUE)
    system(paste("sed -i 's/-GAIN-/-GAINED-/g'",GAINED_region_maternal_fasta_filename),wait=TRUE)
    system(paste("cat",GAIN_region_paternal_fasta_filename,GAIN_region_maternal_fasta_filename,GAINED_region_maternal_fasta_filename,">",GAIN_region_clone_fasta_filename), wait = TRUE)
    system(paste("rm",GAINED_region_maternal_fasta_filename),wait=TRUE)
  } else {stop("Wrong GAIN haplotype designation")}
  system(paste("rm",GAIN_region_paternal_fasta_filename,GAIN_region_maternal_fasta_filename),wait=TRUE)
  
  fasta2bam_simulate(clone_fasta_filename = GAIN_region_clone_fasta_filename,
                     art_bin = art_bin,
                     haploid_cov = haploid_cov,
                     read_length = read_length,
                     fragment_size = fragment_size,
                     fragment_size_sd = fragment_size_sd,
                     tmp_dir = tmp_dir,
                     fq1 = paste0(gsub(".fa","_",GAIN_region_clone_fasta_filename),"1.fq"),
                     fq2 = paste0(gsub(".fa","_",GAIN_region_clone_fasta_filename),"2.fq"),
                     bwa = bwa,
                     refseq = refseq,
                     bam = bam,
                     ncores = ncores,
                     samtools = samtools,
                     logfile = logfile,
                     skip_art = skip_art,
                     skip_bwa = skip_bwa)
}
