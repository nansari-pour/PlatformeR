# Simulate normal region

simulate_normal_region <- function(chrom,region,diploid_regions,clone_fasta_filename,art_bin,haploid_cov,read_length,fragment_size,fragment_size_sd,tmp_dir,fq1,fq2,bwa,refseq,bam,ncores,samtools,logfile,skip_art=FALSE,skip_bwa=FALSE){
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

fasta2bam_simulate(clone_fasta_filename = clone_fasta_filename,
                   art_bin = art_bin,
                   haploid_cov = haploid_cov,
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
                   logfile = logfile)
}

# Simulate LOH region

simulate_loh_region <- function(chrom,loh_haplotype,loh_region,clone_fasta_filename,art_bin,haploid_cov,read_length,fragment_size,fragment_size_sd,tmp_dir,fq1,fq2,bwa,refseq,bam,ncores,samtools,logfile,skip_art=FALSE,skip_bwa=FALSE){
  
  loh_ccf=LOH_simulate$CCF[loh_region]
  
  if (loh_haplotype=="maternal"){
    # LOH_region_maternal_fasta is NULL - because of LOH, so only generate LOH_region paternal fasta
    LOH_region_clone_fasta_filename=paste0("chr",chrom,".LOH_",loh_region,".fa")
    LOH_region_clone_fasta=substr(chrom_fasta_ref,LOH_simulate$startpos[loh_region],LOH_simulate$endpos[loh_region])
    writeLines(paste0(">chr",chrom,"-paternal","-LOH-",loh_region),LOH_region_clone_fasta_filename)
    cat(LOH_region_clone_fasta, sep = "\n", file = LOH_region_clone_fasta_filename, append = TRUE)
    print(paste0(">chr",chrom,"-paternal","-LOH-",loh_region," length = ",nchar(LOH_region_clone_fasta)))
  } else if (loh_haplotype=="paternal"){
    # LOH_region_paternal_fasta is NULL - because of LOH, so only generate LOH_region maternal fasta 
    LOH_region_clone_fasta_filename=paste0("chr",chrom,".LOH_",loh_region,".fa")
    LOH_region_clone_fasta=substr(chrom_fasta_alt,LOH_simulate$startpos[loh_region],LOH_simulate$endpos[loh_region])
    writeLines(paste0(">chr",chrom,"-maternal","-LOH-",loh_region),LOH_region_clone_fasta_filename)
    cat(LOH_region_clone_fasta, sep = "\n", file = LOH_region_clone_fasta_filename, append = TRUE)
    print(paste0(">chr",chrom,"-maternal","-LOH-",loh_region," length = ",nchar(LOH_region_clone_fasta)))
  } else {
    stop("Wrong LOH haplotype designation")
  }
  
  fasta2bam_simulate(clone_fasta_filename = clone_fasta_filename,
                     art_bin = art_bin,
                     haploid_cov = haploid_cov,
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
                     logfile = logfile)
  
}

# Simulate GAIN region

simulate_gain_region <- function(chrom,gain_haplotype,gain_region,clone_fasta_filename,art_bin,haploid_cov,read_length,fragment_size,fragment_size_sd,tmp_dir,fq1,fq2,bwa,refseq,bam,ncores,samtools,logfile,skip_art=FALSE,skip_bwa=FALSE){
  
  gain_ccf=GAIN_simulate$CCF[gain_region]
  
  GAIN_region_clone_fasta_filename=paste0("chr",chrom,".GAIN_",gain_region,".fa")
  
  GAIN_region_maternal_fasta_filename=paste0("chr",chrom,".maternal_GAIN_",gain_region,".fa")
  GAIN_region_maternal_fasta=substr(chrom_fasta_alt,GAIN_simulate$startpos[gain_region],GAIN_simulate$endpos[gain_region])
  writeLines(paste0(">chr",chrom,"-maternal","-GAIN-",gain_region),GAIN_region_maternal_fasta_filename)
  cat(GAIN_region_maternal_fasta, sep = "\n", file = GAIN_region_maternal_fasta_filename, append = TRUE)
  print(paste0(">chr",chrom,"-maternal","-GAIN-",gain_region," length = ",nchar(GAIN_region_maternal_fasta)))
  
  GAIN_region_paternal_fasta_filename=paste0("chr",chrom,".paternal_GAIN_",gain_region,".fa")
  GAIN_region_paternal_fasta=substr(chrom_fasta_ref,GAIN_simulate$startpos[gain_region],GAIN_simulate$endpos[gain_region])
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
  
  fasta2bam_simulate(clone_fasta_filename = clone_fasta_filename,
                     art_bin = art_bin,
                     haploid_cov = haploid_cov,
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
                     logfile = logfile)
}
