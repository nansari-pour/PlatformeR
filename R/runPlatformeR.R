# run PlatformeR per chromosome

runPlatformeR <- function(cna_simulate,chrom,fasta_dir,phase_dir,art_bin,haploid_coverage,read_length,fragment_size,fragment_size_sd,tmp_dir,bwa,ncores,samtools,loh_haplotype,gain_haplotype,generate_normal=TRUE,haploid_coverage_normal){
  
  chrom_simulate=cna_simulate[cna_simulate$chr==chrom,]
  chrom_simulate$length=chrom_simulate$endpos-chrom_simulate$startpos
  print(chrom_simulate)
  
  # STEP 1 - generate maternal and paternal chromosome fasta file
  prepare_chrom_fasta(chrom=chrom,
                      fasta_dir = fasta_dir,
                      phase_dir = phase_dir)
  
  
  # STEP 2 - identify diploid regions and generate their bams if any
  
  diploid_regions=alternate_segments(chr_length,chrom_simulate[,c("startpos","endpos")])
  print(diploid_regions)
  
  if (!is.null(diploid_regions)){
  if (nrow(diploid_regions)>0){
    
    for (region in 1:nrow(diploid_regions)){
      
      simulate_normal_region(chrom = chrom,
                             region = region,
                             diploid_regions = diploid_regions,
                             art_bin = art_bin,
                             haploid_cov = haploid_coverage,
                             read_length = read_length,
                             fragment_size = fragment_size,
                             fragment_size_sd = fragment_size_sd,
                             tmp_dir = tmp_dir,
                             bwa = bwa,
                             refseq = paste0(fasta_dir,"chr",chrom,".fa"),
                             bam = paste0("chr",chrom,".normal.clone_region",region,".bam"),
                             ncores = ncores,
                             samtools = samtools,
                             logfile = paste0("BWA_normal_clone_region",region,".log"))
    }
  } else {print("no diploid regions to be simulated")}
  } else {print("CNA across entire chromosome")}
  
  # STEP 3 - identify and generate BAMs for LOH regions if present in chrom_simulate
  
  seeds=1:100
  
  loh_simulate=chrom_simulate[which(!is.na(match(chrom_simulate$cna,c("LOH","loh","Loh")))),]
  
  if (nrow(loh_simulate)>0){
    
    print(loh_simulate)
    
    
    for (loh_region in 1:nrow(loh_simulate)){
      
      loh_ccf=loh_simulate$CCF[loh_region]
      
      simulate_loh_region(chrom = chrom,
                          loh_haplotype = loh_haplotype,
                          loh_region = loh_region,
                          loh_simulate = loh_simulate,
                          art_bin = art_bin,
                          haploid_cov = haploid_coverage,
                          read_length = read_length,
                          fragment_size = fragment_size,
                          fragment_size_sd = fragment_size_sd,
                          tmp_dir = tmp_dir,
                          bwa = bwa,
                          refseq = paste0(fasta_dir,"chr",chrom,".fa"),
                          bam = paste0("chr",chrom,".LOH.clone_region",loh_region,".bam"),
                          ncores = ncores,
                          samtools = samtools,
                          logfile = paste0("BWA_LOH_clone_region",loh_region,".log"))
      
      if (loh_ccf*simulated_purity<1){
        
        LOH_normal_fasta_filename=paste0("chr",chrom,".LOH.normal_region",loh_region,".fa")
        
        LOH_normal_paternal_fasta=substr(chrom_fasta_ref,loh_simulate$startpos[loh_region],loh_simulate$endpos[loh_region])
        writeLines(paste0(">chr",chrom,"-paternal","-LOH-",loh_region,"_normal"),LOH_normal_fasta_filename)
        cat(LOH_normal_paternal_fasta, sep = "\n", file = LOH_normal_fasta_filename, append = TRUE)
        LOH_normal_maternal_fasta=substr(chrom_fasta_alt,loh_simulate$startpos[loh_region],loh_simulate$endpos[loh_region])
        cat(paste0(">chr",chrom,"-maternal","-LOH-",loh_region,"_normal"),sep = "\n", file = LOH_normal_fasta_filename, append = TRUE)
        cat(LOH_normal_maternal_fasta, sep = "\n", file = LOH_normal_fasta_filename, append = TRUE)
        
        fasta2bam_simulate(clone_fasta_filename = LOH_normal_fasta_filename,
                           art_bin = art_bin,
                           haploid_cov = haploid_coverage,
                           read_length = read_length,
                           fragment_size = fragment_size,
                           fragment_size_sd = fragment_size_sd,
                           tmp_dir = tmp_dir,
                           fq1 = paste0(gsub(".fa","_",LOH_normal_fasta_filename),"1.fq"),
                           fq2 = paste0(gsub(".fa","_",LOH_normal_fasta_filename),"2.fq"),
                           bwa = bwa,
                           refseq = paste0(fasta_dir,"chr",chrom,".fa"),
                           bam = paste0("chr",chrom,".LOH.normal_region",loh_region,".bam"),
                           ncores = ncores,
                           samtools = samtools,
                           logfile = paste0("BWA_LOH_normal_region",loh_region,".log"),
                           skip_art = FALSE,
                           skip_bwa = FALSE)
        
        down_seed=sample(seeds,1,replace = FALSE)
        
        mix_bams_loh(simulated_purity = simulated_purity,
                     cna_ccf = loh_ccf,
                     chrom_length = chr_length,
                     cna_length = loh_simulate$endpos[loh_region]-loh_simulate$startpos[loh_region],
                     downsampling_seed = down_seed,
                     samtools = samtools,
                     clone0 = paste0("chr",chrom,".LOH.clone_region",loh_region),
                     clone1 = paste0("chr",chrom,".LOH.normal_region",loh_region),
                     clone01 = ifelse (loh_ccf<1,paste0("chr",chrom,".sLOH.clone_region",loh_region),paste0("chr",chrom,".LOH.clone_region",loh_region))
        )
        
        seeds=setdiff(seeds,down_seed)
        if (loh_ccf<1){
          system(paste("rm",paste0("chr",chrom,".LOH.clone_region",loh_region,".bam*")))
        }
        system(paste("rm",paste0("chr",chrom,".LOH.*_region",loh_region,"*down.bam")))
        
      }
    }
  } else {print("no LOH to be simulated")}
  
  
  
  # STEP 4 - identify and generate BAMs for GAIN regions if present in chrom_simulate
  
  gain_simulate=chrom_simulate[which(!is.na(match(chrom_simulate$cna,c("GAIN","Gain","gain")))),]
  
  if (nrow(gain_simulate)>0){
    
    print(gain_simulate)
    
    for (gain_region in 1:nrow(gain_simulate)){
      
      gain_ccf=gain_simulate$CCF[gain_region]
      
      simulate_gain_region(chrom = chrom,
                           gain_haplotype = gain_haplotype,
                           gain_region = gain_region,
                           gain_simulate = gain_simulate,
                           art_bin = art_bin,
                           haploid_cov = haploid_coverage,
                           read_length = read_length,
                           fragment_size = fragment_size,
                           fragment_size_sd = fragment_size_sd,
                           tmp_dir = tmp_dir,
                           bwa = bwa,
                           refseq = paste0(fasta_dir,"chr",chrom,".fa"),
                           bam = paste0("chr",chrom,".GAIN.clone_region",gain_region,".bam"),
                           ncores = ncores,
                           samtools = samtools,
                           logfile = paste0("BWA_GAIN_clone_region",gain_region,".log"),
                           skip_art = FALSE,
                           skip_bwa = FALSE)
      
      
      if (gain_ccf*simulated_purity<1){
        
        GAIN_normal_fasta_filename=paste0("chr",chrom,".GAIN.normal_region",gain_region,".fa")
        
        GAIN_normal_paternal_fasta=substr(chrom_fasta_ref,gain_simulate$startpos[gain_region],gain_simulate$endpos[gain_region])
        writeLines(paste0(">chr",chrom,"-paternal","-GAIN-",gain_region,"_normal"),GAIN_normal_fasta_filename)
        cat(GAIN_normal_paternal_fasta, sep = "\n", file = GAIN_normal_fasta_filename, append = TRUE)
        GAIN_normal_maternal_fasta=substr(chrom_fasta_alt,gain_simulate$startpos[gain_region],gain_simulate$endpos[gain_region])
        cat(paste0(">chr",chrom,"-maternal","-GAIN-",gain_region,"_normal"),sep = "\n", file = GAIN_normal_fasta_filename, append = TRUE)
        cat(GAIN_normal_maternal_fasta, sep = "\n", file = GAIN_normal_fasta_filename, append = TRUE)
        
        fasta2bam_simulate(clone_fasta_filename = GAIN_normal_fasta_filename,
                           art_bin = art_bin,
                           haploid_cov = haploid_coverage,
                           read_length = read_length,
                           fragment_size = fragment_size,
                           fragment_size_sd = fragment_size_sd,
                           tmp_dir = tmp_dir,
                           fq1 = paste0(gsub(".fa","_",GAIN_normal_fasta_filename),"1.fq"),
                           fq2 = paste0(gsub(".fa","_",GAIN_normal_fasta_filename),"2.fq"),
                           bwa = bwa,
                           refseq = paste0(fasta_dir,"chr",chrom,".fa"),
                           bam = paste0("chr",chrom,".GAIN.normal_region",gain_region,".bam"),
                           ncores = ncores,
                           samtools = samtools,
                           logfile = paste0("BWA_GAIN_normal_region",gain_region,".log"),
                           skip_art = FALSE,
                           skip_bwa = FALSE)
        
        down_seed=sample(seeds,1,replace = FALSE)
        
        mix_bams_gain(simulated_purity = simulated_purity,
                      cna_ccf = gain_ccf,
                      chrom_length = chr_length,
                      cna_length = gain_simulate$endpos[gain_region]-gain_simulate$startpos[gain_region],
                      downsampling_seed = down_seed,
                      samtools = samtools,
                      clone0 = paste0("chr",chrom,".GAIN.clone_region",gain_region),
                      clone1 = paste0("chr",chrom,".GAIN.normal_region",gain_region),
                      clone01 = ifelse (gain_ccf<1,paste0("chr",chrom,".sGAIN.clone_region",gain_region),paste0("chr",chrom,".GAIN.clone_region",gain_region))
        )
        if (gain_ccf<1){
          system(paste("rm",paste0("chr",chrom,".GAIN.clone_region",gain_region,".bam*")))
        }
        system(paste("rm",paste0("chr",chrom,".GAIN.*_region",gain_region,"*down.bam")))
        
      }
    }
  } else {print("no GAIN to be simulated")}
  
  # STEP 5 - merge BAMs of all segments to generate a full chromosome containing CNAs in chrom_simulate
  
  system(paste("ls", paste0("chr",chrom,"*clone_region*bam >"),paste0("chr",chrom,"_mini_bams_list.txt")))
  
  merge_mini_bams(samtools = samtools,
                  mini_bams_list_file = paste0("chr",chrom,"_mini_bams_list.txt"),
                  ncores = ncores,
                  clone_bam_filename = ifelse(simulated_purity==1,paste0("chr",chrom,".clone.bam"),paste0("chr",chrom,".clone_purity_",simulated_purity,".bam")))
  
  
  # STEP 6 - OPTIONAL - Create Normal BAM if paired analysis
  
  if (generate_normal){
    
    generate_diploid_bam(chrom = chrom,
                        fasta_dir = fasta_dir,
                        phase_dir = phase_dir,
                        art_bin = art_bin,
                        haploid_cov = haploid_coverage_normal,
                        read_length = read_length,
                        fragment_size = fragment_size,
                        fragment_size_sd = fragment_size_sd,
                        tmp_dir = tmp_dir,
                        bwa = bwa,
                        refseq = paste0(fasta_dir,"chr",chrom,".fa"),
                        bam = paste0("chr",chrom,".normal.bam"),
                        ncores = ncores,
                        samtools = samtools,
                        logfile = "BWA_normal_clone.log",
                        skip_art = FALSE,
                        skip_bwa = FALSE)
  }
}
