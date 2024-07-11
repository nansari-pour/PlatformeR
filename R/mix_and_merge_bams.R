#' Function for mixing BAMs for subclonal LOH events
#' 
#' @param simulated_purity Full path to the chromosome fasta directory
#' @param cna_ccf The tumour purity to be simulated (numeric, ranging (0,1])
#' @param downsampling_seed The seed value to be used for downsampling with samtools
#' @param samtools Full path to the samtools bin
#' @param clone0 Name of the clone with the simulated LOH
#' @param clone1 Name of the clone with normal diploid copy number
#' @param clone01 Name of the mixed clone to be used for the output BAM
#' @author naser.ansari-pour
#' @export

mix_bams_loh <- function(simulated_purity,cna_ccf,downsampling_seed,samtools,clone0,clone1,clone01){
  print(paste("LOH CCF =",cna_ccf))
  print(paste("NORMAL CCF =",1-cna_ccf))
  loh_proportion=round(simulated_purity*cna_ccf,digits = 2)
  loh_prop_input=paste0(downsampling_seed,".",ifelse(loh_proportion>=0.1,loh_proportion*100,gsub("0\\.","",as.character(loh_proportion))))
  print(loh_prop_input)
  non_loh_prop_input=paste0(downsampling_seed,".",ifelse(loh_proportion>=0.1,(1-loh_proportion)*100,gsub("0\\.","",as.character(1-loh_proportion))))
  print(non_loh_prop_input)
  # merge bams with loh_prop and non_loh_prop
  system(paste(samtools,"view -bs",loh_prop_input,paste0(clone0,".bam >") ,paste0(clone0,"_down.bam")))
  system(paste(samtools,"view -bs",non_loh_prop_input,paste0(clone1,".bam >") ,paste0(clone1,"_down.bam")))
  system(paste(samtools,"merge -h",paste0(clone0,"_down.bam"),"-c -p -f",paste0(clone01,".bam"),paste0(clone0,"_down.bam"),paste0(clone1,"_down.bam")))
  system(paste(samtools,"sort",paste0(clone01,".bam -o"),paste0(clone01,"_sorted.bam")))
  system(paste("mv",paste0(clone01,"_sorted.bam"),paste0(clone01,".bam")))
  system(paste(samtools,"index",paste0(clone01,".bam")))
}

#' Function for mixing BAMs for subclonal GAIN events
#' 
#' @param simulated_purity Full path to the chromosome fasta directory
#' @param cna_ccf The tumour purity to be simulated (numeric, ranging (0,1])
#' @param downsampling_seed The seed value to be used for downsampling with samtools
#' @param samtools Full path to the samtools bin
#' @param clone0 Name of the clone with the simulated GAIN
#' @param clone1 Name of the clone with normal diploid copy number
#' @param clone01 Name of the mixed clone to be used for the output BAM
#' @author naser.ansari-pour
#' @export

mix_bams_gain <- function(simulated_purity,cna_ccf,downsampling_seed,samtools,clone0,clone1,clone01){
  print(paste("GAIN CCF =",cna_ccf))
  print(paste("NORMAL CCF =",1-cna_ccf))
  gain_proportion=round(simulated_purity*cna_ccf,digits = 2)
  gain_prop_input=paste0(downsampling_seed,".",ifelse(gain_proportion>=0.1,gain_proportion*100,gsub("0\\.","",as.character(gain_proportion))))
  print(gain_prop_input)
  non_gain_prop_input=paste0(downsampling_seed,".",ifelse(gain_proportion>=0.1,(1-gain_proportion)*100,gsub("0\\.","",as.character(1-gain_proportion))))
  print(non_gain_prop_input)
  # merge bams with gain_prop and non_gain_prop
  system(paste(samtools,"view -bs",gain_prop_input,paste0(clone0,".bam >") ,paste0(clone0,"_down.bam")))
  system(paste(samtools,"view -bs",non_gain_prop_input,paste0(clone1,".bam >") ,paste0(clone1,"_down.bam")))
  system(paste(samtools,"merge -h",paste0(clone0,"_down.bam"),"-c -p -f",paste0(clone01,".bam"),paste0(clone0,"_down.bam"),paste0(clone1,"_down.bam")))
  system(paste(samtools,"sort",paste0(clone01,".bam -o"),paste0(clone01,"_sorted.bam")))
  system(paste("mv",paste0(clone01,"_sorted.bam"),paste0(clone01,".bam")))
  system(paste(samtools,"index",paste0(clone01,".bam")))
}

#' Function for mixing BAMs for chromosome X CNA events
#' 
#' @param simulated_purity Full path to chromosome X fasta directory
#' @param cna_ccf The tumour purity to be simulated (numeric, ranging (0,1])
#' @param downsampling_seed The seed value to be used for downsampling with samtools
#' @param samtools Full path to the samtools bin
#' @param clone0 Name of the clone with the simulated CNA
#' @param clone1 Name of the clone with normal copy number
#' @param clone01 Name of the mixed clone to be used for the output BAM
#' @author naser.ansari-pour
#' @export

mix_bams_chrX <- function(simulated_purity,cna_ccf,downsampling_seed,samtools,clone0,clone1,clone01){
  print(paste("CNA CCF =",cna_ccf))
  print(paste("NORMAL CCF =",1-cna_ccf))
  cna_proportion=round(simulated_purity*cna_ccf,digits = 2)
  cna_prop_input=paste0(downsampling_seed,".",ifelse(cna_proportion>=0.1,cna_proportion*100,gsub("0\\.","",as.character(cna_proportion))))
  print(cna_prop_input)
  non_cna_prop_input=paste0(downsampling_seed,".",ifelse(cna_proportion>=0.1,(1-cna_proportion)*100,gsub("0\\.","",as.character(1-cna_proportion))))
  print(non_cna_prop_input)
  # merge bams with cna_prop and non_cna_prop
  system(paste(samtools,"view -bs",cna_prop_input,paste0(clone0,".bam >") ,paste0(clone0,"_down.bam")))
  system(paste(samtools,"view -bs",non_cna_prop_input,paste0(clone1,".bam >") ,paste0(clone1,"_down.bam")))
  system(paste(samtools,"merge -h",paste0(clone0,"_down.bam"),"-c -p -f",paste0(clone01,".bam"),paste0(clone0,"_down.bam"),paste0(clone1,"_down.bam")))
  system(paste(samtools,"sort",paste0(clone01,".bam -o"),paste0(clone01,"_sorted.bam")))
  system(paste("mv",paste0(clone01,"_sorted.bam"),paste0(clone01,".bam")))
  system(paste(samtools,"index",paste0(clone01,".bam")))
}

#' Function for merging BAM segments/regions to form a full chromosome
#' 
#' This functions merges the pieces/segments of the chromosome based on the bam list file
#' that need to be stitched back together some of which will contain CNA and 
#' some are diploid regions (if any) in between
#' @param samtools Full path to the samtools bin
#' @param mini_bams_list_file Full path to a one-column file listing the bam names to be merged (including the .bam file extension)
#' @param ncores Number of cores to be used in running BWA and samtools
#' @param clone_bam_filename The file name (including the .bam file extension) to be used for the final merged bam output
#' @author naser.ansari-pour
#' @export

merge_mini_bams <- function(samtools,mini_bams_list_file,ncores,clone_bam_filename){
  print(clone_bam_filename)
  first_mini_bam_filename=readLines(mini_bams_list_file)[1]
  system(paste(samtools,"merge -@",ncores,"-h",first_mini_bam_filename,"-c -p -f",clone_bam_filename,"-b",mini_bams_list_file))
  system(paste(samtools,"sort -@",ncores,clone_bam_filename,"-o",paste0(gsub(".bam","",clone_bam_filename),"_sorted.bam")))
  system(paste("mv",paste0(gsub(".bam","",clone_bam_filename),"_sorted.bam"),clone_bam_filename))
  system(paste(samtools,"index -@",ncores,clone_bam_filename))
}

#' Deprecated - Function for mixing BAMs for subclonal LOH events
#' 
#' @param simulated_purity Full path to the chromosome fasta directory
#' @param cna_ccf The tumour purity to be simulated (numeric, ranging (0,1])
#' @param chrom_length Chromosome sequence length (integer value)
#' @param cna_length CNA sequence length (integer value)
#' @param downsampling_seed The seed value to be used for downsampling with samtools
#' @param samtools Full path to the samtools bin
#' @param clone0 Name of the clone with the simulated GAIN
#' @param clone1 Name of the clone with normal diploid copy number
#' @param clone01 Name of the mixed clone to be used for the output BAM
#' @author naser.ansari-pour
#' @keywords internal


mix_bams_loh_legacy <- function(simulated_purity,cna_ccf,chrom_length,cna_length,downsampling_seed,samtools,clone0,clone1,clone01){
  # Calculate loh_prop with precision
  loh_prop=loh_mix_prop(cna_ccf = cna_ccf,
                        cna_length = cna_length,
                        chrom_length = chrom_length)
  loh_prop_input=paste0(downsampling_seed-1,".",loh_prop*1000)
  print(loh_prop_input)
  non_loh_prop_input=paste0(downsampling_seed+1,".",(1-loh_prop)*1000)
  print(non_loh_prop_input)
  #loh_prop_input="--subsample 0.500"
  #non_loh_prop_input="--subsample 0.500"
  # merge bams with loh_prop and non_loh_prop
  system(paste(samtools,"view -bs",loh_prop_input,paste0(clone0,".bam >") ,paste0(clone0,"_down.bam")))
  system(paste(samtools,"view -bs",non_loh_prop_input,paste0(clone1,".bam >") ,paste0(clone1,"_down.bam")))
  # system(paste(samtools,"view -b",loh_prop_input,paste0(clone0,".bam >") ,paste0(clone0,"_down.bam")))
  # system(paste(samtools,"view -b",non_loh_prop_input,paste0(clone1,".bam >") ,paste0(clone1,"_down.bam")))
  system(paste(samtools,"merge -h",paste0(clone0,"_down.bam"),"-c -p -f",paste0(clone01,".bam"),paste0(clone0,"_down.bam"),paste0(clone1,"_down.bam")))
  system(paste(samtools,"sort",paste0(clone01,".bam -o"),paste0(clone01,"_sorted.bam")))
  system(paste("mv",paste0(clone01,"_sorted.bam"),paste0(clone01,".bam")))
  system(paste(samtools,"index",paste0(clone01,".bam")))
}

#' Deprecated - Function for mixing BAMs for subclonal GAIN events
#' 
#' @param cna_ccf The tumour purity to be simulated (numeric, ranging (0,1])
#' @param chrom_length Chromosome sequence length (integer value)
#' @param cna_length CNA sequence length (integer value)
#' @param downsampling_seed The seed value to be used for downsampling with samtools
#' @param samtools Full path to the samtools bin
#' @param clone0 Name of the clone with the simulated GAIN
#' @param clone1 Name of the clone with normal diploid copy number
#' @param clone01 Name of the mixed clone to be used for the output BAM
#' @author naser.ansari-pour
#' @keywords internal

mix_bams_gain_legacy <- function(cna_ccf,chrom_length,cna_length,downsampling_seed,samtools,clone0,clone1,clone01){
  # Calculate loh_prop with precision
  gain_prop=gain_mix_prop(cna_ccf = cna_ccf,
                          cna_length = cna_length,
                          chrom_length = chrom_length)
  gain_prop_input=paste0(downsampling_seed-1,".",cna_ccf*1000)
  print(gain_prop_input)
  non_gain_prop_input=paste0(downsampling_seed+1,".",(1-cna_ccf)*1000)
  print(non_gain_prop_input)
  # gain_prop_input=paste("--subsample",cna_ccf)
  print(paste("GAIN CCF =",cna_ccf))
  # non_gain_prop_input=paste("--subsample",1-cna_ccf)
  print(paste("NORMAL CCF =",1-cna_ccf))
  # merge bams with gain_prop and non_gain_prop
  system(paste(samtools,"view -bs",gain_prop_input,paste0(clone0,".bam >") ,paste0(clone0,"_down.bam")))
  system(paste(samtools,"view -bs",non_gain_prop_input,paste0(clone1,".bam >") ,paste0(clone1,"_down.bam")))
  system(paste(samtools,"merge -h",paste0(clone0,"_down.bam"),"-c -p -f",paste0(clone01,".bam"),paste0(clone0,"_down.bam"),paste0(clone1,"_down.bam")))
  system(paste(samtools,"sort",paste0(clone01,".bam -o"),paste0(clone01,"_sorted.bam")))
  system(paste("mv",paste0(clone01,"_sorted.bam"),paste0(clone01,".bam")))
  system(paste(samtools,"index",paste0(clone01,".bam")))
}

#' Deprecated - Function for calculating clone proportion for normal contamination in non-pure simulated tumour
#' 
#' @param simulated_purity Full path to the chromosome fasta directory
#' @param clone_length_ratio Ratio of the genome length of the clone with CNA to normal diploid genome
#' @author naser.ansari-pour
#' @keywords internal

clone_proportion <- function(simulated_purity,clone_length_ratio){
  clone_prop = (simulated_purity*clone_length_ratio)/(simulated_purity*clone_length_ratio + (1-simulated_purity))
  return(round(clone_prop,digits = 3))
}

#' Deprecated - Function for adding normal contamination to a GAIN clone
#' 
#' @param simulated_purity Full path to the chromosome fasta directory
#' @param downsampling_seed The seed value to be used for downsampling with samtools
#' @param samtools Full path to the samtools bin
#' @param clone Name of the clone with the simulated GAIN
#' @param normal Name of the clone with normal diploid copy number
#' @param clone_normal Name of the mixed clone to be used for the output BAM
#' @author naser.ansari-pour
#' @keywords internal

contaminate_gain_clone_bam <- function(simulated_purity,downsampling_seed,samtools,clone,normal,clone_normal){
  # Calculate clone proportion with precision
  #if (cna_ccf==1){
    clone_proportion=simulated_purity
  #} else if (cna_ccf<1){
  #clone_proportion=clone_prop(simulated_purity = simulated_purity,
  #                      clone_length_ratio = clone_length_ratio)
  #} else {print("CCF not defined properly")}
  clone_prop_input=paste0(downsampling_seed-1,".",clone_proportion*1000)
  print(clone_prop_input)
  normal_prop_input=paste0(downsampling_seed+1,".",(1-clone_proportion)*1000)
  print(normal_prop_input)
  # clone_prop_input=paste0("--subsample 0.",clone_proportion*1000)
  # normal_prop_input=paste0("--subsample 0.",(1-clone_proportion)*1000)
  # merge bams with clone_prop and normal_prop
  system(paste(samtools,"view -bs",clone_prop_input,paste0(clone,".bam >") ,paste0(clone,"_down.bam")))
  system(paste(samtools,"view -bs",normal_prop_input,paste0(normal,".bam >") ,paste0(normal,"_down.bam")))
  # system(paste(samtools,"view -b",clone_prop_input,paste0(clone,".bam >") ,paste0(clone,"_down.bam")))
  # system(paste(samtools,"view -b",normal_prop_input,paste0(normal,".bam >") ,paste0(normal,"_down.bam")))
  system(paste(samtools,"merge -h",paste0(clone,"_down.bam"),"-c -p -f",paste0(clone_normal,".bam"),paste0(clone,"_down.bam"),paste0(normal,"_down.bam")))
  system(paste(samtools,"sort",paste0(clone_normal,".bam -o"),paste0(clone_normal,"_sorted.bam")))
  system(paste("mv",paste0(clone_normal,"_sorted.bam"),paste0(clone_normal,".bam")))
  system(paste(samtools,"index",paste0(clone_normal,".bam")))
}

#' Deprecated - Function for adding normal contamination to an LOH clone
#' 
#' @param simulated_purity Full path to the chromosome fasta directory
#' @param cna_ccf The cancer cell fraction of the copy number to be simulated
#' @param clone_length_ratio Ratio of the genome length of the clone with CNA to normal diploid genome
#' @param downsampling_seed The seed value to be used for downsampling with samtools
#' @param samtools Full path to the samtools bin
#' @param clone Name of the clone with the simulated GAIN
#' @param normal Name of the clone with normal diploid copy number
#' @param clone_normal Name of the mixed clone to be used for the output BAM
#' @author naser.ansari-pour
#' @keywords internal

contaminate_loh_clone_bam <- function(simulated_purity,cna_ccf,clone_length_ratio,downsampling_seed,samtools,clone,normal,clone_normal){
  # Calculate clone proportion with precision
  if (cna_ccf==1){
  loh_clone_proportion=simulated_purity
  } else if (cna_ccf<1){
  loh_clone_proportion=clone_proportion(simulated_purity = simulated_purity,
                        clone_length_ratio = clone_length_ratio)
  } else {print("CCF not defined properly")}
  clone_prop_input=paste0(downsampling_seed,".",loh_clone_proportion*1000)
  print(clone_prop_input)
  normal_prop_input=paste0(downsampling_seed,".",(1-loh_clone_proportion)*1000)
  print(normal_prop_input)
  # clone_prop_input=paste0("--subsample 0.",clone_proportion*1000)
  # normal_prop_input=paste0("--subsample 0.",(1-clone_proportion)*1000)
  # merge bams with clone_prop and normal_prop
  system(paste(samtools,"view -bs",clone_prop_input,paste0(clone,".bam >") ,paste0(clone,"_down.bam")))
  system(paste(samtools,"view -bs",normal_prop_input,paste0(normal,".bam >") ,paste0(normal,"_down.bam")))
  # system(paste(samtools,"view -b",clone_prop_input,paste0(clone,".bam >") ,paste0(clone,"_down.bam")))
  # system(paste(samtools,"view -b",normal_prop_input,paste0(normal,".bam >") ,paste0(normal,"_down.bam")))
  system(paste(samtools,"merge -h",paste0(clone,"_down.bam"),"-c -p -f",paste0(clone_normal,".bam"),paste0(clone,"_down.bam"),paste0(normal,"_down.bam")))
  system(paste(samtools,"sort",paste0(clone_normal,".bam -o"),paste0(clone_normal,"_sorted.bam")))
  system(paste("mv",paste0(clone_normal,"_sorted.bam"),paste0(clone_normal,".bam")))
  system(paste(samtools,"index",paste0(clone_normal,".bam")))
}
