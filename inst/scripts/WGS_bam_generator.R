# PlatformeR genome-wide/WGS merge bam run

# add full Library path of PlatformeR installation to .libPaths()
.libPaths("/full/path/to/R/Library/where/PlatformeR/is/installed")

library(PlatformeR)

Args=commandArgs(trailingOnly = TRUE)
chrom_names = unlist(strsplit(toString(Args[1]), split = " "))
print(chrom_names)
simulated_purity = as.numeric(toString(Args[2]))
print(paste("Simulated purity =",simulated_purity))
out_dir = toString(Args[3])
print(out_dir)
tumourname = toString(Args[4])
samtools = toString(Args[5])
ncores = as.integer(toString(Args[6]))
# Optional - paired analysis
merge_normal = toString(Args[7])
# if merge_normal = TRUE, then:
if (merge_normal){
normalname = toString(Args[8])
print(normalname)
}

setwd(out_dir)
print(getwd())

genomewide_merge_bam_list(chrom_names = chrom_names,
                          simulated_purity = simulated_purity,
                          out_dir = out_dir,
                          samplename = tumourname,
                          is.normal = FALSE)

merge_mini_bams(samtools = samtools,
                mini_bams_list_file = paste0(out_dir,tumourname,"_genomewide_bam_list.txt"),
                ncores = ncores,
                clone_bam_filename = ifelse(simulated_purity==1,paste0(tumourname,".bam"),paste0(tumourname,"_purity_",simulated_purity,".bam")))

if (merge_normal){
  
  genomewide_merge_bam_list(chrom_names = chrom_names,
                            simulated_purity = NULL,
                            out_dir = out_dir,
                            samplename = normalname,
                            is.normal = TRUE)

  merge_mini_bams(samtools = samtools,
                  mini_bams_list_file = paste0(out_dir,normalname,"_genomewide_bam_list.txt"),
                  ncores = ncores,
                  clone_bam_filename = paste0(normalname,".bam"))
}
