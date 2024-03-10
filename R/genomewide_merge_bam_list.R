#' Function to generate a list of bams to be merged for a genomewide level bam
#' 
#' @param chrom_names Name of the chromosomes to be used for merging their respective chromosome-level bams 
#' @param simulated_purity The tumour purity to be simulated (numeric, ranging (0,1])
#' @param out_dir Full path to the out directory where genomewide files will be written
#' @param samplename Sample name to be used for the final genomewide bam file with simulated CNA
#' @param is.normal set to TRUE if the bam list is to be generated for normal germline chromosome-level bams
#' @author naser.ansari-pour
#' @export

genomewide_merge_bam_list <- function(chrom_names,simulated_purity,out_dir,samplename,is.normal){
  chrom_names=gsub("chr","",chrom_names)
  if (is.normal){
    for (chrom in chrom_names){
      bam_name=paste0(out_dir,"chr",chrom,"/chr",chrom,".normal.bam")
      # write bam_name to a single file - input to merge_mini_bams function in out_dir
      cat(bam_name, sep = "\n", file = paste0(out_dir,samplename,"_genomewide_bam_list.txt"), append = TRUE)
    }
  } else {
  for (chrom in chrom_names){
    if (simulated_purity<1){
      bam_name=paste0(out_dir,"chr",chrom,"/chr",chrom,".clone_purity_",simulated_purity,".bam")
    } else {
      bam_name=paste0(out_dir,"chr",chrom,"/chr",chrom,".clone.bam")
    }
    # write bam_name to a single file - input to merge_mini_bams function in out_dir
    cat(bam_name, sep = "\n", file = paste0(out_dir,samplename,"_genomewide_bam_list.txt"), append = TRUE)
  }
  }
}
