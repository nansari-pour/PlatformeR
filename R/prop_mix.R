#' Non-piecewise function for calculating the proportion of downsampling from LOH and normal BAMs
#' 
#' Correct function to calculate the prop for the LOH bam for subclonal events
#' @param cna_ccf The tumour purity to be simulated (numeric, ranging (0,1])
#' @param chrom_length Chromosome sequence length (integer value)
#' @param cna_length CNA sequence length (integer value)
#' @author naser.ansari-pour
#' @keywords internal

loh_mix_prop <- function(cna_ccf,cna_length,chrom_length){
  loh_frac=(cna_ccf)*((2*chrom_length-cna_length)/(2*chrom_length))
  non_loh_frac=(1-cna_ccf)*1
  return(round(loh_frac/(loh_frac+non_loh_frac),digits = 3))
}

#' Non-piecewise function for calculating the proportion of downsampling from GAIN and normal BAMs
#' 
#' Correct function to calculate the prop for the GAIN bam for subclonal events
#' @param cna_ccf The tumour purity to be simulated (numeric, ranging (0,1])
#' @param chrom_length Chromosome sequence length (integer value)
#' @param cna_length CNA sequence length (integer value)
#' @author naser.ansari-pour
#' @keywords internal

gain_mix_prop <- function(cna_ccf,cna_length,chrom_length){
  gain_frac=(cna_ccf)*((2*chrom_length+cna_length)/(2*chrom_length))
  non_gain_frac=(1-cna_ccf)*1
  return(round(gain_frac/(gain_frac+non_gain_frac),digits = 3))
}
