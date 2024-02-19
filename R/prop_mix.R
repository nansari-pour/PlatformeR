# Non-piecewise functions for calculating the proportion of downsampling from LOH/GAIN and normal BAMs

# Correct function to calculate the prop for the LOH bam for subclonal events

loh_mix_prop <- function(cna_ccf,cna_length,chrom_length){
  loh_frac=(cna_ccf)*((2*chrom_length-cna_length)/(2*chrom_length))
  non_loh_frac=(1-cna_ccf)*1
  return(round(loh_frac/(loh_frac+non_loh_frac),digits = 3))
}

# Correct function to calculate the prop for the GAIN bam for subclonal events

gain_mix_prop <- function(cna_ccf,cna_length,chrom_length){
  gain_frac=(cna_ccf)*((2*chrom_length+cna_length)/(2*chrom_length))
  non_gain_frac=(1-cna_ccf)*1
  return(round(gain_frac/(gain_frac+non_gain_frac),digits = 3))
}
