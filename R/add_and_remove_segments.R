#' Function for removing a segment from the chromosome fasta sequence 
#'
#' Non-piecewise function to remove segments due to LOH from a chromosome based on start and end positions (using alternate_segments)
#' @param chromosome_string Chromosome fasta sequence string (one contiguous string)
#' @param loh_segments A two column dataframe with headers startpos and endpos providing regions to be simulated for LOH
#' @author naser.ansari-pour
#' @export


remove_segments <- function(chromosome_string, loh_segments) {
  
  chromosome_length=nchar(chromosome_string)
  # Ensure the input string length is greater than 0
  if (chromosome_length == 0) {
    stop(print("Chromosome string is empty"))  # stop if chromosome string is empty
  }
  
  # Sort segments by start position in descending order
  loh_segments <- loh_segments[order(loh_segments$startpos),]
  
  # Identify segments to keep and generate retained string
  segments_to_keep=alternate_segments(chromosome_length = chromosome_length,
                                      segments_to_remove = loh_segments)
  
  if (!is.null(segments_to_keep)){
    post_loh_string=NULL
    for (i in 1:nrow(segments_to_keep)){
      post_loh_string=paste0(post_loh_string,substr(chromosome_string,segments_to_keep$startpos[i],segments_to_keep$endpos[i]),collapse = "")
    }
  } else {post_loh_string=NULL} # post_loh_string=NULL i.e. entire chromosome is lost
  return(post_loh_string)
}

#' Function for adding a segment from the chromosome fasta sequence 
#'
#' Non-piecewise function to add segment due to GAIN at the end of a chromosome fasta sequence based on start and end positions (using alternate_segments)
#' @param chromosome_string Chromosome fasta sequence string (one contiguous string)
#' @param gain_segments A two column dataframe with headers startpos and endpos providing regions to be simulated for GAIN
#' @author naser.ansari-pour
#' @export

add_segments <- function(chromosome_string, gain_segments) {
  
  chromosome_length=nchar(chromosome_string)
  # Ensure the input string length is greater than 0
  if (chromosome_length == 0) {
    stop(print("Chromosome string is empty"))  # stop if chromosome string is empty
  }
  
  # Sort segments by start position in descending order
  gain_segments <- gain_segments[order(gain_segments$startpos),]
  
  # add gain segments to chromosome string
  post_gain_string=chromosome_string
  for (i in 1:nrow(gain_segments)){
    post_gain_string=paste0(post_gain_string,paste0(rep(substr(chromosome_string,gain_segments$startpos[i],gain_segments$endpos[i]),gain_segments$TCN[i]-1),collapse = ""),collapse = "")
  }
  return(post_gain_string)
}
