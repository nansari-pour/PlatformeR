# Function to identify alternate segments after CNA events - updated version to allow completely adjacent CNA

alternate_segments <- function(chromosome_length,segments_to_remove){
  segments_to_remove=segments_to_remove[order(segments_to_remove$startpos),]
  non_CNA=data.frame(startpos=numeric(),endpos=numeric()) # get all non_CNA regions  
  for (seg in 1:(nrow(segments_to_remove)+1)){
    # start <- segments_to_remove$startpos[seg]
    # end <- segments_to_remove$endpos[seg]
    # if (start <= end && start >= 1 && end <= chromosome_length) {
    if (seg == 1 & segments_to_remove$startpos[seg]>1){
      non_cna=data.frame(startpos=1,endpos=segments_to_remove$startpos[seg]-1)
    } else if (seg == 1 & segments_to_remove$startpos[seg]==1){
      next
    } else if (seg>1 & seg <= nrow(segments_to_remove)){
      if (segments_to_remove$startpos[seg]==segments_to_remove$endpos[seg-1]+1){
        next
      } else {
      non_cna=data.frame(startpos=segments_to_remove$endpos[seg-1]+1,endpos=segments_to_remove$startpos[seg]-1)
      }
    } else if (seg == nrow(segments_to_remove)+1){
      if (segments_to_remove$endpos[seg-1]==chromosome_length){
        print("reached end of chromosome")
        next
      } else {
        non_cna=data.frame(startpos=segments_to_remove$endpos[seg-1]+1,endpos=chromosome_length)
      }
    }
    non_CNA=rbind(non_CNA,non_cna)
    rm(non_cna)
  }
  if (nrow(non_CNA)>0){
    return(non_CNA)
  } else {non_CNA=NULL}
}
