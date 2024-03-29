% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_and_remove_segments.R
\name{remove_segments}
\alias{remove_segments}
\title{Function for removing a segment from the chromosome fasta sequence}
\usage{
remove_segments(chromosome_string, loh_segments)
}
\arguments{
\item{chromosome_string}{Chromosome fasta sequence string (one contiguous string)}

\item{loh_segments}{A two column dataframe with headers startpos and endpos providing regions to be simulated for LOH}
}
\description{
Non-piecewise (legacy) function to remove segments due to LOH from a chromosome based on start and end positions (using alternate_segments)
}
\author{
naser.ansari-pour
}
\keyword{internal}
