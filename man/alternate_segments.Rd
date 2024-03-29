% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alternate_segments.R
\name{alternate_segments}
\alias{alternate_segments}
\title{Function to identify alternate segments after CNA events - updated version to allow completely adjacent CNA}
\usage{
alternate_segments(chromosome_length, segments_to_remove)
}
\arguments{
\item{chromosome_length}{Chromosome sequence length (integer value)}

\item{segments_to_remove}{A two column dataframe (chrom_simulate) with headers startpos and endpos providing regions to be simulated with CNA}
}
\description{
Function to identify chromosome segments that are not to be simulated with CNA
}
\author{
naser.ansari-pour
}
