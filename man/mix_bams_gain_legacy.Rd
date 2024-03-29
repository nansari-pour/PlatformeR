% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mix_and_merge_bams.R
\name{mix_bams_gain_legacy}
\alias{mix_bams_gain_legacy}
\title{Deprecated - Function for mixing BAMs for subclonal GAIN events}
\usage{
mix_bams_gain_legacy(
  cna_ccf,
  chrom_length,
  cna_length,
  downsampling_seed,
  samtools,
  clone0,
  clone1,
  clone01
)
}
\arguments{
\item{cna_ccf}{The tumour purity to be simulated (numeric, ranging (0,1])}

\item{chrom_length}{Chromosome sequence length (integer value)}

\item{cna_length}{CNA sequence length (integer value)}

\item{downsampling_seed}{The seed value to be used for downsampling with samtools}

\item{samtools}{Full path to the samtools bin}

\item{clone0}{Name of the clone with the simulated GAIN}

\item{clone1}{Name of the clone with normal diploid copy number}

\item{clone01}{Name of the mixed clone to be used for the output BAM}
}
\description{
Deprecated - Function for mixing BAMs for subclonal GAIN events
}
\author{
naser.ansari-pour
}
\keyword{internal}
