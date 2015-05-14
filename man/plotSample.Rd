\name{plotSample}
\alias{plotSample}
\title{Plot the data and results for a single sample}
\description{
  Plots copy number log-ratio, variant allele log-odds ratio
} 
\usage{
  plotSample(x, clustered=FALSE, chromlevels=c(1:22,"X"))
}
\arguments{
  \item{x}{output from procSample}
  \item{clustered}{indicator of whether segment or cluster summary plotted}
  \item{chromlevels}{names of chromosomes. Chromosomes will be ordered
    in the order given. Useful for non-human genome.}
}
