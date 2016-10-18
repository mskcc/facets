\name{plotSample}
\alias{plotSample}
\title{Plot the data and results for a single sample}
\description{
  Plots copy number log-ratio, variant allele log-odds ratio as well as 
  the copy number and cellular fraction fits.
} 
\usage{
  plotSample(x, emfit=NULL, clustered=FALSE, plot.type=
                c("em","naive","both","none"), sname=NULL)
}
\arguments{
  \item{x}{output from procSample}
  \item{emfit}{output of emcncf}
  \item{clustered}{indicator of whether segment or cluster summary plotted}
  \item{plot.type}{the type of plot. The default is \code{em} in which
    the logR and logOR data as well as the copy number and cellular 
    fraction fits from EM are graphed. For \code{naive} the naive copy
    number and cellular fraction fits are used instead of EM. For
    \code{none} only the data are shown and for \code{both} both fits
    are shown.}
  \item{sname}{sample name give as a character string.}
}
