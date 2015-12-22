\name{plotSampleCNCF}
\alias{plotSampleCNCF}
\title{Plot the data and results for a single sample}
\description{
  Plots copy number log-ratio, variant allele log-odds ratio as well as
  fitted cellular fraction and integer copy numbers.
} 
\usage{
  plotSampleCNCF(x, fit, plotboth=FALSE)
}
\arguments{
  \item{x}{output of procSample with elements jointseg and out}
  \item{fit}{output of emcncf}
  \item{plotboth}{plot both the initial and the refined estimates}
}
