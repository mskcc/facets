\name{logRlogORspider}
\alias{logRlogORspider}
\title{Plot segment logR-logOR summary statistics}
\description{
  Plots copy number log-ratio, variant allele log-odds ratio as well as 
  the copy number and cellular fraction fits.
} 
\usage{
  logRlogORspider(cncf, dipLogR=0, nfrac=0.005)
}
\arguments{
  \item{cncf}{Copy number and cellular fraction data frame either the
    naive one (out) from procSample or the EM fit (cncf) from emcncf.}
  \item{dipLogR}{the log-ratio value corresponding to the diploid state.}
  \item{nfrac}{a segment is shown if the proportion of loci and het SNPs
    (\code{num.mark} and \code{nhet}) nfrac. Default is 0.01.}
}
\details{
  This is a diagnostic plot to check how well the copy number fits work.
  The estimated segment summaries are plotted as circles where the size
  of the circle increases with the number of loci in the segment. The
  expected value for various integer copy number states are drawn as
  curves for purity ranging from 0 to 0.95. For a good fit, the segment
  summaries should be close to one of the lines.
}
