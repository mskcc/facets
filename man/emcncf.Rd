\name{emcncf}
\alias{emcncf}
\title{EM estimate of copy number and cellular fraction of segment clusters}
\description{
  Fits the cluster specific copy number and cellular fraction based on
  the cnlr.median and mafR values and uses EM to align with other clusters
} 
\usage{
  emcncf(jointseg, out, trace=F, unif=F, dipLogR=NA, dipt=NA, maxiter=10,
            eps=1e-4)
}
\arguments{
  \item{jointseg}{the jointseg element from procSample fit}
  \item{out}{the out element from procSample fit}
  \item{trace}{flag to print the EM criteria at every step}
  \item{unif}{random EM start values of cellular fractions instead of
    clusteredcncf values}
  \item{dipLogR}{log-ratio level corresponding to diploid copy number}
  \item{dipt}{whole genome doubling value. 2 is not doubled. 4 is
    doubled. default is NA and value determined from data.}
  \item{maxiter}{ maximum number of EM iterations}
  \item{eps}{the convergence threshold}
}
\value{
  A list containing:
  \item{loglik}{loglikelihood value of the fitted model}
  \item{purity}{fraction tumor cells in the tumor sample}
  \item{ploidy}{average total copy number of the tumor cells}
  \item{dipLogR}{estimated logR value of diploid segments}
  \item{dipt}{indicator of genome doubling: 2=normal, 4=doubled}
  \item{cncf}{dataframe consisting of the columns of out as well as
    cellular fraction (cf), total (tcn) and lesser (lcn) copy number of
    each segment and their em counterpart (with .em suffix)}
}
