\name{emcncf}
\alias{emcncf}
\title{EM estimate of copy number and cellular fraction of segment clusters}
\description{
  Fits the cluster specific copy number and cellular fraction based on
  the cnlr.median and mafR values and uses EM to align with other clusters
} 
\usage{
  emcncf(x, trace=FALSE, unif=FALSE, maxiter=10, eps=1e-3)
}
\arguments{
  \item{x}{the output from procSample. This function uses the elements
    jointseg, out and dipLogR from the output.}
  \item{trace}{flag to print the EM criteria at every step}
  \item{unif}{random EM start values of cellular fractions instead of
    clusteredcncf values}
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
