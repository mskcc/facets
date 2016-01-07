\name{emcncf}
\alias{emcncf}
\title{EM estimate of copy number and cellular fraction of segment clusters}
\description{
  Uses genotype mixture model to estimate the cluster specific copy number and cellular fraction. 
  Uses estimates based on the cnlr.median and mafR as initial values for the EM iteration.
} 
\usage{
  emcncf(x, trace=FALSE, unif=FALSE, min.nhet=15, maxiter=10, eps=1e-3)
}
\arguments{
  \item{x}{the output from procSample. This function uses the elements
    jointseg, out and dipLogR from the output.}
  \item{trace}{flag to print the EM criteria at every step}
  \item{unif}{random EM start values of cellular fractions instead of
    clusteredcncf values}
  \item{min.nhet}{minimum number of heterozygote snps in a segment used to call minor cn}
  \item{maxiter}{maximum number of EM iterations}
  \item{eps}{the convergence threshold}
}
\value{
  A list containing:
  \item{loglik}{loglikelihood value of the fitted model}
  \item{purity}{fraction tumor cells in the tumor sample}
  \item{ploidy}{average total copy number of the tumor cells}
  \item{dipLogR}{estimated logR value of diploid segments}
  \item{cncf}{dataframe consisting of the columns of segmentation output as well as
    cellular fraction (cf), total (tcn) and lesser (lcn) copy number of
    each segment and their em counterpart (with .em suffix)}
}
