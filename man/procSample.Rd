\name{procSample}
\alias{procSample}
\title{Process a sample}
\description{
  Processes the output from preProcSample for given cval and min.nhet
}
\usage{
  procSample(x, cval=35, min.nhet=25)
}
\arguments{
  \item{x}{the output from preProcSample; should contain seg.tree and jointseg}
  \item{cval}{critical value for segmentation}
  \item{min.nhet}{minimum number of heterozygote snps in a segment used
    for bivariate t-statistic during clustering of segments}
}
\value{
  A list consisting of two elements:
  \item{jointseg}{The data that were segmented. Only the loci that were
    sampled within a snp.nbhd are present}
  \item{out}{data frame of segment summaries pre and post clustering of
    segments}
}
