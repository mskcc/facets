\name{preProcSample}
\alias{preProcSample}
\title{Pre-process a sample}
\description{
  Processes a snp read count matrix and generates a segmentation tree
}
\usage{
  preProcSample(filename, ndepth=35, het.thresh=0.25, snp.nbhd=250,
                cval=25, chromlevels=c(1:22,"X"), hetscale=TRUE,
                unmatched=FALSE)
}
\arguments{
  \item{filename}{input matrix of read depths and variant allele frequencies}
  \item{ndepth}{minimum normal sample depth to keep}
  \item{het.thresh}{vaf threshold to call a SNP heterozygous}
  \item{snp.nbhd}{window size}
  \item{cval}{critical value for segmentation}
  \item{chromlevels}{names of chromosomes. Chromosomes will be ordered
    in the order given. Useful for non-human genome.}
  \item{hetscale}{logical variable to indicate if logOR should get more
    weight in the test statistics for segmentation and clustering. Usually
    only 10\% of snps are hets and hetscale gives the logOR contribution
    to T-square as 0.25/proportion of hets.}
  \item{unmatched}{indicator of whether the normal sample is unmatched.
    When this is TRUE hets are called using tumor reads only and logOR
    calculations are different. Use het.thresh = 0.1 or lower when TRUE.}
}
\value{
  A list consisting of three elements:
  \item{pmat}{Read counts and other elements of all the loci}
  \item{seg.tree}{a list of matrices one for each chromosome. the matrix
    gives the tree structure of the splits. each row corresponds to a
    segment with the parent row as the first element the start-1 and end
    index of each segment and the maximal T^2 statistic. the first row
    is the whole chromosome and its parent row is by definition 0.}
  \item{jointseg}{The data that were segmented. Only the loci that were
    sampled within a snp.nbhd are present. segment results given.}
  \item{hscl}{scaling factor for logOR data.}
}
