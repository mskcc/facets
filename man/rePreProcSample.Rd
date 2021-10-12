\name{rePreProcSample}
\alias{rePreProcSample}
\title{Re-(pre)-process a sample}
\description{
  Re-process jointseg data frame from procSample output to generate a
  segmentation tree for further analysis.
}
\usage{
  rePreProcSample(jseg, cval=25, deltaCN=0, gbuild=c("hg19", "hg38",
       "hg18", "mm9", "mm10"), hetscale=TRUE, unmatched=FALSE)
}
\arguments{
  \item{jseg}{data frame with 8 required columns: \code{chrom},
    \code{maploc}, \code{rCountT}, \code{vafT}, \code{rCountN},
    \code{vafN}, \code{het} and \code{keep}. Rest are ignored.}
  \item{cval}{critical value for segmentation}
  \item{deltaCN}{minimum detectable difference in CN from diploid state}
  \item{gbuild}{genome build used for the alignment of the genome.
    Default value is human genome build hg19. Other possibilities are
    hg38 & hg18 for human and mm9 & mm10 for mouse. Chromosomes used for
    analysis are \code{1-22, X} for humans and \code{1-19} for mouse.
    Option udef in preProcSample not implemented.}
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
\details{
  The output from \code{procSample} can be processed to generate the R
  object needed for further analysis. This function will create what is
  necessary for \code{procSample} can be applied. Usage is  
  \code{xx <- rePreProcSample(zout$jointseg[,1:8])}
  \code{oo <- procSample(xx, ...)}
}
