\name{procSample}
\alias{procSample}
\title{Process a sample}
\description{
  Processes the output from preProcSample for given cval and min.nhet
}
\usage{
  procSample(x, cval=150, min.nhet=15, dipLogR=NULL)
}
\arguments{
  \item{x}{the output from preProcSample; should contain seg.tree and jointseg}
  \item{cval}{critical value for segmentation}
  \item{min.nhet}{minimum number of heterozygote snps in a segment used
    for bivariate t-statistic during clustering of segments}
  \item{dipLogR}{diploid level obtained from a fit, typically using a
    higher cval, can be used with lower cval to recover focal changes}
}
\value{
  A list consisting of two elements:
  \item{jointseg}{The data that were segmented. Only the loci that were
    sampled within a snp.nbhd are present}
  \item{out}{data frame of segment summaries pre and post clustering of
    segments. The columns are: \code{chrom} the chromosome to which the
    segment belongs; \code{seg} the segment number; \code{num.mark} the
    number of SNPs in the segment; \code{nhet} the number of SNPs that
    are deemed heterozygous; \code{cnlr.median} the median log-ratio of
    the segment; \code{mafR} the log-odds-ratio summary for the segment;
    \code{segclust} the segment cluster to which segment belongs;
    \code{cnlr.median.clust} the median log-ratio of the segment
    cluster; \code{mafR.clust} the log-odds-ratio summary for the
    segment cluster; \code{cf} the cellular fraction of the segment;
    \code{tcn} the total copy number of the segment; \code{lcn} the
    minor copy number of the segment.
  }
  \item{dipLogR}{specified or obtained from data}
  \item{...}{other output when findDiploidLogR is used}
}
\details{
  The minor copy number \code{lcn} may not be estimated with confidence
  when a segment has fewer than \code{min.nhet} heterozygous SNPs and 
  hence will return \code{NA}. If there are too few heterozygous SNPs in
  a segment then \code{mafR} and \code{mafR.clust} can be \code{NA}.
}
