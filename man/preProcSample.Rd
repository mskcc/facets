\name{preProcSample}
\alias{preProcSample}
\title{Pre-process a sample}
\description{
  Processes a snp read count matrix and generates a segmentation tree
}
\usage{
  preProcSample(rcmat, ndepth=35, het.thresh=0.25, snp.nbhd=250, cval=25,
       deltaCN=0, gbuild=c("hg19", "hg38", "hg18", "mm9", "mm10", "udef"),
       ugcpct = NULL, hetscale=TRUE, unmatched=FALSE, ndepthmax=1000)
}
\arguments{
  \item{rcmat}{data frame with 6 required columns: \code{Chrom},
    \code{Pos}, \code{NOR.DP}, \code{NOR.RD}, \code{TUM.DP} and
    \code{TUM.RD}. Additional variables are ignored.}
  \item{ndepth}{minimum normal sample depth to keep}
  \item{het.thresh}{vaf threshold to call a SNP heterozygous}
  \item{snp.nbhd}{window size}
  \item{cval}{critical value for segmentation}
  \item{deltaCN}{minimum detectable difference in CN from diploid state}
  \item{gbuild}{genome build used for the alignment of the genome.
    Default value is human genome build hg19. Other possibilities are
    hg38 & hg18 for human and mm9 & mm10 for mouse. Chromosomes used for
    analysis are \code{1-22, X} for humans and \code{1-19} for mouse.
    Option udef can be used to analyze other genomes.}
  \item{ugcpct}{If udef is chosen for gbuild then appropriate GC
    percentage date should be provided through this option. This is a
    list of numeric vectors that gives the GC percentage windows of
    width 1000 bases in steps of 100 i.e. 1-1000, 101-1100 etc. for the
    autosomes and the X chromosome.}
  \item{hetscale}{logical variable to indicate if logOR should get more
    weight in the test statistics for segmentation and clustering. Usually
    only 10\% of snps are hets and hetscale gives the logOR contribution
    to T-square as 0.25/proportion of hets.}
  \item{unmatched}{indicator of whether the normal sample is unmatched.
    When this is TRUE hets are called using tumor reads only and logOR
    calculations are different. Use het.thresh = 0.1 or lower when TRUE.}
  \item{ndepthmax}{loci for which normal coverage exceeds this number
    (default is 1000) will be discarded as PCR duplicates. Fof high
    coverage sample increase this and ndepth commensurately.}
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
  The SNPs in a genome are not evenly spaced. Some regions have multiple
  SNPs in a small neighborhood. Thus using all loci will induce serial
  correlation in the data. To avoid it we sample loci such that only a
  single locus is used in an interval of length \code{snp.nbhd}. So in
  order to get reproducible results use \code{set.seed} to fix the
  random number generator seed.
}
