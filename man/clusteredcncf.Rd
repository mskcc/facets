\name{clusteredcncf}
\alias{clusteredcncf}
\title{Copy number and cellular fraction of segment clusters}
\description{
  Fits the cluster specific copy number and cellular fraction based on
  the cnlr.median and mafR values
} 
\usage{
  clusteredcncf(out, dipLogR=0)
}
\arguments{
  \item{out}{the out element from procSample fit}
  \item{dipLogR}{log-ratio level corresponding to diploid copy number}
}
\value{
  A data frame with the same columns as out and three additional columns
  which give the cellular fraction (cf), total copy number (tcn) and
  lower copy number (lcn) for each segment. If a cluster does not have
  sufficient heterozygotes the lcn value will be NA.
}
