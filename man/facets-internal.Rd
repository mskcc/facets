\name{facets-internal}
\alias{counts2logROR}
\alias{jointsegsummary}
\alias{mafmle}
\alias{procSnps}
\alias{scanSnp}
\alias{fit.cpt.tree}
\alias{prune.cpt.tree}
\alias{segsnps}
\alias{clustersegs}
\alias{clusteredcncf.fit}
\alias{mwstat}
\alias{mergexy}
\alias{hg19gcpct}
\alias{gcpctdb}
\alias{find.mode}
\title{facets functions}
\description{
  Internal functions of package facets.
}
\usage{
jointsegsummary(jointseg)
mafmle(fcount, vaf)
procSnps(filename, ndepth=35, het.thresh=0.25, snp.nbhd=250, chromlevels=c(1:22,"X"))
counts2logROR(mat, f=0.2)
scanSnp(maploc, het, nbhd)
fit.cpt.tree(genomdat, edgelim=10, cval=25)
prune.cpt.tree(seg.tree, cval=25)
segsnps(mat, cval=25)
clustersegs(jointseg, out, min.nhet=25, cval=35)
clusteredcncf.fit(out, dipLogR=0)
mwstat(x,y)
mergexy(x,y)
find.mode(x)
data(hg19gcpct)
}
\details{These are not to be called directly by the user}
\keyword{internal}
