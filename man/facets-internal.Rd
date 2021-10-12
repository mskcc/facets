\name{facets-internal}
\alias{counts2logROR}
\alias{jointsegsummary}
\alias{procSnps}
\alias{scanSnp}
\alias{fit.cpt.tree}
\alias{prune.cpt.tree}
\alias{segsnps}
\alias{clustersegs}
\alias{gcpctdb}
\alias{find.mode}
\alias{findDiploidLogR}
\alias{km2class}
\alias{acnsplit}
\alias{dlrdev}
\alias{fitcncf0}
\alias{optcfutil}
\alias{mergecf}
\title{facets functions}
\description{
  Internal functions of package facets.
}
\usage{
jointsegsummary(jointseg)
procSnps(rcmat, ndepth=35, het.thresh=0.25, snp.nbhd=250, nX=23,
unmatched=FALSE, ndepthmax=1000)
counts2logROR(out, gbuild, unmatched=FALSE, ugcpct = NULL, f=0.2)
scanSnp(maploc, het, nbhd)
fit.cpt.tree(genomdat, edgelim=10, cval=25, hscl=1, delta=0)
prune.cpt.tree(seg.tree, cval=25)
segsnps(mat, cval=25, hetscale = FALSE, delta=0)
clustersegs(out, jointseg, min.nhet=10)
find.mode(x)
findDiploidLogR(out, cnlr)
km2class(cnlr.med, num.mark)
acnsplit(lr0, lr1, lorsq)
dlrdev(lr2, lr4, out)
fitcncf0(out, dipLogR=0)
optcfutil(tcn, ocn, maf)
mergecf(out)
}
\details{These are not to be called directly by the user}
\keyword{internal}
