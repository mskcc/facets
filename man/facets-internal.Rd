\name{facets-internal}
\alias{counts2logROR}
\alias{jointsegsummary}
\alias{mafmle}
\alias{procSnps}
\alias{procSnpsDT}
\alias{scanSnp}
\alias{fit.cpt.tree}
\alias{prune.cpt.tree}
\alias{segsnps}
\alias{clustersegs}
\alias{hg19gcpct}
\alias{gcpctdb}
\alias{onepass}
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
mafmle(fcount, vaf)
procSnps(filename, ndepth=35, het.thresh=0.25, snp.nbhd=250, chromlevels=c(1:22,"X"), unmatched=FALSE)
procSnpsDT(filename, ndepth=35, het.thresh=0.25, snp.nbhd=250, chromlevels=c(1:22,"X"), unmatched=FALSE)
counts2logROR(mat, unmatched=FALSE, f=0.2)
scanSnp(maploc, het, nbhd)
fit.cpt.tree(genomdat, edgelim=10, cval=25, hscl = 1)
prune.cpt.tree(seg.tree, cval=25)
segsnps(mat, cval=25, hetscale = FALSE)
clustersegs(out, jointseg, min.nhet=10)
find.mode(x)
onepass(x, trace, unif, rho, rhov, prior, posterior, sigma, min.nhet, rho.clust, maxiter, eps)
data(hg19gcpct)
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
