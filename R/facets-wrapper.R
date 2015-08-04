preProcSample <- function(filename, ndepth=35, het.thresh=0.25, snp.nbhd=250, cval=25, chromlevels=c(1:22,"X"), hetscale=TRUE) {
    pmat <- procSnpsDT(filename, ndepth, het.thresh, snp.nbhd, chromlevels)
    dmat <- counts2logROR(pmat[pmat$rCountT>0,])
    tmp <- segsnps(dmat, cval, hetscale)
    out <- list(pmat=pmat)
    c(out, tmp)
}

procSample <- function(x, cval=50, min.nhet=15, dipLogR=NULL) {
    # ensure availability of seg.tree
    if (is.null(x$seg.tree)) stop("seg.tree is not available")
    # make sure that original cval is smaller than current one
    cval.fit <- attr(x$seg.tree, "cval")
    if (cval.fit > cval) stop("original fit used cval = ", cval.fit)
    # jointseg etc
    jseg <- x$jointseg
    jseg <- jseg[is.finite(jseg$cnlr),]
    # loop through chromosomes to get segments
    nchr <- length(x$seg.tree)
    # number the segs from 1 to the total number
    nsegs <- 0
    for (i in 1:nchr) {
        seg.widths <- diff(prune.cpt.tree(x$seg.tree[[i]], cval))
        jseg$seg[jseg$chrom==i] <- nsegs + rep(1:length(seg.widths), seg.widths)
        nsegs <- nsegs + length(seg.widths)
    }
    out <- jointsegsummary(jseg)
    out <- clustersegs(out, jseg, min.nhet)
    # put in the clustered values for snps
    jseg$segclust[is.finite(jseg$cnlr)] <- rep(out$segclust, out$num.mark)
    # find dipLogR and fit cncf
    if (is.null(dipLogR)) {
        oo <- findDiploidLogR(out, jseg$cnlr)
    } else {
        oo <- list()
        oo$out0 <- "empty"
        oo$dipLogR <- dipLogR
    }
    out <- fitcncf(out, oo$dipLogR)
    c(list(jointseg=jseg, out=out), oo[-1])
}

plotSample <- function(x, clustered=FALSE, chromlevels=c(1:22,"X")) {
    def.par <- par(no.readonly = TRUE) # save default, for resetting...
    # raw data used for joint segmentation
    jseg <- x$jointseg
    out <- x$out
    # determine which of the cnlr.median & mafR to show
    if (clustered) {
        cnlr.median <- out$cnlr.median.clust
        mafR <- out$mafR.clust
        mafR[is.na(mafR)] <- out$mafR[is.na(mafR)]
    } else {
        cnlr.median <- out$cnlr.median
        mafR <- out$mafR
    }
    mafR <- abs(mafR)
    # chromosome colors
    chrcol <- 1+rep(out$chrom-2*floor(out$chrom/2), out$num.mark)
    nn <- cumsum(table(jseg$chrom[is.finite(jseg$cnlr)]))
    segbdry <- cumsum(c(0,out$num.mark))
    segstart <- segbdry[-length(segbdry)]
    segend <- segbdry[-1]
    # layout of 2 panel figure
    layout(matrix(c(1,1,2,2), ncol=1))
    par(mar=c(0.25,3,0.25,1), mgp=c(2, 0.7, 0), oma=c(3,0,1.25,0))
    # plot the logR data and segment medians
    plot(jseg$cnlr[is.finite(jseg$cnlr)], pch=".", cex=2, col = c("grey","lightblue","azure4","slateblue")[chrcol], ylab="log-ratio", xaxt="n")
    abline(h=median(jseg$cnlr, na.rm=TRUE), col="green4")
    abline(h=x$dipLogR, col="magenta4")
    segments(segstart, cnlr.median, segend, cnlr.median, lwd=1.75, col=2)
    if (missing(chromlevels)) { # human genome
        if (length(nn) == 23) { # chromsome X is present
            mtext(c(1:22,"X"), side=3, line=0, at=(nn+c(0,nn[-23]))/2, cex=0.65)
        } else { # IMPACT assay X not present
            mtext(1:22, side=3, line=0, at=(nn+c(0,nn[-22]))/2, cex=0.65)
        }
    } else {
        mtext(chromlevels, side=3, line=0, at=(nn+c(0,nn[-length(nn)]))/2, cex=0.65)
    }
    # plot the logOR data and mafR
    plot(jseg$valor[is.finite(jseg$cnlr)], pch=".", cex=2.5, col = c("grey","lightblue","azure4","slateblue")[chrcol], ylab="log-odds-ratio", ylim=c(-4,4), xaxt="n")
    segments(segstart, sqrt(mafR), segend, sqrt(mafR), lwd=1.75, col=2)
    segments(segstart, -sqrt(mafR), segend, -sqrt(mafR), lwd=1.75, col=2)
    axis(1)
    mtext(side=1, line=1.75, "Index", cex=0.8)
    par(def.par)  #- reset to default
}
