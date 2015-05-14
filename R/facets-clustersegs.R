clustersegs <- function(jointseg, out, min.nhet=25, cval=35) {
    # data for further comparisons
    jseg <- jointseg[is.finite(jointseg$cnlr),c(1:2,11:14)]
    # create an easily ordered segment indicator
    jseg$segs <- jseg$segs + 1000*jseg$chrom
    # segment summary
    out$seg <- out$seg + 1000*out$chr
    # initialize cluster indicator 
    snpclust <- rep(0, nrow(jseg))    # snp level
    segclust <- rep(0, nrow(out))     # segment level
    # big segments that have sufficient number of het snps
    bigsegs <- out$nhet >= min.nhet
    if (sum(1*bigsegs) == 0) warning(paste("No segment with at least", min.nhet, "hets; try reducing min.nhet"))
    # order the segments
    ii <- which(bigsegs)[order(-out$num.mark[bigsegs])]
    nbseg <- length(ii)
    # need to be careful if length(ii) is 0 or 1
    nclust <- 0
    if (nbseg > 0) {
        segclust1 <- rep(0, nbseg)   # segments with sufficient het count
        # cluster them using both cnlr and valor
        # begin with the biggest segment as a cluster of its own
        nclust <- nclust + 1
        ijj <- jseg$segs == out$seg[ii[1]]
        snpclust[ijj] <- 1
        segclust1[1] <- 1
        # collect data from segments into cluster specific lists
        # copy number log-ratio (pre-sorted)
        ccnlr <- list()
        ccnlr[[1]] <- sort(jseg$cnlr[ijj])
        # log variant allele odds-ratio (pre-sorted)
        finitevalor <- is.finite(jseg$valor)
        cvalor <- list()
        cvalor[[1]] <- sort(abs(jseg$valor[ijj & finitevalor]))
        # now loop through other segments
        if (nbseg > 1) {
            for (i in 2:nbseg) {
                # segment data (pre-sorted)
                ijj <- jseg$segs == out$seg[ii[i]]
                icnlr <- sort(jseg$cnlr[ijj])
                ivalor <- sort(abs(jseg$valor[ijj & finitevalor]))
                # log-ratio wilcoxon statistic
                stat1 <- lapply(ccnlr, function(x, y) {
                    mwstat(x,y)
                }, icnlr)
                # log odds ratio wilcoxon statistic
                stat2 <- lapply(cvalor, function(x, y) {
                    mwstat(x,y)
                }, ivalor)
                # T squared
                tstat <- unlist(stat1) + unlist(stat2)
                # minimum value over all clusters
                mintstat <- min(tstat)
                if (mintstat < cval) {
                    # if minimal value < cval add to cluster with smallest T
                    imin <- which(tstat == mintstat)
                    snpclust[ijj] <- imin
                    ccnlr[[imin]] <- sort(c(ccnlr[[imin]], icnlr))
                    cvalor[[imin]] <- sort(c(cvalor[[imin]], ivalor))
                    segclust1[i] <- imin
                } else {
                    # else create a new cluster
                    nclust <- nclust + 1
                    snpclust[ijj] <- nclust
                    ccnlr[[nclust]] <- icnlr
                    cvalor[[nclust]] <- ivalor
                    segclust1[i] <- nclust
                }
            }
        }
        # put the big segment indicator in it
        segclust[ii] <- segclust1
    }
    # now cluster the other segments using cnlr alone
    ii <- which(!bigsegs)[order(-out$num.mark[!bigsegs])]
    # need to be careful if length(ii) is 0 or 1
    noseg <- length(ii)
    nclust1 <- 0
    if (noseg > 0) {
        segclust2 <- rep(0, noseg)   # segments with sufficient het count
        # begin with the biggest segment as a cluster of its own
        nclust1 <- nclust1 + 1
        ijj <- jseg$segs == out$seg[ii[1]]
        snpclust[ijj] <- nclust + 1
        segclust2[1] <- nclust + 1
        # collect data from segments into cluster specific lists
        # copy number log-ratio
        ccnlr0 <- list()
        ccnlr0[[1]] <- sort(jseg$cnlr[ijj])
        # now loop through other segments
        if (noseg > 1) {
            for (i in 2:noseg) {
                # segment data
                ijj <- jseg$segs == out$seg[ii[i]]
                icnlr <- sort(jseg$cnlr[ijj])
                # log-ratio wilcoxon statistic
                stat1 <- lapply(ccnlr0, function(x, y) {
                    mwstat(x,y)
                }, icnlr)
                # T squared
                tstat <- unlist(stat1)
                # minimum value over all clusters
                mintstat <- min(tstat)
                if (mintstat < cval) {
                    # if minimal value < cval add to cluster with smallest T
                    imin <- which(tstat == mintstat)
                    snpclust[ijj] <- nclust + imin
                    ccnlr0[[imin]] <- sort(c(ccnlr0[[imin]], icnlr))
                    segclust2[i] <- nclust + imin
                } else {
                    # else create a new cluster
                    nclust1 <- nclust1 + 1
                    snpclust[ijj] <- nclust + nclust1
                    ccnlr0[[nclust1]] <- icnlr
                    segclust2[i] <- nclust + nclust1
                }
            }
        }
        # put the big segment indicator in it
        segclust[ii] <- segclust2
    }
    # summarize cluster specific logR and logOR
    num.mark <- tapply(out$num.mark, segclust, sum)
    nhet <- tapply(out$nhet, segclust, sum)
    cnlr.median <- tapply(jseg$cnlr, snpclust, median)
    mafR <- tapply(1:nrow(jseg), snpclust, function(ii, valor, lorvar) {
        sum(((valor[ii])^2 - lorvar[ii])/lorvar[ii], na.rm=T)/sum(1/lorvar[ii], na.rm=T)
    }, jseg$valor, jseg$lorvar)
    # mafR only makes sense for clusters made from segs with sufficient hets
    if (nclust1 > 0) mafR[(nclust+1):(nclust+nclust1)] <- NA
    segclustsummary <- cbind(num.mark, nhet, cnlr.median, mafR)
    # return snpclust, segclust and segclustsummary
    list(snpclust=snpclust, segclust=segclust, segclustsummary=segclustsummary)
}

# Mann-Whitney statistic when x and y are already sorted
mwstat <- function(x,y) {
   n <- length(x)
   m <- length(y)
   zzz <- .Fortran("mwstat",
                   as.double(x),
                   as.integer(n),
                   as.double(y),
                   as.integer(m),
                   ustat=double(1))
   zzz$ustat
}

# mergw two sorted vectors to get a new sorted vector (to replace sort(c(x,y)))
mergexy <- function(x,y) {
   n <- length(x)
   m <- length(y)
   zzz <- .Fortran("mergexy",
                   as.double(x),
                   as.integer(n),
                   as.double(y),
                   as.integer(m),
                   xy=double(m+n),
                   as.integer(m+n))
   zzz$xy
}
