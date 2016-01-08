clustersegs <- function(out, jointseg, min.nhet=10) {
    # number of segments
    nsegs <- nrow(out)
    # logR with NAs removed
    cnlr <- jointseg$cnlr[is.finite(jointseg$cnlr)]
    # logOR and its variance for the het snps (and no NAs from cnlr)
    lor <- jointseg[is.finite(jointseg$cnlr), c("valor","lorvar","het")]
    # segid is the segment number re-ordered by cnlr.median
    segid <- rep(rank(out$cnlr.median, ties.method="random"), out$num.mark)
    # logR, logOR data re-ordered by re-ordering the segments
    cnlr <- cnlr[order(segid)]
    lor <- lor[order(segid),]
    segid <- sort(segid)
    # sort the out dataframe
    out <- out[order(out$cnlr.median),]
    # observed copy number (ocn) clusters
    ocnclust <- 1:nsegs
    ocnlevels <- out$cnlr.median
    while ((length(ocnlevels) > 1) && (min(diff(ocnlevels)) < 0.04)) {
        j <- which.min(diff(ocnlevels))
        ocnlevels <- ocnlevels[-j]
        ocnclust[ocnclust>j] <- ocnclust[ocnclust>j] - 1
        segid[segid > j] <- segid[segid > j] - 1
        ocnlevels[j] <- median(cnlr[segid==j])
    }
    out$segclust <- ocnclust
    out$cnlr.median.clust <- ocnlevels[ocnclust]
    # create the logOR data needed for merging by maf
    segid <- rep(1:nsegs, out$num.mark)[lor$het==1]
    lor <- lor[lor$het==1, c("valor","lorvar")]
    # loop through the ocnclust to cluster based on mafR
    mafR.clust <- mafclust <- rep(NA, nsegs)
    # function to estimate maf from clustered data
    maffun <- function(x) {
        # occasional extreme valor can screw maf. so winsorize the maf
        valor <- abs(x$valor)
        lorvar <- x$lorvar
        # extreme large values
        valor.thresh <- median(valor) + 3*sqrt(quantile(lorvar, 0.8, type=1))
        valor[valor > valor.thresh] <- valor.thresh
        # extreme small values (not likely)
        valor.thresh <- median(valor) - 3*sqrt(quantile(lorvar, 0.8, type=1))
        valor[valor < valor.thresh] <- valor.thresh
        sum(((valor)^2 - lorvar)/lorvar)/sum(1/lorvar)
    }
    # loop through ocn clusters
    for(i in unique(ocnclust)) {
        # segments with at least 10 hets and 2% hets (address male = single X)
        ii <- ocnclust==i & out$nhet >= min.nhet & out$nhet/out$num.mark > 0.01
        segs <- which(ii)
        # if more than one segment start merging from largest
        if (length(segs) > 1) {
            # order the segs by nhet
            segs <- segs[order(out$nhet[ii], decreasing=TRUE)]
            # list to hold the data for clusters
            lorclust <- list()
            # first seg starts a cluster
            mafclust[segs[1]] <- 1
            lorclust[[1]] <- lor[segid==segs[1],]
            # merge only to a cluster with maf closest to it
            # maf of first cluster
            cmaf <- maffun(lorclust[[1]])
            # loop through other segs
            nclust <- 1
            for (j in 2:length(segs)) {
                lorj <- lor[segid==segs[j],]
                mafj <- maffun(lorj)
                # cluster with lower & higher maf
                jlo <- which(cmaf < mafj)[which.max(cmaf[cmaf < mafj])]
                jhi <- which(cmaf > mafj)[which.min(cmaf[cmaf > mafj])]
                # Mann-Whitney p-value compared wrt existing clusters
                mwplo <- mwphi <- 0
                if (length(jlo) > 0) mwplo <- wilcox.test(abs(lorclust[[jlo]]$valor), abs(lorj$valor), exact=FALSE)$p.value
                if (length(jhi) > 0) mwphi <- wilcox.test(abs(lorclust[[jhi]]$valor), abs(lorj$valor), exact=FALSE)$p.value
                # if p-value is not small enough
                if (max(mwplo, mwphi) > 0.001) {
                    if (mwplo > mwphi) {
                        mafclust[segs[j]] <- jlo
                        lorclust[[jlo]] <- rbind(lorclust[[jlo]], lorj)
                        cmaf[jlo] <- maffun(lorclust[[jlo]])
                    } else {
                        mafclust[segs[j]] <- jhi
                        lorclust[[jhi]] <- rbind(lorclust[[jhi]], lorj)
                        cmaf[jhi] <- maffun(lorclust[[jhi]])
                    }
                } else {
                    nclust <- nclust+1
                    mafclust[segs[j]] <- nclust
                    lorclust[[nclust]] <- lorj
                    cmaf <- c(cmaf, mafj)
                }
            }
            mafR.clust[ii] <- cmaf[mafclust[ii]]
        } else {
            mafclust[ii] <- 1
            mafR.clust[ii] <- out$mafR[ii]
        }
    }
    # redo the segclust to include mafclust
    nclust <- 0
    segclust <- ocnclust
    for (i in unique(ocnclust)) {
        ii <- ocnclust==i
        segclust[ii] <- nclust + mafclust[ii]
        # number of clusters so far only if at least one non-NA mafclust
        if (sum(is.finite(mafclust[ii])) > 0) {
            nclust <- nclust + max(mafclust[ii], na.rm=TRUE)
        }
        # change the NA into new cluster
        if (sum(is.na(mafclust[ii])) > 0) {
            nclust <- nclust + 1
            segclust[ii][is.na(mafclust[ii])] <- nclust
        }
    }
    out$segclust <- segclust
    out$mafR.clust <- mafR.clust
    # order out back to the original genomic order
    out[order(out$seg),]
}
