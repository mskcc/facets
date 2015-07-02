# function to find logR level corresponding to the diploid state
# uses the clusters of segments of (logR, logOR) data with high cval
# we then look for clusters with allelic balance (as in normal diploid)
# allelic balance is also in k+k state (k copies of both parental chroms)
# we believe 0+0 can only happen for small focal regions and >2+2 rare
#
findDiploidLogR <- function(out, cnlr) {
    flags <- NULL
    # summarize the segment specific info into clusters
    num.mark <- tapply(out$num.mark, out$segclust, sum)
    segclust <- names(num.mark)
    nhet <- tapply(out$nhet, out$segclust, sum)
    cnlr.median <- tapply(out$cnlr.median.clust, out$segclust, median)
    mafR <- tapply(out$mafR.clust, out$segclust, median)
    out0 <- data.frame(segclust, num.mark, nhet, cnlr.median, mafR)
    # further clustering to get ocn states (e.g. 1+1 and 2+0 at same level)
    # set up the cnlr data and segid (segclust level)
    cnlr <- cnlr[is.finite(cnlr)]
    segid <- rep(out$segclust, out$num.mark)
    # re-order data by cnlr.medians
    segid <- rank(out0$cnlr.median, ties="random")[segid]
    cnlr <- cnlr[order(segid)]
    segid <- sort(segid)
    # re-order the clusters by cnlr.median
    out0 <- out0[order(out0$cnlr.median),]
    # the components need to be re-ordered too
    mafR <- out0$mafR
    cnlr.median <- out0$cnlr.median
    num.mark <- out0$num.mark
    segclust <- out0$segclust
    nhet <- out0$nhet
    # create new clustered medians by joining adjacent cnlr.median levels
    k <- nrow(out0)
    ocnclust <- 1:k
    ocnlevels <- out0$cnlr.median
    while (min(diff(ocnlevels)) < 0.05 & length(ocnlevels)>1) {
        j <- which.min(diff(ocnlevels))
        ocnlevels <- ocnlevels[-j]
        ocnclust[ocnclust>j] <- ocnclust[ocnclust>j] - 1
        segid[segid > j] <- segid[segid > j] - 1
        ocnlevels[j] <- median(cnlr[segid==j])
    }
    # add to out0 data frame
    out0$ocnclust <- ocnclust
    out0$ocnlevels <- ocnlevels[ocnclust]
    # revised ocn
    #
    # mafR < 0.02 allows for single copy change at no more than 15.2%
    # check if there are balanced clusters
    bsegs <- which(out0$mafR < 0.02)
    # if none exists use the cluster with smallest mafR
    if(length(bsegs) == 0) {
        bsegs <- which.min(mafR)
        flags <- "mafR not sufficiently small"
    }
    #
    # do the cnrl.median of balanced segs suggest two possible locations?
    if (length(bsegs) > 1) {
        # find the low and high values by clustering them
        cm2cls <- km2class(cnlr.median[bsegs], num.mark[bsegs])
        # if the centers are separated by < 0.1375 only one solution
        # 0.1375 = log2(2.2) - log2(2); 2.2 implies CN of 4 at 10% cf
        if (diff(cm2cls$centers) < 0.1375) {
            dipLogR <- median(rep(cnlr.median[bsegs], num.mark[bsegs]))
        } else {
            # check if the lower value represent 1+1 or 2+2?
            dipLogR <- cm2cls$centers
            # if the lower value represents small fraction of data can be 0+0
            if (sum(num.mark[bsegs][cm2cls$cluster==1])/sum(num.mark) < 0.01) {
                dipLogR <- cm2cls$centers[2]
            }
        }
    } else {
        dipLogR <- cnlr.median[bsegs]
    }
    # cat("dipLogR =", dipLogR, "\n")
    # print(out0)
    # first remove the balanced segs
    out1 <- out0[-bsegs,]
    # is dipLogR[1] the 1+1 level; if so all lower cnlr.med values are losses
    # since 2+0 will have same true cnlr as 1+1 check only clear losses
    # clear losses defined as ocn=1.85 (loss at 15% cf); log2(1.85/2)=-0.1125
    lsegs <- which(out1$cnlr.median <= dipLogR[1] - 0.1125 & is.finite(out1$mafR))
    not1plus1 <- FALSE
    if (length(lsegs) == 0) {
        # identifiability issues if none exists; call dipLogR[1] as 1+1
        dipLogR <- dipLogR[1]
    } else {
        # check if dipLogR[1] is 1+1 or 2+2
        out1 <- out1[lsegs,]
        # assume dipLogR[1] is clonal 2+2 get acn consistent with cnlr, mafR
        # search acn from 3+0, 2+0 & 1+0; as 2+1 looks like 2+0 with lower cf
        # 0+0 has balanced alleles and so won't be in lsegs
        out2 <- t(sapply(1:nrow(out1), function(i, dlr, out1) {
                             acnsplit(dlr, out1$cnlr.median[i], out1$mafR[i])
                         }, dipLogR[1], out1))
        # print(out2)
        out1$acn <- apply(out2[,2*(1:3), drop=FALSE], 1, which.min)
        # all segments 2+0 can be because of dipLogR[1] being 1+1
        ii <- which(out1$acn != 2)
        # if the clusters affected account for more than 2% of data
        if (sum(out1$num.mark[ii])/sum(num.mark) > 0.02) {
            # dipLogR[1] is not 1+1
            not1plus1 <- TRUE
        } else {
            # 2+2 --> {2+1, 2+0} indistinguishable from 1+1 --> 1+0
            # check that it's not 2+2 with 2+1 and 2+0 states only
            # ***** can only be done if more than one cluster exists *****
            # cf2 is cellular fraction of 2+0 segment (delta-logR from 2+2)
            # cf3 is cellular fraction of 2+1 calculated as if its is 2+0
            # then cf3 = cf2/(2+cf2)
            if (sum(out1$acn==2) > 1) {
                # get optimal cellular fraction 
                rho <- seq(0.01, 0.99, by=0.01)
                # deviance
                cf2 <- out2[out1$acn==2,3]
                nm <- out1$num.mark[out1$acn==2]
                dev <- sapply(rho, function(x, cf2, nm) {
                                  sum(nm*pmin((cf2-x)^2, (cf2-x/(2+x))^2))
                              }, cf2, nm)
                cf0 <- rho[which.min(dev)]
                # are both cf0 and cf0/(2+cf0) represented among segments
                acn <- sapply(cf2, function(x, cf0) {
                                  which.min((c(x,x)-c(cf0, cf0/(1+cf0)))^2)
                              }, cf0)
                #cat("acn =", acn, "\n")
                #cat("cf2 =", cf2, "\n")
                # are the two levels represented in sufficient proportions
                if (min(sum(nm[acn==1]), sum(nm[acn==2]))/sum(nm) > 0.25) {
                    not1plus1 <- TRUE
                    # if the discrepancy from cf large flag it
                    if (max(abs(cf2-c(cf0, cf0/(1+cf0))[acn])) > 0.1) {
                        flags <- c(flags, "could be polyclonal 1 copy loss")
                    }
                }
            }
        }
        #print(out1)
    }
    # if dipLogR[1] determined to be not 1+1 estimate CN=2 level
    if (not1plus1) {
        # find deviance for each ocnlevel
        out1 <- out0[out0$cnlr.median <= max(dipLogR) & is.finite(out0$mafR),]
        # ocn levels cannot be any lower than lr4-1
        ocnlevels0 <- ocnlevels[ocnlevels > dipLogR[1]-1 & ocnlevels < dipLogR[1]]
        dev1 <- sapply(ocnlevels0, dlrdev, dipLogR[1], out1)
        colr <- rep(1, length(dev1))
        if (length(dipLogR) == 2) {
            ocnlevels1 <- ocnlevels[ocnlevels > dipLogR[2]-1 & ocnlevels < dipLogR[2]]
            dev2 <- sapply(ocnlevels1, dlrdev, dipLogR[2], out1)
            ocnlevels0 <- c(ocnlevels0, ocnlevels1)
            dev1 <- c(dev1, dev2)
            colr <- c(colr, rep(2, length(dev2)))
        }
        #plot(ocnlevels0, dev1, pch=16, col=colr)
        #print(dev1)
        cn2logR <- ocnlevels0[which.min(dev1)]
        # return the estimate 1+1 dipLogR
    } else {
        cn2logR <- dipLogR[1]
    }
    names(cn2logR) <- NULL # remove names it may have acquired
    out0 <- out0[order(out0$segclust),] # reorder by segclust
    list(out0=out0, dipLogR=cn2logR, flags=flags)
}

# split the segclust medians into two groups high and low
km2class <- function(cnlr.med, num.mark) {
    n <- length(cnlr.med)
    ii <- order(cnlr.med)
    csx <- cumsum((cnlr.med*num.mark)[ii])
    csn <- cumsum(num.mark[ii])
    jj <- which.max(csx[-n]^2/csn[-n] + (csx[n]-csx[-n])^2/(csn[n]-csn[-n]))
    out <- list()
    out$centers <- c(csx[jj]/csn[jj], (csx[n]-csx[jj])/(csn[n]-csn[jj]))
    out$cluster <- rep(2,n)
    out$cluster[ii][1:jj] <- 1
    out
}


# redo cf using log-logOR
acnsplit <- function(lr0, lr1, lorsq) {
    deltalr <- (lr0-lr1)
    llor <- log(sqrt(lorsq))
    rho <- seq(0.01, 0.99, by=0.01)
    # log-logOR has a much wider range so weigh logR contribution more
    util3.0 <- 5*(deltalr-(log2(2+2*rho) - log2(2+rho)))^2 + (llor-log(log((1+2*rho)/(1-rho))))^2
    util2.0 <- 5*(deltalr-(log2(2+2*rho) - log2(2)))^2 + (llor-log(log((1+rho)/(1-rho))))^2
    util1.0 <- 5*(deltalr-(log2(2+2*rho) - log2(2-rho)))^2 + (llor-log(log(1/(1-rho))))^2
    optrho <- matrix(0,2,3)
    i <- which.min(util1.0); optrho[,1] <- c(rho[i], util1.0[i])
    i <- which.min(util2.0); optrho[,2] <- c(rho[i], util2.0[i])
    i <- which.min(util3.0); optrho[,3] <- c(rho[i], util3.0[i])
    rownames(optrho) <- c("cf","util")
    colnames(optrho) <- c("acn1.0","acn2.0","acn3.0")
    optrho
}

# deviance metric for candidate dipLogR
#  lr2 is candidate logR value at CN=2
#  lr4 is the balanced allele logR value (CN=4)
#  (lr, lor) are the logR logOR for the clusters
dlrdev <- function(lr2, lr4, out) {
    # use cnlr.median, mafR etc from out
    k <- nrow(out)
    lr <- out$cnlr.median
    # change mafR to log-odds ratio
    lor <- sqrt(pmax(0,out$mafR))
    # cellular fraction
    rho <- 2^(lr4-lr2) - 1
    # 1 and 3 copy logR levels
    lr1 <- lr2 - log2(2/(2-rho))
    lr3 <- lr2 + log2((2+rho)/2)
    dev <- sapply(1:k, function(i) {
               # closest copy number state
               cn <- which.min(abs(c(lr1,lr2,lr3,lr4) - lr[i]))
               # log-ratio deviance
               dev <- (lr[i] - c(lr1,lr2,lr3,lr4)[cn])^2
               # allele ratio
               if (cn == 1) dev <- dev + (lor[i] - log(1/(1-rho)))^2
               if (cn == 2) dev <- dev + min(lor[i]^2, (lor[i] - log((1+rho)/(1-rho)))^2)
               if (cn == 3) dev <- dev + min((lor[i]-log(1+rho))^2, (lor[i] - log((1+rho)/(1-rho)))^2)
               if (cn == 4) dev <- dev + min(lor[i]^2, (lor[i]-log(1+2*rho))^2, (lor[i] - log((1+3*rho)/(1-rho)))^2)
               dev
           })
    sum(dev*out$num.mark)
}

# wrapper for using findDiploidLogR and then clusteredcncf.fit
fitcncf <- function(out, cnlr) {
    oo <- findDiploidLogR(out, cnlr)
    out1 <- clusteredcncf.fit(oo$out0, oo$dipLogR)
    cncf <- out
    ii <- match(cncf$segclust, out1$segclust)
    cncf$cf <- out1$cf[ii]
    cncf$tcn <- out1$tcn[ii]
    cncf$lcn <- out1$lcn[ii]
    list(out=cncf, dipLogR=oo$dipLogR, flags=oo$flags)
}
