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
    segclust <- as.numeric(names(num.mark))
    nhet <- tapply(out$nhet, out$segclust, sum)
    cnlr.median <- tapply(out$cnlr.median.clust, out$segclust, median)
    mafR <- tapply(out$mafR.clust, out$segclust, median)
    out0 <- data.frame(segclust, num.mark, nhet, cnlr.median, mafR)
    # further clustering to get ocn states (e.g. 1+1 and 2+0 at same level)
    # set up the cnlr data and segid (segclust level)
    cnlr <- cnlr[is.finite(cnlr)]
    segid <- rep(out$segclust, out$num.mark)
    # re-order data by cnlr.medians
    segid <- rank(out0$cnlr.median, ties.method="random")[segid]
    cnlr <- cnlr[order(segid)]
    segid <- sort(segid)
    # re-order the clusters by cnlr.median
    out0 <- out0[order(out0$cnlr.median),]
    # the components need to be re-ordered too
    mafR <- out0$mafR
    cnlr.median <- out0$cnlr.median
    num.mark <- out0$num.mark
    nsnps <- sum(num.mark)
    segclust <- out0$segclust
    nhet <- out0$nhet
    # 
    ocnlevels <- unique(out0$cnlr.median)
    # revised ocn
    #
    # mafR < 0.025 allows for single copy change at no more than 17.13%
    # check if there are balanced clusters
    mafR.thresh <- 0.025
    bsegs <- which(out0$mafR < mafR.thresh)
    # if none exists of if balanced segs span less than 10% of genome
    # use clusters with smallest mafR
    if (sum(num.mark[bsegs])/nsnps < 0.1) {
        # set flag
        flags <- "mafR not sufficiently small"
        # get bsegs with mafR < 0.05; single copy change at 25%
        mafR.thresh <- 0.05
        bsegs <- which(out0$mafR < mafR.thresh)
        # if still not 10% use mafR.thresh of 0.09 (single copy gain at 35%)
        if (sum(num.mark[bsegs])/nsnps < 0.1) {
            mafR.thresh <- 0.09
            bsegs <- which(out0$mafR < mafR.thresh)
            if (sum(num.mark[bsegs])/nsnps < 0.1) {
                flags <- c(flags, "mafR<0.09 in less than 10% genome")
            }
        }
    }
    #
    # do the cnrl.median of balanced segs suggest two possible locations?
    if (length(bsegs) > 1) {
        # find the low and high values by clustering them
        cm2cls <- km2class(cnlr.median[bsegs], num.mark[bsegs])
        # if the centers are separated by < 0.1375 only one solution
        # 0.2016 = log2(2.3) - log2(2); 2.3 implies CN of 4 at 15% cf
        if (diff(cm2cls$centers) < 0.2) {
            dipLogR <- median(rep(cnlr.median[bsegs], num.mark[bsegs]))
            nbal <- sum(cm2cls$nbal)
        } else {
            # check if the lower value represent 1+1 or 2+2?
            dipLogR <- cm2cls$centers
            nbal <- cm2cls$nbal
            # if the lower value represents small fraction of data can be 0+0
            if (nbal[1]/nsnps < 0.01) {
                dipLogR <- cm2cls$centers[2]
                nbal <- cm2cls$nbal[2]
            }
        }
    } else {
        # make sure bsegs is not empty
        if (length(bsegs) == 0) {
            # if no balanced segs set dipLogR at the median of cnlr
            dipLogR <- median(cnlr)
            nbal <- 0
        } else {
            dipLogR <- cnlr.median[bsegs]
            nbal <- num.mark[bsegs]
        }
    }
    names(dipLogR) <- NULL
    names(nbal) <- NULL
    # flag if dipLogR has two values separated by < log2(1.25) (CN=4 at 25%)
    if (length(dipLogR) == 2)
        if (diff(dipLogR) < log2(1.25)) 
            flags <- c(flags, "possibly subclonal 2+2 states present")
    # cat("dipLogR =", dipLogR, "\n")
    # print(out0)
    # first remove the balanced segs
    out1 <- out0[-bsegs,]
    # is dipLogR[1] the 1+1 level; if so all lower cnlr.med values are losses
    # since 2+0 will have same true cnlr as 1+1 check only clear losses
    # clear losses defined as ocn=1.85 (loss at 15% cf); log2(1.85/2)=-0.1125
    lsegs <- which(out1$cnlr.median <= dipLogR[1] - 0.1125 & is.finite(out1$mafR))
    not1plus1 <- FALSE # indicator whether dipLogR[1] is not 1+1 state
    wgd.likely <- NULL # placeholder variable for ambiguos wgd call
    # identifiability issues when lsegs empty; call dipLogR[1] as 1+1
    if (length(lsegs) > 0) {
        # check if dipLogR[1] is 1+1 or 2+2
        out1 <- out1[lsegs,]
        # assume dipLogR[1] is 1+1. compute cf using logR and logOR data
        cflr <- pmin(2 - 2^(1 + out1$cnlr.median - dipLogR[1]), 1)
        cflor <- 1 - exp(-sqrt(pmax(0,out1$mafR)))
        # if segments where cflor > cflr + 0.1 has >5% of snps
        # cflor > cflr+0.1 because mafR can be low with 1+0 & 2+0 mixture
        if (sum(out1$num.mark[cflor > cflr+0.1])/nsnps > 0.05) {
            not1plus1 <- TRUE
            flags <- c(flags, paste("mafR larger than expected if", dipLogR[1], "is diploid level"))
        }
        # assume dipLogR[1] is clonal 2+2 get acn consistent with cnlr, mafR
        # search acn from 3+0, 2+0 & 1+0; as 2+1 looks like 2+0 with lower cf
        # 0+0 has balanced alleles and so won't be in lsegs
        out2 <- t(sapply(1:nrow(out1), function(i, dlr, out1) {
                             acnsplit(dlr, out1$cnlr.median[i], out1$mafR[i])
                         }, dipLogR[1], out1))
        # print(out2)
        out1$acn <- apply(out2[,2*(1:3), drop=FALSE], 1, which.min)
        # proportion of genome that fits 1+0, 2+0 & 3+0
        acn1prop <- sum(out1$num.mark[out1$acn == 1])/nsnps
        acn2prop <- sum(out1$num.mark[out1$acn == 2])/nsnps
        acn3prop <- sum(out1$num.mark[out1$acn == 3])/nsnps
        # mix of 2+0 and 1+0 from 1+1 (on same segs) can mimic 3+0 from 2+2
        # that is, it looks like a loss but has a large mafR
        # so allow small % if 1+0 and a slightly bigger % of 3+0
        if (!not1plus1 & acn1prop < 0.005 & acn3prop < 0.05) {
            # if acn3prop > 0 flag it
            if (acn3prop > 0) 
                flags <- c(flags, paste("likely mixture of 1+0 & 2+0 in segclust:", paste(out1$segclust[out1$acn==3], collapse=", ")))
            # most of the losses are either 2+0 or 2+1
            # 2+2 --> {2+1, 2+0} indistinguishable from 1+1 --> 1+0
            # check that it's not 2+2 with 2+1 and 2+0 states only
            # ***** can only be done if more than one cluster exists *****
            # cf2 is cellular fraction of 2+0 segment (delta-logR from 2+2)
            # cf3 is cellular fraction of 2+1 calculated as if its is 2+0
            # then cf3 = cf2/(2+cf2)
            # need 2 or more clusters with 2+0 call; otherwise it is 1+1
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
                                  which.min((c(x,x)-c(cf0, cf0/(2+cf0)))^2)
                              }, cf0)
                #cat("acn =", acn, "\n")
                #cat("cf2 =", cf2, "\n")
                # are the two levels represented in sufficient proportions
                nm2 <- sum(nm[acn==1]) # this is the 2+0 fit
                nm3 <- sum(nm[acn==2]) # this is the 2+1 fit
                if (min(nm3, nm2)/(nm3+nm2) > 0.25) {
                    # if the discrepancy from cf large call it polyclonal
                    if (max(abs(cf2-c(cf0, cf0/(2+cf0))[acn])) > 0.05) {
                        flags <- c(flags, "polyclonal 1 copy loss fits better")
                    } else {
                        # polyclonal 1 copy loss not obvious & wgd ambiguity
                        # so calculate alternate dipLogR and present both
                        not1plus1 <- TRUE
                        flags <- c(flags, "alternate dipLogR is possible")
                        # acn2prop >1/3 and 2+1 as likely as 2+0 call it WGD
                        if (acn2prop > 1/3 & nm3 > nm2) {
                            wgd.likely <- TRUE
                        } else {
                            wgd.likely <- FALSE
                        }
                    }
                }
            }
        } else {
            not1plus1 <- TRUE
        }
        if (not1plus1 & acn1prop+acn3prop > 0) {
            flags <- c(flags, paste("not consistent for 1 copy loss from diploid in segclust:", paste(out1$segclust[out1$acn!=2], collapse=", ")))
        }
    }
    # if dipLogR[1] is possibly not 1+1 estimate CN=2 level
    if (not1plus1) {
        # find deviance for each ocnlevel
        out1 <- out0[out0$cnlr.median <= max(dipLogR) & is.finite(out0$mafR),]
        # ocn levels cannot be any lower than lr4-1
        ocnlevels0 <- ocnlevels[ocnlevels > dipLogR[1]-1 & ocnlevels < dipLogR[1]]
        dev1 <- sapply(ocnlevels0, dlrdev, dipLogR[1], out1)
        #colr <- rep(1, length(dev1))
        if (length(dipLogR) == 2) {
            # if dipLogR[1] has >2 (1+1 & 2+0) segclust at >3% don't do this
            if (sum(out0$num.mark[out0$cnlr.median==dipLogR[1]]/nsnps > 0.03) > 2) {
                flags <- c(flags, "multiple mafR levels at bal allele level 1")
            } else {
                ocnlevels1 <- ocnlevels[ocnlevels > dipLogR[2]-1 & ocnlevels < dipLogR[2]]
                dev2 <- sapply(ocnlevels1, dlrdev, dipLogR[2], out1)
                ocnlevels0 <- c(ocnlevels0, ocnlevels1)
                dev1 <- c(dev1, dev2)
                #colr <- c(colr, rep(2, length(dev2)))
            }
        }
        #plot(ocnlevels0, dev1, pch=16, col=colr)
        #print(dev1)
        cn2logR <- ocnlevels0[which.min(dev1)]
        # if wgd.likely is non-null add dipLogR[1]
        if (!is.null(wgd.likely)) {
            # if wgd.likely set dipLogR[1] as altDipLogR
            altDipLogR <- dipLogR[1]
            # o/w change altDipLogR to cn2logR to dipLogR[1] & 
            if (!wgd.likely) {
                altDipLogR <- cn2logR
                cn2logR <- dipLogR[1]
            }
            names(altDipLogR) <- NULL # remove names it may have acquired
        }
    } else {
        cn2logR <- dipLogR[1]
    }
    names(cn2logR) <- NULL    # remove names it may have acquired
    # 2+2 in more than half the genome seems a stretch
    if (any(nbal/nsnps > 0.5)) {
        i <- which(nbal/nsnps > 0.5)
        # if chosen dipLogR is below that level
        if (cn2logR < dipLogR[i] - 0.05) {
            flags <- c(flags, ">1/2 genome balanced at chosen dipLogR")
            if (exists("altDipLogR")) {
                altDipLogR <- c(altDipLogR, cn2logR)
            } else {
                altDipLogR <- cn2logR
                flags <- c(flags, "alternate dipLogR is possible")
            }
            cn2logR <- dipLogR[i]
        }
    }
    # alternate check if "mafR not sufficiently small"
    if ("mafR not sufficiently small" %in% flags) {
        pbal <- sapply(dipLogR, function(x, y, n) {
                           sum(abs(y-x) < 0.1375)/n
                       }, rep(cnlr.median, num.mark), nsnps)
        i <- which(pbal > 2/3)
        if (length(i) > 0) {
            if (cn2logR < dipLogR[i] - 0.05) {
                flags <- c(flags, ">2/3 genome balanced at chosen dipLogR")
                if (exists("altDipLogR")) {
                    altDipLogR <- c(altDipLogR, cn2logR)
                } else {
                    altDipLogR <- cn2logR
                    flags <- c(flags, "alternate dipLogR is possible")
                }
            }
            cn2logR <- dipLogR[i]
        }
    }
    
    out0 <- out0[order(out0$segclust),] # reorder by segclust
    if (is.null(wgd.likely)) {
        list(out0=out0, dipLogR=cn2logR, alBalLogR=cbind(dipLogR,nbal/nsnps), mafR.thresh=mafR.thresh, flags=flags)
    } else {
        list(out0=out0, dipLogR=cn2logR, altDipLogR=altDipLogR, alBalLogR=cbind(dipLogR,nbal/nsnps), mafR.thresh=mafR.thresh, flags=flags)
    }
}

# split the segclust medians into two groups high and low
km2class <- function(cnlr.med, num.mark) {
    n <- length(cnlr.med)
    ii <- order(cnlr.med)
    csx <- rep(cnlr.med[ii], num.mark[ii])
    iseg <- rep(1:n, num.mark[ii])
    gg <- rep(1:2, c(1,n-1))[iseg]
    l1norm <- sapply(1:(n-1), function(i, n, csx, iseg) {
                         gg <- rep(1:2, c(i,n-i))[iseg]
                         med2 <- tapply(csx, gg, median)
                         sum(abs(csx - med2[gg]))
                     }, n, csx, iseg)
    j0 <- which.min(l1norm)
    med2 <- tapply(csx, rep(1:2, c(j0,n-j0))[iseg], median)
    out <- list()
    out$centers <- med2
    out$cluster <- rep(2,n)
    out$cluster[ii][1:j0] <- 1
    out$nbal <- tapply(num.mark, out$cluster, sum)
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
    # if optrho is low weigh the utilities in favor of 2
    if (optrho[1,1] < 0.25) optrho[2,] <- optrho[2,]*(1+c(50,0,50)*(0.25-optrho[1,1]))
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
