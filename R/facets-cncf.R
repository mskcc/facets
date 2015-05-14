clusteredcncf <- function(out, dipLogR=0) {
    out1 <- out[,7:9]
    names(out1)[2:3] <- c("cnlr.median", "mafR")
    out1 <- out1[!duplicated(out1$segclust), , drop=FALSE]
    out1 <- clusteredcncf.fit(out1, dipLogR)
    cncf <- out
    ii <- match(cncf$segclust, out1$segclust)
    cncf$cf <- out1$cf[ii]
    cncf$tcn <- out1$tcn[ii]
    cncf$lcn <- out1$lcn[ii]
    cncf
}

clusteredcncf.fit <- function(out, dipLogR=0) {
    out$lcn <- out$tcn <- out$cf <- out$ocn <- out$util <- rep(0, nrow(out))
    # setup functions for optimization
    utility <- function(cf, tcn, lcn, ocn, maf) {
        ((tcn-lcn)*cf + (1-cf) - ocn*(1-maf))^2 + (lcn*cf + (1-cf) - ocn*maf)^2
    }
    optcf <- function(tcn, lcn, ocn, maf) {
        ocf <- ifelse (tcn == 2 & lcn == 1,
                        1, 
                        (ocn*(tcn-lcn-1) + ocn*maf*(2*lcn-tcn) - (tcn-2))/((tcn-lcn-1)^2 + (lcn-1)^2))
        # if optimal cf is less than 5% set it to 5%
        if (ocf < 0) ocf <- 0
        if (ocf > 1) ocf <- 1
        c(ocf, utility(ocf, tcn, lcn, ocn, maf))
    }
    # loop through the segment clusters
    for(seg in 1:nrow(out)) {
        # observed copy number
        ocn <- 2^(1 + out$cnlr.median[seg] - dipLogR)
        out$ocn[seg] <- ocn
        # estimate copy number ratio from log odds ratio estimate
        if (is.finite(out$mafR[seg])) {
            maf <- 1/(1+exp(sqrt(max(0,out$mafR[seg]))))
            if (ocn > 2.25) {
                # if ocn is large choose the next integer (greedy)
                # initialize genotype vectors
                tcn <- ceiling(ocn)
                # obtain the optimal cf and corresponding utility
                maxlcn <- floor(tcn/2)
                zzz <- sapply(0:maxlcn, function(i) {optcf(tcn, i, ocn, maf)})
                icf <- which.min(zzz[2,])
                # if utility for closest tcn is too high bump tcn to tcn+1
                if (zzz[2, icf] > 0.1) {
                    # tcn+1
                    tcn1 <- tcn+1
                    maxlcn <- floor(tcn1/2)
                    zzz1 <- sapply(0:maxlcn, function(i) {optcf(tcn1, i, ocn, maf)})
                    icf1 <- which.min(zzz1[2,])
                    # if the utility is reduced 
                    if (zzz1[2, icf1] < zzz[2, icf]) {
                        zzz <- zzz1
                        icf <- icf1
                        tcn <- tcn1
                    }
                }
                # print(c(zzz1[,icf1],zzz2[,icf2]))
                out$cf[seg] <- zzz[1,icf]
                out$tcn[seg] <- tcn
                out$lcn[seg] <- icf-1
                out$util[seg] <- zzz[2,icf]
            } else {
                if (ocn < 1.9) {
                    # if ocn is small choose between 0 and 1
                    tcn <- c(0,1)
                    zzz <- sapply(1:2, function(i) {optcf(tcn[i], 0, ocn, maf)})
                    icf <- which.min(zzz[2,])
                    out$cf[seg] <- zzz[1,icf]
                    out$tcn[seg] <- tcn[icf]
                    out$lcn[seg] <- 0
                    out$util[seg] <- zzz[2,icf]
                } else {
                    # if ocn close to 2 choose between 1, 2 & 3
                    tcn <- c(1,2,2,3,3)
                    lcn <- c(0,0,1,0,1)
                    zzz <- sapply(1:5, function(i) {optcf(tcn[i], lcn[i], ocn, maf)})
                    icf <- which.min(zzz[2,])
                    out$util[seg] <- zzz[2,icf]
                    if (zzz[1,icf] < 0.05 | out$mafR[seg] < 0.015) {
                        # if optimal cellular fraction is < 5% or low allelic 
                        # imbalance set it as normal diploid
                        out$cf[seg] <- 1
                        out$tcn[seg] <- 2
                        out$lcn[seg] <- 1
                    } else {
                        out$cf[seg] <- zzz[1,icf]
                        out$tcn[seg] <- tcn[icf]
                        out$lcn[seg] <- lcn[icf]
                    }
                }
            }
        } else {
            # mafR is NA so just calculate TCN greedily if ocn <1.75 or >2.25
            ocn2 <- ocn-2
            if (abs(ocn2) > 0.25) {
                ocn2 <- ocn-2
                out$tcn[seg] <- sign(ocn2)*ceiling(abs(ocn2)) + 2
                out$cf[seg] <- ocn2/(out$tcn[seg]-2)
                out$lcn[seg] <- NA
            } else {
                out$cf[seg] <- 1
                out$tcn[seg] <- 2
                out$lcn[seg] <- NA
            }
        }

    }
    out
}
