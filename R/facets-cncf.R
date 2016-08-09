# fitting copy number and cellular fractions (replaces clusteredcncf)
fitcncf <- function(out, dipLogR=0, nX=23) {
    # save the original out
    cncf <- out
    # get the segclust version of out
    num.mark <- tapply(out$num.mark, out$segclust, sum)
    segclust <- as.numeric(names(num.mark))
    nhet <- tapply(out$nhet, out$segclust, sum)
    cnlr.median <- tapply(out$cnlr.median.clust, out$segclust, median)
    mafR <- tapply(out$mafR.clust, out$segclust, median)
    out0 <- data.frame(segclust, num.mark, nhet, cnlr.median, mafR)
    # fit the copy numbers and cellular fraction
    out1 <- fitcncf0(out0, dipLogR)
    ii <- match(cncf$segclust, out1$segclust)
    cncf$cf <- out1$cf[ii]
    cncf$tcn <- out1$tcn[ii]
    cncf$lcn <- out1$lcn[ii]
    # revise the copy numbers for X if subject is male
    ii <- which(cncf$chrom==nX)
    if (length(ii) > 0) {
        # male if nhets on X chromosome is less than 1% of num.mark
        if (sum(cncf$nhet[ii])/sum(cncf$num.mark[ii]) < 0.01) {
            cncf$tcn[ii] <- ceiling(cncf$tcn[ii]/2)
            cncf$lcn[ii] <- 0
            # set cf to be 1 if tcn==1
            cncf$cf[ii][cncf$tcn[ii]==1] <- 1
        }
    }
    # if tcn is 0 or 1 lcn has to be zero
    cncf$lcn[cncf$tcn<=1] <- 0
    # return result
    cncf
}

fitcncf0 <- function(out, dipLogR=0) {
    # initialize vectors
    out$lcn <- out$tcn <- out$cf <- out$ocn <- rep(NA_real_, nrow(out))
    out$ocn <- 2^(1 + out$cnlr.median - dipLogR)
    ocn0 <- out$ocn
    maf0 <- 1/(1+exp(sqrt(pmax(0,out$mafR))))
    # loop through the segment clusters
    for(seg in 1:nrow(out)) {
        ocn <- ocn0[seg]
        maf <- maf0[seg]
        if (is.finite(maf)) {
            if (ocn > 2.15) {
                tcn <- ceiling(ocn)
                uu <- optcfutil(tcn, ocn, maf)
                # it distance measure is large - 0.0225 when diff CN of 0.15
                # make sure the fitted cf is not >0.99
                if (uu[3] > 0.0225 | uu[2]>0.99) {
                    tcn1 <- tcn+1
                    uu1 <- optcfutil(tcn1, ocn, maf)
                    # if distance measure is reduced by at least 20%
                    if (uu1[3] < 0.8*uu[3]) {
                        # if distance is still large - 0.01 when diff CN of 0.1
                        if (uu1[3] > 0.01) {
                            tcn2 <- tcn+2
                            uu2 <- optcfutil(tcn2, ocn, maf)
                            if (uu2[3] < 0.8*uu1[3]) {
                                uu1 <- uu2
                                tcn1 <- tcn2
                            }
                        }
                        uu <- uu1
                        tcn <- tcn1
                    }
                }
                out$tcn[seg] <- tcn
                out$lcn[seg] <- uu[1]
                out$cf[seg] <- uu[2]
            } else {
                if (ocn < 1.85) {
                    # choose between 0 & 1
                    uu0 <- optcfutil(0, ocn, maf)
                    uu1 <- optcfutil(1, ocn, maf)
                    if (uu0[3] < uu1[3]) {
                        out$tcn[seg] <- 0
                        out$lcn[seg] <- uu0[1]
                        out$cf[seg] <- uu0[2]
                    } else {
                        out$tcn[seg] <- 1
                        out$lcn[seg] <- uu1[1]
                        out$cf[seg] <- uu1[2]
                    }
                } else {
                    # ocn between 1.85 & 2.15; choose between 1+1 & 2+0
                    uu2 <- optcfutil(2, ocn, maf)
                    # if optimal cf is > 15% call it 2+0 o.w 1+1
                    out$tcn[seg] <- 2
                    if (uu2[2] > 0.15) {
                        out$lcn[seg] <- uu2[1]
                        out$cf[seg] <- uu2[2]
                    } else {
                        out$lcn[seg] <- 1
                        out$cf[seg] <- 1
                    }
                }
            }
        }
    }
    # merge thc cf levels
    out <- mergecf(out)
    # now fill in the tcn for clusters with NA for mafR using 3rd quartile cf
    ii <- which(out$cf < 1)
    # don't want to call deletions and amplications in low purity samples
    if (length(ii) > 0) {
        cf3q <- max(quantile(rep(out$cf[ii], out$num.mark[ii]), 0.75), 0.3)
    } else {
        cf3q <- 0.3
    }
    for(seg in 1:nrow(out)) {
        if (is.na(maf0[seg])) {
            ocn <- ocn0[seg]
            # divide ocn by cf3q and round it
            tcn <- round((ocn-2)/cf3q) + 2
            if (tcn < 0) tcn <- 0
            out$tcn[seg] <- tcn
            # set the segment cf to be cf3q
            out$cf[seg] <- cf3q
        }
    }
    out
}

# minimizes the distance(1+l*cf - ocn1)^2 + (1+k*cf - ocn2)^2 to get "cf"
# where ocn1 = ocn*maf, ocn2=ocn*(1-maf), l <= k true copy numbers in tumor
optcfutil <- function(tcn, ocn, maf) {
    cn1 <- ocn*maf
    cn2 <- ocn*(1-maf)
    maxlcn <- floor(tcn/2)
    oo <- sapply(0:maxlcn, function(l, tcn, cn1, cn2) {
                     k <- tcn - l
                     if (k==1 & l==1) {
                         cf = 1
                     } else {
                         cf = ((l-1)*(cn1-1) + (k-1)*(cn2-1))/((l-1)^2+(k-1)^2)
                     }
                     # if the optimal cf outside (0,1) set it at boundary
                     if (cf < 0) cf <- 0
                     if (cf > 1) cf <- 1
                     # calculate utility
                     util <- (1+(l-1)*cf - cn1)^2 + (1+(k-1)*cf - cn2)^2
                     c(l,cf,util)
                 }, tcn, cn1, cn2)
    oo[,which.min(oo[3,])]
}

# merging cellular fractions
mergecf <- function(out) {
    # clusters with defined cf and less than 1
    out0 <- out[which(out$cf < 1),]
    if (nrow(out0) > 0) { 
        # order out0 by cf
        out0 <- out0[order(out0$cf),]
        # prepare the variables needed for deviance
        maf <- 1/(1+exp(sqrt(pmax(out0$mafR,0))))
        acn1o <- out0$ocn*maf
        acn2o <- out0$ocn*(1-maf)
        acn1t <- out0$lcn - 1
        acn2t <- out0$tcn-out0$lcn -1
        cf <- out0$cf
        nm <- out0$num.mark
        # merge cf values closest to one another
        cfclust0 <- cfclust <- 1:nrow(out0)
        cflevels0 <- cflevels <- cf
        dev1 <- dev0 <- sum(((1+cf*acn1t - acn1o)^2 + (1+cf*acn2t - acn2o)^2)*nm)
        while (length(cflevels)>1 & dev1 <= 1.1*dev0) {
            # save the previous levels
            cflevels0 <- cflevels
            cfclust0 <- cfclust
            # dev1a <- dev1
            # update
            j <- which.min(diff(cflevels))
            cflevels <- cflevels[-j]
            cfclust[cfclust>j] <- cfclust[cfclust>j] - 1
            cflevels[j] <- sum((nm*cf)[cfclust==j])/sum(nm[cfclust==j])
            dev1 <- sum(((1+cflevels[cfclust]*acn1t - acn1o)^2 + (1+cflevels[cfclust]*acn2t - acn2o)^2)*nm)
        }
        # replace the existing cf with merged cf
        ii <- match(out0$segclust, out$segclust)
        out$cf[ii] <- cflevels0[cfclust0]
    }
    out
}
