# mle of maf per segment
mafmle <- function(fcount, vaf) {
    nn <- fcount
    xx <- round(vaf*nn)
    plo <- 0.01
    phi <- 0.5
    # if nn is too large log-likelihood can be -Inf for some p
    llk01 <- sum(log(dbinom(xx, nn, plo) + dbinom(xx, nn, 1-plo)))
    llk50 <- sum(log(dbinom(xx, nn, phi) + dbinom(xx, nn, 1-phi)))
    if (!is.finite(llk01) | !is.finite(llk50)) {
        pp <- seq(0.01, 0.5, by=0.005)
        llk <- sapply(pp, function(x, xx, nn) {sum(log(dbinom(xx, nn, x) + dbinom(xx, nn, 1-x)))}, xx, nn)
        kk <- range(which(is.finite(llk)))
        plo <- pp[kk[1]]
        phi <- pp[kk[2]]
    }
    while(phi - plo > 0.0001) {
        pmid <- (plo+phi)/2
        llk <- sapply(pmid+c(-0.00005,0.00005), function(x, xx, nn) {sum(log(dbinom(xx, nn, x) + dbinom(xx, nn, 1-x)))}, xx, nn)
        ifelse(diff(llk) > 0, plo <- pmid, phi <- pmid)
    }
    c(pmid, mean(pmin(xx,nn-xx)/nn))
}

# segment summary
jointsegsummary <- function(jointseg) {
    # remove snps with NA in segs (due to NA in cnlr)
    jointseg <- jointseg[is.finite(jointseg$seg),]
    # initialize output table
    nsegs <- max(jointseg$seg)
    out <- as.data.frame(matrix(0, nsegs, 6))
    names(out) <- c("chrom","seg","num.mark","nhet","cnlr.median","mafR")
    # loop over the segments
    for(seg in 1:nsegs) {
        zz <- jointseg[jointseg$seg==seg, c("chrom","cnlr","het","valor","lorvar")]
        zz1 <- zz[zz$het==1,]
        # output
        out[seg,1] <- zz$chrom[1]
        out[seg,2] <- seg
        out[seg,3] <- nrow(zz)
        out[seg,4] <- nrow(zz1)
        out[seg,5] <- median(zz$cnlr)
        if (out[seg,4] > 0) {
            # weighted average of squared log-odds ratio minus variance
            out[seg,6] <- sum(((zz1$valor)^2 - zz1$lorvar)/zz1$lorvar)/sum(1/zz1$lorvar)
        }
    }
    out
}
