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
    jointseg <- jointseg[is.finite(jointseg$segs),]
    # initialize output table
    out <- as.data.frame(matrix(0, length(table(paste(jointseg$chrom, jointseg$segs))), 6))
    names(out) <- c("chr","seg","num.mark","nhet","cnlr.median","mafR")
    l <- 0
    for(chr in unique(jointseg$chrom)) {
        zz <- jointseg[jointseg$chrom==chr, c("maploc","rCountT","vafT","rCountN","vafN","segs","cnlr","het","valor","lorvar")]
        zz1 <- zz[zz$het==1,]
        nseg <- max(zz$segs)
        for(seg in 1:nseg) {
            l <- l+1
            # segment indicator
            iseg <- zz$segs==seg
            iseg1 <- zz1$segs==seg
            # output
            out[l,1] <- chr
            out[l,2] <- seg
            out[l,3] <- sum(1*iseg)
            out[l,4] <- sum(zz$het[iseg])
            out[l,5] <- median(zz$cnlr[iseg])
            if (out[l,4] > 0) {
                # weighted average of squared log-odds ratio minus variance
                out[l,6] <- sum(((zz1$valor[iseg1])^2 - zz1$lorvar[iseg1])/zz1$lorvar[iseg1])/sum(1/zz1$lorvar[iseg1])
            }
        }
    }
    out
}
