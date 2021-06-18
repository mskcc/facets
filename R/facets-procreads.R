# heterozygous and keep flags of the SNPs
procSnps <- function(rcmat, ndepth=35, het.thresh=0.25, snp.nbhd=250, nX=23, unmatched=FALSE, ndepthmax=1000) {
    # keep only chromsomes 1-22 & X for humans and 1-19, X for mice
    # for other genomes (gbuild = udef) nX is number of autosomes plus 1
    rcmat <- rcmat %>%
        tibble::add_column(vafT = 1 - rcmat$TUM.RD/rcmat$TUM.DP) %>%
        tibble::add_column(vafN = 1 - rcmat$NOR.RD/rcmat$NOR.DP) %>%
        dplyr::filter(Chromosome %in% c(1:(nX-1), "X")) %>%
        dplyr::filter((NOR.DP >= ndepth) & (NOR.DP < ndepthmax)) %>%
        dplyr::filter(TUM.DP>0) %>%
        dplyr::rename(chrom = Chromosome) %>%
        dplyr::rename(maploc = Position) %>%
        dplyr::rename(rCountT = TUM.DP) %>%
        dplyr::rename(rCountN = NOR.DP) %>%

        dplyr::arrange(chrom)
    

    # call a snp heterozygous if min(vafN, 1-mafN) > het.thresh
    if (unmatched) {
        if (het.thresh == 0.25) het.thresh <- 0.1
        rcmat$het <- 1*(pmin(rcmat$vafT, 1-rcmat$vafT) > het.thresh & rcmat$rCountT >= 50)
    } else {
        rcmat$het <- 1*(pmin(rcmat$vafN, 1-rcmat$vafN) > het.thresh)
    }
    # scan maploc for snps that are close to one another (within snp.nbhd bases)
    # heep all the hets (should change if too close) and only one from a nbhd
    rcmat$keep <- scanSnp(rcmat$maploc, rcmat$het, snp.nbhd)
    rcmat <- rcmat %>%
        dplyr::filter(keep == 1)
    as.data.frame(rcmat)
}

scanSnp <- function(maploc, het, nbhd) {
    n <- length(maploc)
    zzz <- .Fortran("scansnp",
                    as.integer(n),
                    as.double(maploc),
                    as.double(het),
                    keep=double(n),
                    as.double(nbhd))
    zzz$keep
}

# obtain logR and logOR from read counts and GC-correct logR
counts2logROR <- function(mat, gbuild, unmatched=FALSE, ugcpct=NULL, f=0.2, nX=23) {
    # gc percentage
    mat$gcpct <- rep(NA_real_, nrow(mat))
    # get GC percentages from pctGCdata package
    # loop thru chromosomes

    nchr <- nX - 1 # IMPACT doesn't have X so only 22
    for (i in 1:nchr) {
        ii <- which(mat$chrom==i)
        # allow for chromosomes with no SNPs i.e. not targeted
        if (length(ii) > 0) {
            if (gbuild == "udef") {
                mat$gcpct[ii] <- getGCpct(i, mat$maploc[ii], gbuild, ugcpct)
            } else {
                mat$gcpct[ii] <- getGCpct(i, mat$maploc[ii], gbuild)
            }
        }
    }
    ##### log-ratio with gc correction and maf log-odds ratio steps
    chrom <- mat$chrom
    maploc <- mat$maploc
    rCountN <- mat$rCountN
    rCountT <- mat$rCountT
    vafT <- mat$vafT
    vafN <- mat$vafN
    het <- mat$het
    gcpct <- mat$gcpct
    # compute gc bias
    ncount <- tapply(rCountN, gcpct, sum)
    tcount <- tapply(rCountT, gcpct, sum)
    pctgc <- as.numeric(names(ncount))
    tscl <- sum(ncount)/sum(tcount)
    gcb <- lowess(pctgc, log2(tcount*tscl)-log2(ncount), f=f)
    jj <- match(gcpct, gcb$x)
    gcbias <- gcb$y[jj]
    # compute cn log-ratio (gc corrected) and baf log odds-ratio
    cnlr <- log2(1+rCountT*tscl) - log2(1+rCountN) - gcbias
    # minor allele log-odds ratio and weights
    lorvar <- valor <- rep(NA_real_, length(maploc))
    if (unmatched) {
        # read count matrix for odds ratio etc
        rcmat <- round(cbind(vafT[het==1]*rCountT[het==1], (1-vafT[het==1])*rCountT[het==1]))
        # folded log of Tukey (with 1/6 correction)
        valor[het==1] <- log(rcmat[,1]+1/6) - log(rcmat[,2]+1/6)
        # variance - approximation using delta method
        lorvar[het==1] <- 1/(rcmat[,1]+1/6) + 1/(rcmat[,2]+1/6)
    } else {
        # read count matrix for odds ratio etc
        rcmat <- round(cbind(vafT[het==1]*rCountT[het==1], (1-vafT[het==1])*rCountT[het==1], vafN[het==1]*rCountN[het==1], (1-vafN[het==1])*rCountN[het==1]))
        # log-odds-ratio (Haldane correction)
        valor[het==1] <- log(rcmat[,1]+0.5) - log(rcmat[,2]+0.5) - log(rcmat[,3]+0.5) + log(rcmat[,4]+0.5)
        # variance of log-odds-ratio (Haldane; Gart & Zweifel Biometrika 1967)
        lorvar[het==1] <- (1/(rcmat[,1]+0.5) + 1/(rcmat[,2]+0.5) + 1/(rcmat[,3]+0.5) + 1/(rcmat[,4]+0.5))
    }
    # put them together
    mat$lorvar <- mat$valor <- mat$cnlr <- mat$gcbias <- rep(NA_real_, nrow(mat))
    mat$gcbias <- gcbias
    mat$cnlr <- cnlr
    mat$valor <- valor
    mat$lorvar <- lorvar
    mat
}
