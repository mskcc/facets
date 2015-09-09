# heterozygous and keep flags of the SNPs
procSnps <- function(filename, ndepth=35, het.thresh=0.25, snp.nbhd=250, chromlevels=c(1:22,"X"), unmatched=FALSE) {
    # read the SNP matrix
    xx <- scan(filename, what=list(Chrom="", Pos=0, Ref="", Alt="", TUM.DP=0, TUM.Ap=0, TUM.Cp=0, TUM.Gp=0, TUM.Tp=0, TUM.An=0, TUM.Cn=0, TUM.Gn=0, TUM.Tn=0, NOR.DP=0, NOR.Ap=0, NOR.Cp=0, NOR.Gp=0, NOR.Tp=0, NOR.An=0, NOR.Cn=0, NOR.Gn=0, NOR.Tn=0), sep="\t", skip=1)
    # remove chr if present in Chrom
    if (xx$Chrom[1] == "chr1") xx$Chrom <- gsub("chr","",xx$Chrom)
    # data frame
    xx <- as.data.frame(xx)
    # keep only chromsomes 1-22, X
    chr.keep <- xx$Chrom %in% chromlevels
    # keep only snps with normal read depth between ndepth and 1000
    depthN.keep <- (xx$NOR.DP >= ndepth) & (xx$NOR.DP < 1000)
    # reduce the data frame to these snps
    xx <- xx[chr.keep & depthN.keep,]
    # ref allele column
    RefID <- match(xx$Ref, c("A","C","G","T"))
    ii <- 1:nrow(xx)
    # output data frame
    out <- list()
    out$chrom <- xx$Chrom
    out$maploc <- xx$Pos
    out$rCountT <- xx$TUM.DP
    out$rCountN <- xx$NOR.DP
    out$vafT <- 1 - (as.matrix(xx[,6:9])[cbind(ii,RefID)] + as.matrix(xx[,10:13])[cbind(ii,RefID)])/out$rCountT
    out$vafN <- 1 - (as.matrix(xx[,15:18])[cbind(ii,RefID)] + as.matrix(xx[,19:22])[cbind(ii,RefID)])/out$rCountN
    # make chromosome ordered and numeric
    out$chrom <- as.numeric(ordered(out$chrom, levels=chromlevels))
    # call a snp heterozygous if min(vafN, 1-mafN) > het.thresh
    if (unmatched) {
        if (het.thresh == 0.25) het.thresh <- 0.1
        out$het <- 1*(pmin(out$vafT, 1-out$vafT) > het.thresh & out$rCountT >= 50)
    } else {
        out$het <- 1*(pmin(out$vafN, 1-out$vafN) > het.thresh)
    }
    # scan maploc for snps that are close to one another (within snp.nbhd bases)
    # heep all the hets (should change if too close) and only one from a nbhd
    out$keep <- scanSnp(out$maploc, out$het, snp.nbhd)
    as.data.frame(out)
}

# procSnps code above redone to use fread in data.table
procSnpsDT <- function(filename, ndepth=35, het.thresh=0.25, snp.nbhd=250, chromlevels=c(1:22,"X"), unmatched=FALSE) {
    # read the SNP matrix
    xx <- fread(paste("gunzip -c", filename), header=T, sep="\t")
    # remove chr if present in Chrom add it chromlevels
    if (xx$Chrom[1] == "chr1") chromlevels <- paste("chr", chromlevels, sep="")
    # keep only chromsomes 1-22, X
    chr.keep <- xx$Chrom %in% chromlevels
    # keep only snps with normal read depth between ndepth and 1000
    depthN.keep <- (xx$NOR.DP >= ndepth) & (xx$NOR.DP < 1000)
    # reduce the data frame to these snps
    xx <- xx[chr.keep & depthN.keep,]
    # ref allele column
    RefID <- match(xx$Ref, c("A","C","G","T"))
    ii <- 1:nrow(xx)
    # output data frame
    out <- list()
    out$chrom <- xx$Chrom
    out$maploc <- xx$Pos
    out$rCountT <- xx$TUM.DP
    out$rCountN <- xx$NOR.DP    
    out$vafT <- 1 - (as.matrix(xx[,.(TUM.Ap,TUM.Cp,TUM.Gp,TUM.Tp)])[cbind(ii,RefID)] + as.matrix(xx[,.(TUM.An,TUM.Cn,TUM.Gn,TUM.Tn)])[cbind(ii,RefID)])/out$rCountT
    out$vafN <- 1 - (as.matrix(xx[,.(NOR.Ap,NOR.Cp,NOR.Gp,NOR.Tp)])[cbind(ii,RefID)] + as.matrix(xx[,.(NOR.An,NOR.Cn,NOR.Gn,NOR.Tn)])[cbind(ii,RefID)])/out$rCountN
    # make chromosome ordered and numeric
    out$chrom <- as.numeric(ordered(out$chrom, levels=chromlevels))
    # call a snp heterozygous if min(vafN, 1-mafN) > het.thresh
    if (unmatched) {
        if (het.thresh == 0.25) het.thresh <- 0.1
        out$het <- 1*(pmin(out$vafT, 1-out$vafT) > het.thresh & out$rCountT >= 50)
    } else {
        out$het <- 1*(pmin(out$vafN, 1-out$vafN) > het.thresh)
    }
    # scan maploc for snps that are close to one another (within snp.nbhd bases)
    # heep all the hets (should change if too close) and only one from a nbhd
    out$keep <- scanSnp(out$maploc, out$het, snp.nbhd)
    as.data.frame(out)
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
counts2logROR <- function(mat, unmatched=FALSE, f=0.2) {
    out <- mat[mat$keep==1,]
    # gc percentage
    out$gcpct <- rep(NA_real_, nrow(out))
    # load gc data 
    if(!exists("gcpctdb")) data(hg19gcpct, package="facets", envir=environment())
    # loop thru chromosomes
    nchr <- max(mat$chrom) # IMPACT doesn't have X so only 22
    for (i in 1:nchr) {
        ii <- out$chrom==i
        jj <- ceiling((out$maploc[ii]-450)/100)
        jj[jj < 1] <- 1 # fix issues (maploc < 500) with chr17
        out$gcpct[ii] <- gcpctdb[[i]][jj]
    }
    ##### log-ratio with gc correction and maf log-odds ratio steps
    chrom <- out$chrom
    maploc <- out$maploc
    rCountN <- out$rCountN
    rCountT <- out$rCountT
    vafT <- out$vafT
    vafN <- out$vafN
    het <- out$het
    gcpct <- out$gcpct
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
        # folded log ot Tukey (with 1/6 correction)
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
    out$lorvar <- out$valor <- out$cnlr <- out$gcbias <- rep(NA_real_, nrow(out))
    out$gcbias <- gcbias
    out$cnlr <- cnlr
    out$valor <- valor
    out$lorvar <- lorvar
    out
}
