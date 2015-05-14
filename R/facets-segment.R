# recursively split the chromosomes using cval and save the tree structure
fit.cpt.tree <- function(genomdat, edgelim=10, cval=25) {
    # genomdat has 3 columns: logR, logOR and het indicator
    n <- nrow(genomdat)
    seg.end <- c(0,n)
    # this refers to the row in seg.tree
    pnode.vec <- c(0,0)
    # tree structure of splits
    # each row refers to a node and has 4 values -
    # parent node, start and end of segment and maximal statistic
    seg.tree <- NULL
    k <- length(seg.end)
    change.loc <- NULL
    # row number of seg.tree
    current.node <- 0
    while (k > 1) {
        current.n <- seg.end[k]-seg.end[k-1]
        # tree structure
        parent.node <- pnode.vec[k]
        current.node <- current.node + 1
        if (current.n > 1) {
            current.genomdat <- genomdat[(seg.end[k-1]+1):seg.end[k],]
            # heterozygous positions
            het <- current.genomdat[,3]==1
            # number of heterozygous positions
            nhet <- sum(1*het)
            # rank, center and scale (to make it unit variance) the data
            # log-ratio is for all SNPs
            current.genomdat[,1] <- (rank(current.genomdat[,1]) - (current.n+1)/2)/sqrt((current.n+1)*current.n/12)
            # maf is only for germline heterozygous snps; set to 0 in segsnps
            if (nhet > 0) {
                current.genomdat[het,2] <- (rank(current.genomdat[het,2]) - (nhet+1)/2)/sqrt((nhet+1)*max(nhet,1)/12)
            }
            # cumulative heterozygous counts
            current.genomdat[,3] <- cumsum(current.genomdat[,3])
            # call segmentation code
            zzz <- .Fortran("t2maxo",
                            n = as.integer(current.n),
                            sx = as.double(current.genomdat[,1:2]),
                            ihet = as.double(current.genomdat[,3]),
                            iseg = integer(2),
                            ostat = double(1))
            # if ostat > cval, there are 2 changepoints, else 0.
            # 35 seems a reasonable cutoff based on null simulations
            zzz$ncpt <- ifelse(zzz$ostat > cval, 1, 0)
        } else {
            zzz$ncpt <- 0
            zzz$ostat <- 0 # make the statistic 0 for the tree structure
        }
        # add the current node
        seg.tree <- c(seg.tree, parent.node, seg.end[k-(1:0)], zzz$ostat)
        ncpt <- zzz$ncpt
        if (ncpt==1) {
            # segment lengths
            iseg <- diff(c(0, zzz$iseg, current.n))
            # if there are 3 pieces
            if (length(iseg)==3) {
                # if first piece is of length < edgelim merge to the right
                if (iseg[1] < edgelim) iseg[1:2] <- c(0, iseg[1]+iseg[2])
                # if third piece is of length < edgelim merge to the left
                if (iseg[3] < edgelim) iseg[2:3] <- c(iseg[2]+iseg[3], 0)
            }
            # create the change point locations
            iseg <- cumsum(iseg)
            # keep only the interior ones
            iseg <- iseg[iseg>0 & iseg < current.n]
            # if no interior change points after merging ncpt=0
            if(length(iseg)==0) ncpt <- 0
        }        
        if (ncpt==0) {
            change.loc <- c(change.loc,seg.end[k])
            seg.end <- seg.end[-k]
            pnode.vec <- pnode.vec[-k]
        }
        if (ncpt==1) {
            seg.end <- c(seg.end[1:(k-1)],seg.end[k-1]+iseg,seg.end[k])
            pnode.vec <- c(pnode.vec[-k], rep(current.node, length(iseg)+1))
        }
        k <- length(seg.end)
    }
    seg.ends <- unique(c(0,rev(change.loc)))
    list(seg.ends=seg.ends, seg.tree=matrix(seg.tree, ncol=4, byrow=TRUE))
}

# prune an existing seg.tree using a larger clar
prune.cpt.tree <- function(seg.tree, cval=25) {
    k <- nrow(seg.tree)
    keep <- rep(0, k)
    # first row is the whole chromosome; so it is always kept
    keep[1] <- 1
    # now check all daughter nodes and keep them depending on parent node
    for(i in 2:k) {
        parent.node <- seg.tree[i,1]
        # keep a segment if the parent segment is kept and stat > cval
        if (keep[parent.node]==1 & seg.tree[parent.node, 4] > cval) keep[i] <- 1
    }
    # segments that are kept (0 and the end of each kept segment)
    c(0, sort(unique(seg.tree[keep==1,3])))
}

# segment by looping over the chromosomes
segsnps <- function(mat, cval=25) {
    # keep the original data
    mat0 <- mat
    # keep only rows that have finite values for cnlr
    ii <- is.finite(mat$cnlr)
    # keep only the necessary variables for segmentation
    mat <- mat[ii, c("chrom","cnlr","valor","het")]
    # set valor=0 for homozygous snps
    mat$valor[mat$het==0] <- 0
    # take absolute value of valor
    mat$valor <- abs(mat$valor)
    # initialize segment indicator
    mat$segs <- rep(NA_real_, nrow(mat))
    # loop over chromosomes
    nchr <- max(mat$chrom) # IMPACT doesn't have X so only 22
    # initialize segmentation tree
    seg.tree <- list()
    for(i in 1:nchr) {
        genomdat <- as.matrix(mat[mat$chrom==i, c("cnlr","valor","het")])
        # fit segment tree
        tmp <- fit.cpt.tree(genomdat, cval=cval)
        seg.tree[[i]] <- tmp$seg.tree
        # segment indicator
        seg.widths <- diff(tmp$seg.ends)
        mat$segs[mat[,1]==i] <- rep(1:length(seg.widths), seg.widths)
    }
    attr(seg.tree, "cval") <- cval
    # add segs to original matrix
    mat0$segs <- rep(NA_real_, nrow(mat0))
    mat0$segs[ii] <- mat$segs
    # return matrix
    list(seg.tree=seg.tree, jointseg=mat0)
}
