# 12/08/2017 Function modeled after the one contributed by Dario Beraldi
# (https://github.com/dariober) it uses data.table package
readSnpMatrixDT <- function(pileup, err.thresh= Inf, del.thresh= Inf) {
    rcmat<- fread(sprintf('gunzip -c %s', pileup),
                  select= c('Chromosome', 'Position', 'File1R', 'File1A',
                            'File1E', 'File1D', 'File2R', 'File2A', 'File2E',
                            'File2D'))
    setnames(rcmat, c('File1R', 'File1A', 'File2R', 'File2A'),
                    c('NOR.RD', 'NOR.DP', 'TUM.RD', 'TUM.DP'))
    rcmat<- rcmat[File1E <= err.thresh & File2E <= err.thresh &
                  File1D <= del.thresh & File2D <= del.thresh,
                  list(Chromosome, Position, NOR.DP, NOR.RD, TUM.DP, TUM.RD)]
    rcmat[, NOR.DP := NOR.DP + NOR.RD]
    rcmat[, TUM.DP := TUM.DP + TUM.RD]
    rcmat[, Chromosome := sub("chr", "", Chromosome, fixed= TRUE)]
    setcolorder(rcmat, c('Chromosome', 'Position', 'NOR.DP', 'NOR.RD', 'TUM.DP', 'TUM.RD'))
    return(rcmat)
}
