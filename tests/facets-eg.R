######################################################################
# Code to try out code coverage
######################################################################

library("facets")
# check if .Random.seed exists
seedexists <- exists(".Random.seed")
# save seed
if(seedexists) oldSeed <- .Random.seed
# Alway use the same random seed
set.seed(0xfade)

# running the stomach example in the vignette 
datafile = system.file("extdata", "stomach.csv.gz", package="facets")
# read the data
rcmat = readSnpMatrix(datafile)
# fit segmentation tree
xx = preProcSample(rcmat)
# estimate allele specific copy numbers
oo=procSample(xx,cval=150)
# EM fit version 1
fit=emcncf(oo)
# EM fit version 2 (removed in v0.6.0)
# fit2=emcncf2(oo)
# finished

# Reset to previous random seed
if(seedexists) .Random.seed <- oldSeed
