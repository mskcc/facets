# facets
Algorithm to implement Fraction and Allele specific Copy number Estimate from Tumor/normal Sequencing.

**Travis CI** (_Linux & macOS_) [![Build Status](https://travis-ci.com/mskcc/facets.svg?branch=master)](https://travis-ci.com/mskcc/facets)  **Test coverage** [![codecov.io](https://codecov.io/github/mskcc/facets/coverage.svg?branch=master)](https://codecov.io/github/mskcc/facets?branch=master)

You can install the current version (along with the vignette) using the command

```R
devtools::install_github("mskcc/facets", build_vignettes = TRUE)
```

pctGCdata is a required package. So install that also (needs to be done only once)

```R
devtools::install_github("mskcc/pctGCdata")
```
If you get an error message about pctGCdata use

```R
devtools::install_github("veseshan/pctGCdata")
```

## NOTES

### 2016_11_11 (version 0.5.6)

The new version estimates the log-ratio level corresponding to the diploid state. It is embedded into the procSample call.
In terms of using the package you can now do:

```R
rcmat <- readSnpMatrix(filename, ...)
xx <- preProcSample(rcmat, ...)
# specify cval you like
oo <- procSample(xx, cval = 300)
```

And go straight to
```R
emcncf(oo)
```
The ```emcncf2(oo)``` option imposes a clonal cluster structure. This function is currently being reworked. Please use with caution.

The output of procSample now has 4 elements:

* ```jointseg``` – same as before
* ```out``` – as before with 3 additional columns: cf, tcn, lcn
* ```dipLogR``` – the estimated location of diploid log-ratio value
* ```flags``` – this gives an indication of whether the dipLogR is estimated well.

If ```flags``` is NULL then no obvious problem with the dipLogR estimate. It can have two other comments: "mafR not sufficiently small" and "could be polyclonal 1 copy loss". The first one means that there aren't segments with sufficiently balanced alleles and so the estimate may not be great. The second one means that it looks like genome doubling; the only lower level segments are 2+1 and 2+0 (from 2+2); a model without genome doubling but single copy loss (1+0 from 1+1) with two different cellular fraction could fit better.

Note: I am not claiming that I have covered all possible scenarios that lead to bad dipLogR estimate. If you come across anything that doesn't seem right but flags is NULL let me know.

