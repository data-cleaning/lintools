[![Build Status](https://travis-ci.org/markvanderloo/lintools.svg?branch=master)](https://travis-ci.org/markvanderloo/lintools)
[![Coverage Status](https://coveralls.io/repos/markvanderloo/lintools/badge.svg)](https://coveralls.io/r/markvanderloo/lintools) 
[![CRAN](http://www.r-pkg.org/badges/version/lintools)](http://cran.r-project.org/web/packages/lintools/NEWS)
[![Downloads](http://cranlogs.r-pkg.org/badges/lintools)](http://cran.r-project.org/package=lintools/) 


## lintools

Tools for manipulating systems of linear (in)equations.

This is a low-level package, to be used by several packages in the
`validate`-suite of packages.

### Successive Projection Algorithm

A C/R implementation of the Successive Projection Algorithm projects a vector
not satisfying a set of linear (in)equality restrictions onto the convex
polytope described by the restrictions.

The SPA algorithm is used to minimally adjust a vector, in the (weighted) Euclidean
sense, to satisfy a set of linear (in)equality constraints.

This package is under development and not on CRAN yet. You can install the
latest beta version from my
[drat](https://cran.rstudio.com/web/packages/drat/index.html) repo as follows
(first install `drat` if you don't already have it).
```
drat::addRepo("data-cleaning")
install.packages("lintools")
```

To get started
```
library(lintools)
vignette("intro",package="lintools")
```

