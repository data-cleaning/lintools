[![Build Status](https://travis-ci.org/markvanderloo/supral.svg?branch=master)](https://travis-ci.org/markvanderloo/supral)
[![Coverage Status](https://coveralls.io/repos/markvanderloo/supral/badge.svg)](https://coveralls.io/r/markvanderloo/supral) 
[![CRAN](http://www.r-pkg.org/badges/version/supral)](http://cran.r-project.org/web/packages/supral/NEWS)
[![Downloads](http://cranlogs.r-pkg.org/badges/supral)](http://cran.r-project.org/package=supral/) 


## supral

A C/R implementation of the Successive Projection Algorithm

The SPA algorithm is used to minimally adjust a vector, in the (weighted) Euclidean
sense, to satisfy a set of linear (in)equality constraints.

This package is under development and not on CRAN yet. You can install the
latest beta version from my
[drat](https://cran.rstudio.com/web/packages/drat/index.html) repo as follows
(first install `drat` if you don't already have it).
```
drat::addRepo("markvanderloo")
install.packages("supral")
```

To get started
```
library(supral)
vignette("intro",package="supral")
```

