[![Build Status](https://travis-ci.org/data-cleaning/lintools.svg?branch=master)](https://travis-ci.org/data-cleaning/lintools)
[![Coverage Status](https://coveralls.io/repos/data-cleaning/lintools/badge.svg)](https://coveralls.io/r/data-cleaning/lintools) 
[![CRAN](http://www.r-pkg.org/badges/version/lintools)](http://cran.r-project.org/web/packages/lintools/NEWS)
[![Downloads](http://cranlogs.r-pkg.org/badges/lintools)](http://cran.r-project.org/package=lintools/) 


## lintools

Tools for manipulating systems of linear (in)equations.

This is a low-level package, to be used by several packages in the
`validate`-suite of packages. It is not finished yet, but it will at
least offer the following functionality.

#### Reduced row echelon

Bring linear systems in [reduced row echelon](https://en.wikipedia.org/wiki/Row_echelon_form) form.

#### Substitute variables

Simplify systems of (in)equations when one or more of the values is known.

#### Eliminate variables

Rewrite systems by eliminating variables. Gaussian elimination for
equalities or Fourier-Motzkin elimination for inequalities.

#### Successive Projection Algorithm

The SPA is used to minimally adjust a vector, in the (weighted) Euclidean
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

