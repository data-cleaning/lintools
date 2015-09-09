[![Build Status](https://travis-ci.org/data-cleaning/lintools.svg?branch=master)](https://travis-ci.org/data-cleaning/lintools)
[![Coverage Status](https://coveralls.io/repos/data-cleaning/lintools/badge.svg)](https://coveralls.io/r/data-cleaning/lintools) 
[![CRAN](http://www.r-pkg.org/badges/version/lintools)](http://cran.r-project.org/web/packages/lintools/NEWS)
[![Downloads](http://cranlogs.r-pkg.org/badges/lintools)](http://cran.r-project.org/package=lintools/) 


## lintools

Tools for manipulating systems of linear (in)equations.

This is a low-level package, to be used by several packages in the
`validate`-suite of packages. It is not finished yet, but it will at
least offer the following functionality.


#### Blocking

Separate matrices into independent blocks.

#### Reduced row echelon

Bring linear systems in [reduced row echelon](https://en.wikipedia.org/wiki/Row_echelon_form) form.

#### Substitute variables

Simplify systems of (in)equations when one or more of the values is known.

#### Eliminate variables

Rewrite systems by eliminating variables. Gaussian elimination for
equalities or Fourier-Motzkin elimination for inequalities.

#### Feasibility checks

Check whether a system of linear (in)equations has any solution.

#### Project on convex polytope

Given a vector not satisfying a set of (in)equations, project it onto
the convex polytope described by the restrictions.

#### Extract partial solutions

Extract values that are uniquely defined (implied) by a set of (in)equations.



### Installation

This package is under development and not on CRAN yet. You can install the
latest beta version from my
[drat](https://cran.rstudio.com/web/packages/drat/index.html) repo as follows
(first install `drat` if you don't already have it).
```
drat::addRepo("data-cleaning")
install.packages("lintools")
```


