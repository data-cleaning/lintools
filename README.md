[![CRAN](http://www.r-pkg.org/badges/version/lintools)](http://cran.r-project.org/web/packages/lintools/NEWS)
[![Downloads](http://cranlogs.r-pkg.org/badges/lintools)](http://cran.r-project.org/package=lintools/) 
[![status](https://tinyverse.netlify.app/badge/lintools)](https://CRAN.R-project.org/package=lintools)


## lintools

Tools for manipulating systems of linear (in)equations.


The package offers fairly generic functionality for manipulating linear systems.
Some if not all of this is functionality is probably available in R, scattered
accross packages. This package (re)implements such manipulations and offers them
with a basic but consistent interface.

To test the latest beta version, please install it with the
instructions below. We're happy to receive feedback on the [issues
page](https://github.com/data-cleaning/lintools/issues).


#### Compactifying

Simplify sets of equations by removing spurious rows and columns and combining inequations of the form
```
a.x >= 0
a.x <= 0
```
into a single equality
```
a.x == 0
```


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

#### Test for total unimodularity

Given a system of equalities `A.x=b`, there exist integer solutions `x` iff `A` and `b` are integer and `A` is
[totally unimodular](https://en.wikipedia.org/wiki/Unimodular_matrix).


### Installation

To install the latest CRAN version, open an R session and type
```
install.packages("lintools")
```

To install from source:

```bash
git clone https://github.com/markvanderloo/lintools
cd lintools
make install
```
On the OS that shall not be named you need to install
[rtools](https://cran.r-project.org/bin/windows/Rtools/) to build a package.



