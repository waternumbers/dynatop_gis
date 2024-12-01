# dynatopGIS

<!-- badges: start -->
[![R-CMD-check](https://github.com/waternumbers/dynatopGIS/workflows/R-CMD-check/badge.svg)](https://github.com/waternumbers/dynatopGIS/actions)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/waternumbers/dynatopGIS?branch=master&svg=true)](https://ci.appveyor.com/project/waternumbers/dynatopGIS)
[![codecov](https://img.shields.io/codecov/c/gh/waternumbers/dynatopGIS/master)](https://app.codecov.io/gh/waternumbers/dynatopGIS)
[![dynatopGIS status badge](https://waternumbers.r-universe.dev/badges/dynatopGIS)](https://waternumbers.r-universe.dev)
<!-- badges: end -->

This R package contains the code to generate a dynamic TOPMODEL from GIS
data. A related package for computing dynamic TOPMODEL simulations can be found
[here](https://waternumbers.github.io/dynatop/).

These two packages are intended as successors to the original Dynamic TOPMODEL implementation in R;
the [dynatopmodel](https://CRAN.R-project.org/package=dynatopmodel) package. Currently
(as a design choice) not all the features of the original package are
implemented. 

If there is a feature or improvement you want or you experience problems
please raise an [issue](https://github.com/waternumbers/dynatopGIS/issues)
(or better still submit a pull request with alterations or fixes!) rather than contact the
maintainers directly in the first instance.

Thanks goes to:
* Peter Metcalfe who developed the original port of Dynamic TOPMODEL to R
during his PhD sponsored by the [JBA Trust](https://www.jbatrust.org).
* The [NERC Q-NFM project (NE/R004722/1)](https://www.lancaster.ac.uk/lec/sites/qnfm/) at Lancaster University for supporting this update to dynatop.
* The original developers of Dynamic TOPMODEL, [Keith
Beven](https://www.lancaster.ac.uk/lec/about-us/people/keith-beven) and [Jim
Freer](http://www.bristol.ac.uk/geography/people/jim-e-freer/index.html), who
contributed to the original R port of Dynamic TOPMODEL.

## Using the code

Alongside the packages available on [CRAN](https://cran.r-project.org/) the
latest development version can be installed from
the git repository from within R using the devtools package:

```
devtools::install_github("waternumbers/dynatopGIS")
```

Prebuild packages of the development version (which may be ahead of CRAN) are available from the [r-universe](https://waternumbers.r-universe.dev/#package:dynatopGIS)
