
# dsftools

<!-- badges: start -->
<!-- badges: end -->

The goal of `dsftools` is to provide a library of common functions for fisheries 
analyses within the Alaska Department of Fish & Game, Division of Sport Fish.

The (multivariate) hope is that R-users within the Division can:
* More easily share code we're written, and make use of each other's code
* Standardize our methods for common analyses
* Collaborate and improve each other's code
* Get some Git(hub) & package development learning for free!

## Installation

If you just want to *use* the functions in `dsftools`, you can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ADFG-DSF/dsftools")
```

## What is an R package, and why is `dsftools` structured like one?

An R package is a standardized bundle of R functions, function documentation, 
example data, and most importantly, metadata.  It can be loaded into your R package
library, and you will be able to use the functions, access the function documentation
(i.e. `help()`), and run example scripts within R/Rstudio.

You may notice the files DESCRIPTION and NAMESPACE at the root directory; these
are important files that tell R how to interpret the package and its dependencies.
The /man folder contains all function documentation.



