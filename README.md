
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

## A brief overview of commonly-used functions

#### For ASL (age, sex, and length) summaries:

* `ASL_table()` generates a summary table of ASL from vectors of data, and is robust to multiple sampling scenarios (stratified vs. pooled, total abundance known vs. estimated vs. unknown).
  - `verify_ASL_table()` provides a graphical check of the methods used in `ASL_table()` via simulation.
* `ASL_boilerplate()` automagically creates text to accompany an ASL summary analysis or Operational Plan, to be written to Rmarkdown.  This function may be used with inputs to `ASL_table()`, or directly from inputs describing the sampling scenario.
  
#### Miscellaneous functions

* `logit()` and `expit()`: functions for logit and inverse-logit

## What is an R package, and why is `dsftools` structured like one?

An R package is a standardized bundle of R functions, function documentation, 
example data, and most importantly, metadata.  It can be loaded into your R package
library, and you will be able to use the functions, access the function documentation
(i.e. `help()`), and run example scripts within R/Rstudio.

You may notice the files DESCRIPTION and NAMESPACE at the root directory; these
are important files that tell R how to interpret the package and its dependencies.
The /man folder contains all function documentation.

## I've written some cool R functions.  Can I contribute them?

Yes please!!!  That's the whole point of this package.

If you're new to the R package development workflow or just want to learn more, I *highly* recommend checking out https://r-pkgs.org/ because Hadley and Jennifer explain it all much better than I can.

Definitely check out Adam's excellent writeup on the use of Git & Github: https://adfg-dsf.github.io/Git_book/ as well as https://happygitwithr.com/

### A few ground rules so we're all on the same page:

* **Use `devtools`**: The `devtools` package automates a lot of the package development tasks, and integrates well with RStudio.  Make sure you have the latest version.

* **Use `roxygen2`**: The `roxygen2` package automates building package documentation and metadata.  Notice how the lines preceding each function within this package have tags that look like `#' @description`, etc.  These "roxygen tags" create documentation and/or metadata in the appropriate places.  
  - A couple of important ones: `@export` is necessary to export your function to the final package at the user level (an unexported function will be accessible to other `dsftools` functions, but not the user), and `@importFrom` will import functions from external packages.
  - **An annoying step to enable `Roxygen`**: In RStudio, navigate to Tools > Project Options > Build Tools, then check "Generate R Documentation with Roxygen", then "Configure..." and check all the boxes in the next popup.
  
* **Add documentation**: Add a help file with Roxygen for each function you contribute, and (for now) add a blurb in the README file.

* **Add some automated tests**: In tests/testthat/test_dsftools.R there are series of automated tests.  These will run silently if `dsftools`'s functions produce the desired output, and will complain otherwise.  This is not required, but if you added a function, please add whatever tests you can think of!  This provides a safeguard in case errors are inadvertently introduced, or (more likely) some update to R or a package dependency causes something to break.  Tests aren't for now, tests are for the future.

* **Version up**: If you add or change stuff, bump up the Version in the DESCRIPTION file to let R know that this is a new version.  I'm not sure if this is absolutely necessary, but it's a good practice.

* I don't know yet if the **Fork/Pull Request** or **Collaborate** model makes more sense.  For now, let's try **Fork/Pull Request** for practice.

### A few RStudio buttons you'll use

In the Build pane of RStudio:

* **Install** rebuilds the package and re-generates all Roxygen.
* **Test** runs all automated tests in the /tests folder.
* **Check** runs all tests, plus many additional checks to make sure that the package is functioning correctly and is complete in every conceivable way.  It's totally possible to have a working package that does not completely pass Check, but since this package is collaborative, *do your best to address any problems that Check identifies*.


