# GeneCastR

Introduction time and spread of genetic variants under Wright–Fisher dynamics.

Fits a Wright–Fisher model with genic selection and recurrent gene flow to
panels of allele-frequency surveys, with optional per-deme hierarchies over
the growth advantage and the effective size, and reports projected
frequencies, presence and exceedance probabilities, and posteriors over the
introduction time and the founding time of the lineage that was detected.

## Installation

```r
# install.packages("remotes")
remotes::install_github("mrc-ide/genecastR", build_vignettes = TRUE)
```

`build_vignettes = TRUE` costs a few minutes and is worth it: without it the
tutorials listed at the end are not installed and `vignette()` will not find
them. Compilation needs a C++ toolchain, as below.

## Getting started

Needs R >= 4.1 and a C++ compiler. No Fortran toolchain is required --
the package links no BLAS or LAPACK, deliberately (see `src/Makevars`).

```r
install.packages(c("R6", "Rcpp", "roxygen2", "testthat", "pkgbuild",
                   "pkgload", "posterior", "ggplot2", "scales"))

Rcpp::compileAttributes(".")                      # regenerates RcppExports
roxygen2::roxygenise(".", roclets = c("rd", "namespace"))   # then reads them
pkgload::load_all(".")
testthat::test_local(".")
```

Regenerating the stored fit (`data-raw/uganda_fit.R`) writes into the source
tree, while `gc_example_fit()` reads from the installed library -- reinstall
after running it.

Run those two in that order after changing any `// [[Rcpp::export]]`
signature: `compileAttributes()` rewrites `R/RcppExports.R`, and roxygen
reads the result. Installing without them leaves the `.Rd` describing the old
argument list, which `R CMD check` reports as a code/documentation mismatch.

**Do not benchmark a `load_all()` build.** `pkgbuild::compile_dll()`, which
`load_all()` calls, appends `-O0` after R's own `-O2`, so the compiled code
runs roughly six times slower than the installed package. Time things against
`R CMD INSTALL` instead:

```sh
R CMD INSTALL --preclean .
```

`man/` is empty in the repository because roxygen generates it; `document()`
is the first thing to run.

Check the compiled toolchain:

```r
gc_toolchain()
#> $openmp       TRUE / FALSE
#> $max_threads   n
#> $cxx_standard 202002
```

`openmp = FALSE` is expected on macOS with Apple clang alone, and means the
per-deme loop runs serially. That is correct but slower; phase 7 needs libomp
from <https://mac.r-project.org/openmp/> first.

## Tutorials

Start here for a quickstart:

```r
vignette("GeneCastRQuickstart", package = "GeneCastR")
```

Then, for the the full model and parameterisation, the hierarchical options
and the full list of settings:

```r
vignette("GeneCastRTutorial", package = "GeneCastR")
```