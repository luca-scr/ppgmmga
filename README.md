# ppgmmga

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/ppgmmga)](https://cran.r-project.org/package=ppgmmga)
[![CRAN\_MonthlyDownloads](http://cranlogs.r-pkg.org/badges/ppgmmga)](https://cran.r-project.org/package=ppgmmga)


An R package accompanying the paper *Projection pursuit based on Gaussian mixtures and evolutionary algorithms* by Luca Scrucca and Alessio Serafini (2018).

## Installation

You can install the released version of `ppgmmga` from CRAN:

```{r}
install.packages("ppgmmga")
```

or the development version from GitHub:

```{r}
# install.packages("devtools")
devtools::install_github("luca-scr/ppgmmga")
```

## Usage

Usage of the main functions and several examples are included in the
papers shown in the references section below.

For an intro see the vignette **A quick tour of ppgmmga**, which is available
as

```{r}
vignette("ppgmmga")
```

Note that if the package is installed from GitHub the vignette is not
automatically created. However, it can be created when installing from
GitHub with the code:

```{r}
devtools::install_github("luca-scr/ppgmmga", build_vignettes = TRUE)
```

## References

Scrucca, L. and Serafini, A. (2018) Projection pursuit based on Gaussian mixtures and evolutionary algorithms. *Under review*.