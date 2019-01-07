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
devtools::install_github("luca-scr/ppgmmga", build_vignettes = TRUE)
```

## Usage

The methodology implemented in the package is describe in the paper referenced below.

Usage of the main functions and some examples are included in the vignette **A quick tour of ppgmmga**, which is available as

```{r}
vignette("ppgmmga")
```


## References

Scrucca, L. and Serafini, A. (2018) Projection pursuit based on Gaussian mixtures and evolutionary algorithms. To appear in *Journal of Computational and Graphical Statistics*.