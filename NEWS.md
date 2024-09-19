# ppgmmga 1.3 (2023-11)

- Package maintainer change.
- Removed package 'ggthemes' from being a dependency.
- Added closed-form formula for computing the Negentropy of GMMs. This 
  can be used in ppgmmga() function call by specifying `approx = "none"`.
  This is at experimental stage.

# ppgmmga 1.2 (2019-07)

- Updated references to JCGS paper.
- Added examples in ppgmmaga.Rd and vignette of possible usage of 
  `spinplot()` from 'msir' package.
- Updated vignette.

# ppgmmga 1.1 (2019-01)

- Added website via pkgdown.
- Fix a warning in `plot.ppgmmga()` when d > 2 and `dim` is set to two 
  coordinates.
- Updated vignette.
  
# ppgmmga 1.0.1 (2018-10)

- Fix a C++ issue. 

# ppgmmga 1.0.0 (2018-10)

- Initial release on CRAN.
