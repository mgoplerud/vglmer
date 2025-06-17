# vglmer: Variational Generalized Linear Mixed Effects Regression  

Note this version (Poisson branch) is experimental and subject to change. It is used for the results in Goplerud and Auslen (2025); for the results in Goplerud (2022, 2024), please install the master branch; this branch will be eventually be integrated into the main/CRAN version.

[![CRAN status](https://www.r-pkg.org/badges/version/vglmer)](https://CRAN.R-project.org/package=vglmer) [![R-CMD-check](https://github.com/mgoplerud/vglmer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mgoplerud/vglmer/actions/workflows/R-CMD-check.yaml) [![codecov](https://codecov.io/gh/mgoplerud/vglmer/branch/master/graph/badge.svg?token=L8C4260BUW)](https://app.codecov.io/gh/mgoplerud/vglmer)

A package to estimate non-linear hierarchical models using the variational algorithms described in [Goplerud (2022)](https://doi.org/10.1214/21-BA1266) and in [Goplerud (2024)](https://doi.org/10.1017/S0003055423000035). It also provides the option to improve an initial approximation using marginally augmented variational Bayes (MAVB) also described in [Goplerud (2022)](https://doi.org/10.1214/21-BA1266). It can be installed from CRAN or the most-to-update version can be installed using `devtools`.  

```
# CRAN
install.packages("vglmer")
# Up-to-Date GitHub Version
library(devtools)
devtools::install_github("mgoplerud/vglmer", dependencies = TRUE)
```

At present, `vglmer` can fit logistic, linear, Poisson, and negative binomial outcomes with an arbitrary number of random effects. Details on negative binomial inference can be found [here](https://github.com/mgoplerud/vglmer/blob/master/.github/model_addendum.pdf) and are more experimental at the moment.

This package accepts "standard" glmer syntax of the form:

```
vglmer(formula = y ~ x + (x | g), data = data, family = 'binomial')
```

Splines can be estimated using `v_s(x)`, similar to the functionality in `mgcv`, although with many fewer options.

```
vglmer(formula = y ~ v_s(x) + (x | g), data = data, family = 'binomial')
```

Fixed effects can be estimated using the following syntax; this syntax and their estimation is under development.

```
vglmer(formula = y ~ v_fe(id) + x + (x | g), data = data, family = 'binomial')
```

