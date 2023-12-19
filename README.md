# vglmer: Variational Generalized Linear Mixed Effects Regression   
[![CRAN status](https://www.r-pkg.org/badges/version/vglmer)](https://CRAN.R-project.org/package=vglmer) [![R-CMD-check](https://github.com/mgoplerud/vglmer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mgoplerud/vglmer/actions/workflows/R-CMD-check.yaml) [![codecov](https://codecov.io/gh/mgoplerud/vglmer/branch/master/graph/badge.svg?token=L8C4260BUW)](https://app.codecov.io/gh/mgoplerud/vglmer)

A package to estimate non-linear hierarchical models using the variational algorithms described in [Goplerud (2022)](https://arxiv.org/abs/2007.12300) and in [Goplerud (2023)](https://doi.org/10.1017/S0003055423000035). It also provides the option to improve an initial approximation using marginally augmented variational Bayes (MAVB) also described in [Goplerud (2022)](https://arxiv.org/abs/2007.12300). It can be installed from CRAN or the most-to-update version can be installed using `devtools`.  

```
# CRAN
install.packages("vglmer")
# Up-to-Date GitHub Version
library(devtools)
devtools::install_github("mgoplerud/vglmer", dependencies = TRUE)
```

If you are interested in using partially factorized variational inference (Goplerud, Papaspiliopoulos, and Zanella 2023), please switch to the `collapsed` branch and install that version of the package. There are some important differences with this main branch, especially in terms of some `vglmer_control` naming conventions.

At present, `vglmer` can fit logistic, linear, and negative binomial outcomes with an arbitrary number of random effects. Details on negative binomial inference can be found [here](https://github.com/mgoplerud/vglmer/blob/master/.github/model_addendum.pdf) and are more experimental at the moment.

This package accepts "standard" glmer syntax of the form:

```
vglmer(formula = y ~ x + (x | g), data = data, family = 'binomial')
```

Splines can be estimated using `v_s(x)`, similar to the functionality in `mgcv`, although with many fewer options.

```
vglmer(formula = y ~ v_s(x) + (x | g), data = data, family = 'binomial')
```

Many standard methods from `lme4` work, e.g. `fixef`, `coef`, `vcov`, `ranef`, `predict`. Use `format_vglmer` to parse all parameters into a single data.frame. Estimation can be controlled via the numerous arguments to `control` using `vglmer_control`. At the moment, Schemes I, II, and III in Goplerud (2022) correspond to `strong`, `partial`, and `weak`. The default is `strong` which correspond to the strongest (worst) approximation. If the variance of the parameters is of interest, then `weak` will return better results.

Please make an issue on GitHub with any concerns you have.
