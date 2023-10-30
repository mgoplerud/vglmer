# vglmer: Variational Generalized Linear Mixed Effects Regression   
[![CRAN status](https://www.r-pkg.org/badges/version/vglmer)](https://CRAN.R-project.org/package=vglmer) [![R-CMD-check](https://github.com/mgoplerud/vglmer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mgoplerud/vglmer/actions/workflows/R-CMD-check.yaml) [![codecov](https://codecov.io/gh/mgoplerud/vglmer/branch/master/graph/badge.svg?token=L8C4260BUW)](https://app.codecov.io/gh/mgoplerud/vglmer)

*Note this version (collapsed branch) is preliminary and used for the results in
Goplerud, Papaspiliopoulos, and Zanella (2023); for the results in Goplerud
(2022, 2023), please install the master branch; this branch will be eventually
be integrated into the main/CRAN version*

A package to estimate non-linear hierarchical models using the variational algorithms described in [Goplerud (2022)](https://doi.org/10.1214/21-BA1266) and subsequent developments ([Goplerud 2023](https://doi.org/10.1017/S0003055423000035), Goplerud, Papaspiliopoulos, and Zanella 2023). It can be installed from CRAN or the most-to-update version can be installed using `devtools`.

```
# CRAN
install.packages("vglmer")
# Up-to-Date GitHub Version
library(devtools)
devtools::install_github("mgoplerud/vglmer", dependencies = TRUE)
```

At present, it can fit logistic, linear, and negative binomial outcomes with an arbitrary number of random effects. Details on negative binomial inference can be found [here](https://github.com/mgoplerud/vglmer/blob/master/.github/model_addendum.pdf) and are more experimental at the moment.

This package accepts "standard" glmer syntax of the form:

```
vglmer(formula = y ~ x + (x | g), data = data, family = 'binomial')
```

Splines can be estimated using `v_s(x)`, similar to the functionality in `mgcv`, although with many fewer options.

```
vglmer(formula = y ~ v_s(x) + (x | g), data = data, family = 'binomial')
```

Many standard methods from `lme4` work, e.g. `fixef`, `coef`, `vcov`, `ranef`, `predict`. Use `format_vglmer` to parse all parameters into a single data.frame. Estimation can be controlled via the numerous arguments to `control` using `vglmer_control`. At the moment, Schemes I, II, and III in Goplerud (2022) correspond to `strong`, `intermediate`, and `weak`. The "partially factorized" scheme explored in Goplerud, Papaspiliopoulos, and Zanella (2023) can be used with `partially_factorized`, although some features such as splines are currently not enabled for this factorization scheme. 

The default is `strong` which correspond to the strongest (worst) approximation. If the variance of the parameters is of interest, then `partially_factorized` or `weak` will return better results.

Please make an issue on GitHub with any concerns you have.
