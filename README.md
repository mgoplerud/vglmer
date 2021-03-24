# vglmer: Variational Generalized Linear Mixed Effects Regression [![Travis build status](https://travis-ci.com/mgoplerud/vglmer.svg?token=xHM2cTJdHAzcsxnP4SwG&branch=master)](https://travis-ci.com/mgoplerud/vglmer) [![codecov](https://codecov.io/gh/mgoplerud/vglmer/branch/master/graph/badge.svg?token=L8C4260BUW)](https://codecov.io/gh/mgoplerud/vglmer)   [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)


A package to estimate non-linear hierarchical models using the variational algorithms described in [Goplerud (2020)](https://arxiv.org/abs/2007.12300). It also provides the option to improve an initial approximation using marginally augmented variational Bayes (MAVB) also described in the same paper. It can be installed using `devtools`
```
library(devtools); devtools::install_github("mgoplerud/vglmer", dependencies = TRUE)
```

At present, it can fit logistic and negative binomial outcomes with an arbitrary number of random effects. Details on negative binomial inference can be found [here](https://j.mp/goplerud_MAVB_extra) and are more experimential at the moment.

This model accepts "standard" glmer syntax of the form:

```
vglmer(formula = y ~ x + (x | g), data = data, family = 'binomial')

vglmer(formula = y ~ x + (x | g), data = data, family = 'negbin')
```

Many standard methods from `lme4` work, e.g. `fixef`, `coef`, `vcov`, `ranef`, `predict`. Use `format_vglmer` to parse all parameters into a single data.frame. Estimation can be controlled via the numerous arguments to `control` using `vglmer_control`. At the moment, Schemes I, II, and III in Goplerud (2020) correspond to `strong`, `partial`, and `weak`. Unless the model is extremely large, the default (`weak`) should be used.

The package is still "experimental" so some of the inputs and outputs may change! Please make an issue with any concerns you have. 

One known issue is that it is not possible to include a random effect without a main effect (e.g. ``y ~ x + (b | g)`` will fail). This will be fixed shortly.
