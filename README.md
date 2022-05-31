# vglmer: Variational Generalized Linear Mixed Effects Regression [![R-CMD-check](https://github.com/mgoplerud/vglmer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mgoplerud/vglmer/actions/workflows/R-CMD-check.yaml) [![codecov](https://codecov.io/gh/mgoplerud/vglmer/branch/master/graph/badge.svg?token=L8C4260BUW)](https://codecov.io/gh/mgoplerud/vglmer)

A package to estimate non-linear hierarchical models using the variational algorithms described in [Goplerud (2022a)](https://arxiv.org/abs/2007.12300) and in [Goplerud (2022b)](https://www.dropbox.com/s/wb98qa7dqledw0i/Goplerud_MMM.pdf?dl=0). It also provides the option to improve an initial approximation using marginally augmented variational Bayes (MAVB) also described in Goplerud (2022a). It can be installed using `devtools`.

```
library(devtools)
devtools::install_github("mgoplerud/vglmer", dependencies = TRUE)
```

At present, it can fit logistic, linear, and negative binomial outcomes with an arbitrary number of random effects. Details on negative binomial inference can be found [here](https://j.mp/goplerud_MAVB_extra) and are more experimential at the moment.

This model accepts "standard" glmer syntax of the form:

```
vglmer(formula = y ~ x + (x | g), data = data, family = 'binomial')
```

Splines can be estimated using `v_s(x)`, similar to the functionality in `mgcv`, although with many fewer options.

```
vglmer(formula = y ~ v_s(x) + (x | g), data = data, family = 'binomial')
```

Many standard methods from `lme4` work, e.g. `fixef`, `coef`, `vcov`, `ranef`, `predict`. Use `format_vglmer` to parse all parameters into a single data.frame. Estimation can be controlled via the numerous arguments to `control` using `vglmer_control`. At the moment, Schemes I, II, and III in Goplerud (2022) correspond to `strong`, `partial`, and `weak`. The default is `strong` which correspond to the strongest (worst) approximation. If the variance of the parameters is of interest, then `weak` will return better results.

The package is still "experimental" so some of the inputs and outputs may change! Please make an issue with any concerns you have. One known issue is that it is not possible to include a random effect without a corresponding main effect (e.g. ``y ~ x + (b | g)`` will fail).
