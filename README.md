# vglmer: Variational Generalized Linear Mixed Effects Regression [![Travis build status](https://travis-ci.com/mgoplerud/vglmer.svg?token=xHM2cTJdHAzcsxnP4SwG&branch=master)](https://travis-ci.com/mgoplerud/vglmer) [![codecov](https://codecov.io/gh/mgoplerud/vglmer/branch/master/graph/badge.svg?token=L8C4260BUW)](https://codecov.io/gh/mgoplerud/vglmer)

A package to estimate non-linear hierarchical models using the variational algorithms described in [Goplerud (2020)](https://j.mp/goplerud_MAVB). It also provides the option to improve an initial approximation using marginally augmented variational Bayes (MAVB) also described in the same paper. 

At present, it can be fit on logistic and negative binomial outcomes with an arbitrary number of random effects. Details on negative binomial inference can be found [here] and are more experimential at the moment.

Syntax
