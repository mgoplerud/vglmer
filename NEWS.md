# vglmer 1.0.7

*** IMPORTANT: `summary(object)` reported the incorrect variational mean of q(Sigma_j) [variance components]. This has been fixed; nothing changes in estimation and `object$cov` has always extracted the correct variational distribution of q(Sigma_j). Thank you to Ruggero Bellio for bringing this to my attention.

** Added convergence criterion using "relative" ELBO, i.e. stopping when ELBO/N increases by a small amount where N is the number of observations.

** Adjustments to "print", e.g. only showing change in q(alpha,beta) at convergence

* Small adjustments to sparse matrix algebra to improve speed in problems with high-dimensional REs. Thank you to Ruggero Bellio and Hanna Niwinska for finding datasets that illustrate this problem.

* Updates to README to reference new papers and experimental package branches.

* Address certain warnings raised by `Matrix`-package updates.

# vglmer 1.0.6

* Removes unnecessary model preparation steps for `parameter_expansion="translation"` and `factorization_method="strong"`. Improves speed on default settings for models with many random effects.

* Updated references in documentation.

* Small adjustment to tests to prevent failure from update to `waldo` package

# vglmer 1.0.5

** Add gKRLS as an option for smoothing multiple (continuous) covariates. Chang and Goplerud (2024; https://doi.org/10.1017/pan.2023.27) provides more details.

# vglmer 1.0.4

* Adjust `predict.vglmer` to allow for faster predictions on large datasets by not copying and filling in a large sparse matrix. Thank you to Michael Auslen for pointing out this issue.

* Add the option for `terms` to `predict` to allow for predictions for each random effect separately

* Address a bug where predictions with `NA` in new data.frame would fail for certain splines or for cases where `newdata` had a single row.

# vglmer 1.0.3

* Adjust `vglmer` to not throw deprecation messages with Matrix 1.5. Thank you to Mikael Jagan for suggestions on how to adapt the code.

# vglmer 1.0.2

* IMPORTANT: Fixes bug where prediction with only one spline  (and no random effects) was wrong; the non-linear part of the spline was ignored.
* Smaller bug fixes around splines (e.g., for using a single knot) have been added as well as updated tests.

# vglmer 1.0.1

* Patch to address compiler issues on CRAN
* Add links to GitHub to description

# vglmer 1.0.0

* Initial submission to CRAN. Estimates linear, binomial, and negative binomial (experimental) models.
