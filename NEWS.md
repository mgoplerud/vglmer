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
