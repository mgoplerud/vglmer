## Resubmission

This is an update on vglmer 1.0.5 that adds some additional, minor,
functionalities. It considerably improves the speed of
`parameter_expansion="translation"` for `factorization_method="strong"` in large
models by speeding up the model preparation.

It also addresses a failure in the tests from updates to `testthat` and `waldo`.

## R CMD check results

There were no ERRORs, WARNINGs or NOTES. 

## Downstream dependencies

I have checked the downstream dependencies using
https://github.com/r-devel/recheck and all pass R CMD CHECK.