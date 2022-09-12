## Resubmission

This is a resubmission. In this version I have:
  
  - Addressed the issues with the two links noted in the previous submission.
  - Changed `quiet=T` to `quiet=TRUE` in the function `vglmer_control`
  - The documentation for five functions has been adjusted to include a `value` tag as noted by the CRAN team: `sl_vglmer`, `posterior_samples.vglmer`, `v_s`, `vglmer_control`, `vglmer-class`. All exported functions in the manual now contain a "Values" section.
  - The commented out code lines in examples for `predict.vglmer` have been removed.
  
## R CMD check results

There were no ERRORs or WARNINGs. 

There was one note. This is an initial submission.

  * checking CRAN incoming feasibility ... [12s] NOTE
  Maintainer: 'Max Goplerud <mgoplerud@pitt.edu>'
  
  New submission

## Downstream dependencies

There are no downstream dependencies.