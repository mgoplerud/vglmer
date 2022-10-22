context("C++ Verification")

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(111)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that("Calculate E[alpha alpha^T] comparing cpp and base R", {
  loop_outer_alpha <- function(vi_alpha_mean, vi_alpha_decomp, outer_alpha_RE_positions) {
    # Must be such that t(vi_alpha_decomp) %*% vi_alpha_decomp = VAR
    store_oa <- as.list(rep(NA, length(outer_alpha_RE_positions)))
    store_alpha_outer <- as.list(rep(NA, length(outer_alpha_RE_positions)))
    counter_j <- 1
    for (j in outer_alpha_RE_positions) {
      # cat('.')
      summed_oa <- summed_alpha_outer <- array(0, dim = rep(length(j[[1]]), 2))
      for (g in j) {
        summed_oa <- summed_oa + as.matrix(crossprod(vi_alpha_decomp[, g]))
        summed_alpha_outer <- summed_alpha_outer + as.matrix(tcrossprod(vi_alpha_mean[g]))
      }
      store_oa[[counter_j]] <- summed_oa
      store_alpha_outer[[counter_j]] <- summed_alpha_outer
      counter_j <- counter_j + 1
    }
    return(list(outer_alpha = as.matrix(store_oa), alpha_mu_outer = as.matrix(store_alpha_outer)))
  }

  dta <- data.frame(y = 1, g = letters, g2 = LETTERS, g3 = 1:26, a = rnorm(26), x = rnorm(26), z = rnorm(26))
  formula <- y ~ x + (1 + z | g) + (1 + z | g2) + (0 + a | g3)

  mk_Z <- lme4::mkReTrms(lme4::findbars(formula), model.frame(lme4::subbars(formula), data = dta), reorder.terms = FALSE, reorder.vars = FALSE)

  breaks_for_RE <- c(0, cumsum(diff(mk_Z$Gp)))
  d_j <- c(2, 2, 1)
  # Number of GROUPs for each random effect.
  g_j <- diff(mk_Z$Gp) / d_j


  outer_alpha_RE_positions <- mapply(d_j, g_j, breaks_for_RE[-length(breaks_for_RE)], SIMPLIFY = FALSE, FUN = function(a, b, m) {
    split(m + seq(1, a * b), rep(1:b, each = a))
  })

  vi_alpha_var <- drop0(rWishart(n = 1, df = nrow(mk_Z$Zt) + 10, Sigma = diag(nrow(mk_Z$Zt)))[, , 1])
  vi_alpha_chol <- as(drop0((chol(vi_alpha_var))), "generalMatrix")
  expect_equal(as.matrix(t(vi_alpha_chol) %*% vi_alpha_chol), as.matrix(vi_alpha_var))

  vi_alpha_mean <- Matrix(rnorm(nrow(vi_alpha_var)))

  legacy_method <- loop_outer_alpha(vi_alpha_mean, vi_alpha_chol, outer_alpha_RE_positions)
  legacy_method <- mapply(legacy_method[[1]], legacy_method[[2]], SIMPLIFY = FALSE, FUN = function(a, b) {
    a + b
  })

  cpp_method <- calculate_expected_outer_alpha(L = vi_alpha_chol, alpha_mu = as.vector(vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
  cpp_method <- cpp_method$outer_alpha


  comp <- mapply(legacy_method, cpp_method, FUN = function(i, j) {
    all.equal(as.vector(i), as.vector(j))
  })

  expect_equal(unlist(legacy_method), unlist(cpp_method))
})
