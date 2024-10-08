% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MAVB_functions.R
\name{MAVB}
\alias{MAVB}
\title{Perform MAVB after fitting vglmer}
\usage{
MAVB(object, samples, verbose = FALSE, var_px = Inf)
}
\arguments{
\item{object}{Model fit using \code{vglmer}.}

\item{samples}{Number of samples to draw.}

\item{verbose}{Show progress in drawing samples.}

\item{var_px}{Variance of working prior for marginal augmentation. Default
(\code{Inf}) is a flat, improper, prior.}
}
\value{
This function returns a matrix with \code{samples} rows and columns
  for each fixed and random effect.
}
\description{
Given a model estimated using \code{vglmer}, this function
  performs marginally augmented variational Bayes (MAVB) to improve the
  approximation quality.
}
\details{
This function returns the improved estimates of the
  \emph{parameters}. To use MAVB when generating predictions, one should use
  \link{predict_MAVB}. At present, MAVB is only enabled for binomial models.
}
\references{
Goplerud, Max. 2022. "Fast and Accurate Estimation of Non-Nested Binomial
Hierarchical Models Using Variational Inference." \emph{Bayesian Analysis}. 17(2):
623-650.
}
