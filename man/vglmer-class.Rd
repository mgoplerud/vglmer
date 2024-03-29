% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_functions.R
\name{vglmer-class}
\alias{vglmer-class}
\alias{fixef.vglmer}
\alias{sigma.vglmer}
\alias{ranef.vglmer}
\alias{coef.vglmer}
\alias{vcov.vglmer}
\alias{fitted.vglmer}
\alias{print.vglmer}
\alias{summary.vglmer}
\alias{formula.vglmer}
\alias{format_vglmer}
\alias{format_glmer}
\alias{ELBO}
\title{Generic Functions after Running vglmer}
\usage{
\method{fixef}{vglmer}(object, ...)

\method{sigma}{vglmer}(object, ...)

\method{ranef}{vglmer}(object, ...)

\method{coef}{vglmer}(object, ...)

\method{vcov}{vglmer}(object, ...)

\method{fitted}{vglmer}(object, ...)

\method{print}{vglmer}(x, ...)

\method{summary}{vglmer}(object, display_re = TRUE, ...)

\method{formula}{vglmer}(x, form = "original", ...)

format_vglmer(object)

format_glmer(object)

ELBO(object, type = c("final", "trajectory"))
}
\arguments{
\item{object}{Model fit using \code{vglmer}.}

\item{...}{Not used; included to maintain compatibility with existing
methods.}

\item{x}{Model fit using \code{vglmer}.}

\item{display_re}{Default (\code{TRUE}) prints a summary of the
random effects alongside the fixed effects.}

\item{form}{Describes the type of formula to report:
\code{"original"} returns the user input, \code{"fe"} returns the fixed
effects only, \code{"re"} returns the random effects only.}

\item{type}{Default (\code{"final"}) gives the ELBO at convergence.
\code{"trajectory"} gives the ELBO estimated at each iteration. This is
used to assess model convergence.}
}
\value{
The functions here return a variety of different objects depending on
  the specific function. "Details" describes the behavior of each one. Their
  output is similar to the typical behavior for the corresponding generic
  functions.
}
\description{
\code{vglmer} uses many standard methods from \code{lm} and \code{lme4} with
limited changes. These provide summaries of the estimated variational
distributions.
}
\details{
The accompanying functions are briefly described below. 

\code{coef} and \code{vcov} return the mean and variance of the fixed effects
(\eqn{\beta}). \code{fixef} returns the mean of the fixed effects.

\code{ranef} extracts the random effects (\eqn{\alpha}) in a similar,
although slightly different format, to \code{lme4}. It includes the estimated
posterior mean and variance in a list of data.frames with one entry per
random effect \eqn{j}.

\code{fitted} extracts the estimated expected \emph{linear predictor}, i.e.
\eqn{E_{q(\theta)}[x_i^T \beta + z_i^T \alpha]} at convergence.

\code{summary} reports the estimates for all fixed effects as in \code{lm} as
well as some summaries of the random effects (if \code{display_re=TRUE}).

\code{format_vglmer} collects the mean and variance of the fixed and random
effects into a single data.frame. This is useful for examining all of the
posterior estimates simultaneously. \code{format_glmer} converts an object
estimated with \code{[g]lmer} into a comparable format.

\code{ELBO} extracts the ELBO from the estimated model. \code{type} can be
set equal to \code{"trajectory"} to get the estimated ELBO at each iteration
and assess convergence.
  
\code{sigma} extracts the square root of the posterior mode of
\eqn{q(\sigma^2)} if a linear model is used.

\code{formula} extracts the formula associated with the \code{vglmer} object.
By default, it returns the formula provided. The fixed and random effects
portions can be extracted separately using the \code{form} argument.
}
