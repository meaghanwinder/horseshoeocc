% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.summary.horseshoeocc.R
\name{plot.summary.horseshoeocc}
\alias{plot.summary.horseshoeocc}
\title{Credible interval plots for occupancy model parameters fit using the regularized horseshoe prior}
\usage{
\method{plot}{summary.horseshoeocc}(fit_sum, which = c("beta"), median = TRUE, hdi = TRUE, equal = NULL, ...)
}
\arguments{
\item{fit_sum}{an object of class \code{summary.horseshoeocc}, typically the
result of a call to \code{\link{summary.horseshoeocc}}.}

\item{which}{the parameters for which to generate credibility interval plots. The value can
be a specific parameter (i.e. \code{"beta[1]"}), or you can plot the credibility intervals
for all alpha, beta, lambda, or kappa coefficients using \code{"alpha"},
\code{"beta"}, \code{"lambda"}, or \code{"kappa"}, respectively; defaults to \code{"beta"}.}

\item{median}{logical; defaults to \code{TRUE}. If \code{FALSE}, the point
corresponds to the mean rather than the median.}

\item{hdi}{logical; defaults to \code{TRUE}. If \code{FALSE}, the quantile based
credibility intervals are used.}

\item{equal}{optional value to compare to each credibility interval.}

\item{...}{Other arguments.}
}
\value{
Returns a ggplot object that is a credible interval plot of the
chosen model parameter(s).
}
\description{
Plot method for objects of class \code{summary.horseshoeocc}.
}
\details{
You may only select one group of parameters (i.e. \code{"beta"},
\code{"alpha"}, \code{"lambda"}, or \code{"kappa"}) to plot at a time.
}
\examples{
# continuing the summary.horseshoeocc example:
## plot HDI credibility intervals for beta
plot(fit_summary)

## plot quantile based credibility intervals for alpha
plot(fit_summary, "alpha")

## plot the HDI credibility interval for meff
plot(fit_summary, "meff")

## plot non-zero beta credible intervals
plot(fit_summary, c("beta0", "beta[1]", "beta[2]"), equal = 0)

}
