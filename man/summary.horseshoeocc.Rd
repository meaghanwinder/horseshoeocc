% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.horseshoeocc.R
\name{summary.horseshoeocc}
\alias{summary.horseshoeocc}
\alias{print.summary.horseshoeocc}
\title{Summarize occupancy model fit using the regularized horseshoe prior}
\usage{
\method{summary}{horseshoeocc}(fit)

\method{print}{summary.horseshoeocc}(fit_sum, hdi = TRUE)
}
\arguments{
\item{fit}{an object of class \code{horseshoeocc}, typically a result of
a model fit using \code{\link{horseshoeocc}}.}

\item{fit_sum}{an object of class \code{summary.horseshoeocc}, typically the
result of a call to \code{summary.horseshoeocc}.}

\item{hdi}{logical; defaults to \code{TRUE}. If \code{FALSE}, the quantile based
credibility intervals are used.}
}
\value{
\code{summary()} computes and returns an object of class \code{list}
containing the elements:
\itemize{
  \item \code{mcmc_summary} \verb{ }matrix containing posterior summary measures
  calculated for the model parameters.

  \item \code{model} \verb{ }object of class \code{list} containing:
    \itemize{
      \item \code{det_model} \verb{ }formula describing the sample-level model.
      \item \code{occ_model} \verb{ }formula describing the site-level model.
    }

  \item \code{meff} \verb{ }matrix containing posterior summary measures calculated
  for the effective number of coefficients in the model.

  \item \code{opts} \verb{ } object of class \code{list} containing the MCMC options:
  \code{niter}, \code{nchain}, \code{nburnin}, and \code{thin}.
}
}
\description{
Summary method for objects of class \code{horseshoeocc}.
}
\details{
\code{print.summary.horseshoeocc} formats output from \code{summary.horseshoeocc}.
The output includes the model used and coefficient summaries for both occupancy
and detection. Additionally, the output includes a summary of the effective
number of beta coefficients in the model as well the MCMC information.
}
\examples{
# continuing the horseshoeocc example:
fit_summary <- summary(fit)

print(fit_summary)
print(fit_summary, hdi = F)

}
