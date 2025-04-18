% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.horseshoeocc_derived.R
\name{summary.horseshoeocc_derived}
\alias{summary.horseshoeocc_derived}
\alias{print.summary.horseshoeocc_derived}
\title{Summarize occupancy model derived parameters fit using the regularized horseshoe prior}
\usage{
\method{summary}{horseshoeocc_derived}(fit_derived)

\method{print}{summary.horseshoeocc_derived}(sum_fit_derived, hdi = TRUE)
}
\arguments{
\item{fit_derived}{an object of class \code{horseshoeocc_derived}, typically a result of
a calculating the derived parameters from the model fit using \code{\link{derived_parameters}}.}

\item{sum_fit_derived}{An object of class \code{summary.horseshoeocc}, typically the
result of a call to \code{summary.horseshoeocc}}

\item{hdi}{logical; defaults to \code{TRUE}. If \code{FALSE}, the quantile based
credibility intervals are used.}
}
\value{
\code{summary()} computes and returns an object of class \code{list}
containing the elements:
\itemize{
  \item \code{psi} \verb{ }matrix containing posterior summary measures for \code{psi}.

  \item \code{p} \verb{ }matrix containing posterior summary measures for \code{p}.

  \item \code{model} \verb{ }object of class \code{list} containing:
    \itemize{
      \item \code{det_model} \verb{ }formula describing the sample-level model.
      \item \code{occ_model} \verb{ }formula describing the site-level model.
    }

  \item \code{opts} \verb{ } object of class \code{list} containing the MCMC options:
  \code{niter}, \code{nchain}, \code{nburnin}, and \code{thin}.
}
}
\description{
Summary method for objects of class \code{horseshoeocc_derived}.
}
\details{
\code{print.summary.horseshoeocc_derived} formats output from
\code{summary.horseshoeocc_derived}. The output includes the model used and
summaries for both occupancy and detection probabilities, as well as MCMC
information.
}
\examples{
# continuing the derived_parameters example:
derived_summary <- summary(fit_derived)

print(derived_summary)
print(derived_summary, hdi = F)


}
