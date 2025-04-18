% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derived_parameters.R
\name{derived_parameters}
\alias{derived_parameters}
\title{Compute derived parameters psi and p}
\usage{
derived_parameters(fit)
}
\arguments{
\item{fit}{object of class \code{horseshoeocc}.}
}
\value{
An object of class \code{list} containing the elements:
\itemize{
  \item \code{psi} \verb{ }an object of class \code{mcmc} containing the
  derived posterior samples for the occupancy probabilities (\code{psi}).

  \item \code{p} \verb{ }an object of class \code{mcmc} containing the
  derived posterior samples for the detection probabilities (\code{p}).

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
This function computes the derived parameters \code{psi}
and \code{p}, representing the occupancy and detection probabilities respectively.
}
\examples{
# continuing the horseshoeocc example:
fit_derived <- derived_parameters(fit)

# calculate and print a summary of the derived parameters
derived_summary <- summary(fit_derived)

print(derived_summary)

# plot the credibility intervals for psi
plot(derived_summary, "psi")

# plot the credibility intervals for p
p_plots <- plot(derived_summary, "p")

}
