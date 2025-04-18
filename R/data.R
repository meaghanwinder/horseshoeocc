#' Simulated occupancy data
#'
#' @description A simulated occupancy data set where there are 100 sites with 10 samples
#' collected per site. The data set is designed to illustrate the use of the
#' regularized horseshoe prior with a mixture prior on the prior guess for the
#' number of relevant coefficients (\code{p0}).
#'
#' This data set includes three true site-level covariates (\code{x1}:\code{x3})
#' that influence occupancy probability, along with 97 spurious site-level covariates
#' (\code{x4}:\code{x100}) to demonstrate shrinkage under the regularized horseshoe prior.
#' Additionally, the data include a single sample-level covariate (\code{w1}) for
#' the detection probability of each sample.
#'
#' @format A list with 6 elements:
#' \describe{
#'   \item{site_df}{A data frame with site information, containing 100 rows and
#'   101 columns. Columns include the site identifier (\code{site}) and the
#'   site-level covariates (\code{x1}:\code{x100}).}
#'   \item{samp_df}{A data frame with sample information, containing 1000 rows
#'   and 4 columns. Columns include site identifier (\code{site}), sample
#'   identifier (\code{samp}), the sample-level covariate (\code{w1}),
#'   and detection/non-detection data (\code{y}, where 1 indicates detection and
#'   0 indicates non-detection).}
#'   \item{z}{A vector containing the true occupancy state for each site (1 for
#'   occupied, 0 for unoccupied).}
#'   \item{x}{A data frame with 100 rows and columns \code{x1}:\code{x100},
#'   representing the site-level covariates for each of the 100 sites.}
#'   \item{w}{A data frame with 1000 rows and one column (\code{w1}), representing
#'   the sample-level covariate for each sample from each site.}
#'   \item{truth}{A list containing true parameter values for the models,
#'   assuming a logit link:
#'   \describe{
#'      \item{\code{beta0}}{Intercept for the occupancy probability model.}
#'      \item{\code{beta}}{Coefficients for the site-level covariates.}
#'      \item{\code{alpha0}}{Intercept for the detection probability model.}
#'      \item{\code{alpha}}{Coefficients for the sample-level covariate.}
#'   }}
#'
#' }
#'
#'
#' @examples
#' data(sim_data100)
#'
#' @name sim_data100
#' @docType data
#'
"sim_data100"


#' Koala occupancy data
#'
#' @description We use koala occupancy data from Mornington Peninsula, Victoria,
#' Australia collected by
#' \href{https://doi.org/10.1007/s10980-023-01620-2}{Whisson et al. (2023)}.
#' Whisson et al. (2023) investigate koala occupancy of 123 sites using acoustic surveys,
#' and found that tree cover (LTC) and road density within 1 km of the
#' recording device (SRD) were associated with changes in site occupancy rates;
#' they estimated that site occupancy increased as tree cover increased,
#' and declined as road density increased. Additionally, Whisson et al. (2023) found
#' that detection probability decreased with increased inclement weather. We use the
#' data from Whisson et al. (2023) and add spurious covariates to create a synthetic data
#' examples with 100 potential occupancy covariates. The synethic data includes both
#' tree cover and road density; the remaining spurious covariates are generated
#' from a standard normal distribution.
#'
#' @details
#' This data set includes the site-level covariates found to be associated with changes
#' in site occupancy probabilities by Whisson et al. (\code{LTC} and \code{SRD}),
#' along with 98 simulated spurious site-level covariates (\code{x3}:\code{x100})
#' to demonstrate shrinkage under the regularized horseshoe prior. Additionally,
#' the data include a single sample-level covariate (\code{weather}) for the detection
#' probability of each sample.
#'
#' @format A list with 2 elements:
#' \describe{
#'   \item{site_df}{A data frame with site information, containing 123 rows and
#'   101 columns. Columns include the site identifier (\code{site}) and the
#'   site-level covariates (\code{LTC}, \code{SRD}, and \code{x3}:\code{x100}).}
#'   \item{samp_df}{A data frame with sample information, containing 911 rows
#'   and 6 columns. Columns include site identifier (\code{site}), sample
#'   identifier (\code{samp}), the sample-level covariate (\code{weather}), in addition
#'   to the two site-level covariates (\code{LTC} and \code{SRD}). The data frame
#'   also includes detection/non-detection data (\code{y}, where 1 indicates detection and
#'   0 indicates non-detection).}
#'
#' }
#'
#'


#' @examples
#' data(koala_data100)
#'
#' @name koala_data100
#' @docType data
"koala_data100"
