#' Summarize occupancy model derived parameters fit using the regularized horseshoe prior
#'
#' @description Summary method for objects of class \code{horseshoeocc_derived}.
#'
#' @param fit_derived an object of class \code{horseshoeocc_derived}, typically a result of
#' a calculating the derived parameters from the model fit using \code{\link{derived_parameters}}.
#'
#' @returns \code{summary()} computes and returns an object of class \code{list}
#' containing the elements:
#' \itemize{
#'   \item \code{psi} \verb{ }matrix containing posterior summary measures for \code{psi}.
#'
#'   \item \code{p} \verb{ }matrix containing posterior summary measures for \code{p}.
#'
#'   \item \code{model} \verb{ }object of class \code{list} containing:
#'     \itemize{
#'       \item \code{det_model} \verb{ }formula describing the sample-level model.
#'       \item \code{occ_model} \verb{ }formula describing the site-level model.
#'     }
#'
#'   \item \code{opts} \verb{ } object of class \code{list} containing the MCMC options:
#'   \code{niter}, \code{nchain}, \code{nburnin}, and \code{thin}.
#' }
#'
#' @examples
#' # continuing the derived_parameters example:
#' derived_summary <- summary(fit_derived)
#'
#' print(derived_summary)
#' print(derived_summary, hdi = F)
#'
#'
#' @export
#'
summary.horseshoeocc_derived <- function(fit_derived){
  # occupancy probability summary
  psi <- fit_derived$psi

  sum_psi <- summary(psi)
  tmp_psi <- cbind(sum_psi$statistics[, 1:2], sum_psi$quantiles)

  hdi_2.5 <- apply(psi, 2, HDInterval::hdi)[1, ]
  hdi_97.5 <- apply(psi, 2, HDInterval::hdi)[2, ]

  # get median
  median <- apply(psi, 2, median)

  psi_summary <- cbind(tmp_psi,
                       Median = median,
                       "HDI 2.5%" = hdi_2.5,
                       "HDI 97.5%" = hdi_97.5)

  # detection probability summary
  p <- fit_derived$p

  sum_p <- summary(p)
  tmp_p <- cbind(sum_p$statistics[, 1:2], sum_p$quantiles)

  hdi_2.5 <- apply(p, 2, HDInterval::hdi)[1, ]
  hdi_97.5 <- apply(p, 2, HDInterval::hdi)[2, ]

  # get median
  median <- apply(p, 2, median)

  p_summary <- cbind(tmp_p,
                     Median = median,
                     "HDI 2.5%" = hdi_2.5,
                     "HDI 97.5%" = hdi_97.5)

  out <- list(psi = psi_summary,
              p = p_summary,
              model = fit_derived$model,
              opts = fit_derived$opts)

  class(out) <- c('summary.horseshoeocc_derived', class(out))

  return(out)
}

#' Printing summary of occupancy model derived parameters fit using the
#' regularized horseshoe prior
#'
#' @details
#' \code{print.summary.horseshoeocc_derived} formats output from
#' \code{summary.horseshoeocc_derived}. The output includes the model used and
#' summaries for both occupancy and detection probabilities, as well as MCMC
#' information.
#'
#' @param sum_fit_derived An object of class \code{summary.horseshoeocc}, typically the
#' result of a call to \code{summary.horseshoeocc}
#' @param hdi logical; defaults to \code{TRUE}. If \code{FALSE}, the quantile based
#' credibility intervals are used.
#'
#' @export
#'
#' @rdname summary.horseshoeocc_derived
#'
print.summary.horseshoeocc_derived <- function(sum_fit_derived, hdi = TRUE){
  occ_model <- sum_fit_derived$model$occ_model
  det_model <- sum_fit_derived$model$det_model

  # occ estimates
  psi <- sum_fit_derived$psi

  if(hdi == T){
    psi <- psi %>%
      round(3) %>%
      data.frame() %>%
      dplyr::select(-X25., -X50., -X75.) %>%
      dplyr::select(-X2.5., -X97.5.)

    colnames(psi) <- c("Mean",
                        "SD",
                        "Median",
                        "HDI 2.5%",
                        "HDI 97.5%")
  } else if(hdi == F){
    psi <- psi %>%
      round(3) %>%
      data.frame() %>%
      dplyr::select(-X25., -X50., -X75.) %>%
      dplyr::select(-HDI.2.5., -HDI.97.5.)

    colnames(psi) <- c("Mean",
                       "SD",
                       "2.5%",
                       "97.5%",
                       "Median")
  }

  p <- sum_fit_derived$p

  if(hdi == T){
    p <- p %>%
      round(3) %>%
      data.frame() %>%
      dplyr::select(-X25., -X50., -X75.) %>%
      dplyr::select(-X2.5., -X97.5.)

    colnames(psi) <- c("Mean",
                       "SD",
                       "2.5%",
                       "97.5%",
                       "Median")
  } else if(hdi == F){
    p <- p %>%
      round(3) %>%
      data.frame() %>%
      dplyr::select(-X25., -X50., -X75.) %>%
      dplyr::select(-HDI.2.5., -HDI.97.5.)

    colnames(psi) <- c("Mean",
                       "SD",
                       "2.5%",
                       "97.5%",
                       "Median")
  }

  psi_tbl <- psi %>%
    dplyr::select(Mean, SD, Median, dplyr::everything())

  p_tbl <- p %>%
    dplyr::select(Mean, SD, Median, dplyr::everything())

  # print
  cat("Model for occupancy probability:\n")
  print(occ_model)
  cat("\n")

  cat("Occupancy probability summary:\n")
  print(psi_tbl)
  cat("\n")

  cat("Model for detection probability:\n")
  print(det_model)
  cat("\n")

  cat("Detection probability summary:\n")
  print(p_tbl)
  cat("\n")

  cat(paste("Model fit using", sum_fit_derived$opts$nchain, "chains of", sum_fit_derived$opts$niter, "iterations\n"))
  cat(paste("with", sum_fit_derived$opts$nburnin, "samples discarded as burn-in;\n"))
  cat(paste("samples were thinned by", sum_fit_derived$opts$thin))
}
