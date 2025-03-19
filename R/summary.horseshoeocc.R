#' Summarize occupancy model fit using the regularized horseshoe prior
#'
#' @description Summary method for objects of class \code{horseshoeocc}.
#'
#' @param fit an object of class \code{horseshoeocc}, typically a result of
#' a model fit using \code{\link{horseshoeocc}}.
#'
#' @returns \code{summary()} computes and returns an object of class \code{list}
#' containing the elements:
#' \itemize{
#'   \item \code{mcmc_summary} \verb{ }matrix containing posterior summary measures
#'   calculated for the model parameters.
#'
#'   \item \code{model} \verb{ }object of class \code{list} containing:
#'     \itemize{
#'       \item \code{det_model} \verb{ }formula describing the sample-level model.
#'       \item \code{occ_model} \verb{ }formula describing the site-level model.
#'     }
#'
#'   \item \code{meff} \verb{ }matrix containing posterior summary measures calculated
#'   for the effective number of coefficients in the model.
#'
#'   \item \code{opts} \verb{ } object of class \code{list} containing the MCMC options:
#'   \code{niter}, \code{nchain}, \code{nburnin}, and \code{thin}.
#' }
#'
#'
#'
#' @examples
#' # continuing the horseshoeocc example:
#' fit_summary <- summary(fit)
#'
#' print(fit_summary, hdi = T)
#' print(fit_summary, hdi = F)
#'
#' @export
#'
summary.horseshoeocc <- function(fit){
  mcmc <- object$mcmc
  model <- object$model

  # user has the ability to remove using burnin and thin in horseshoeocc()
  warmup = 0; thin = 1

  # convert to coda for summary
  fit_warmup <- lapply(mcmc, function(x) x[(warmup+1):nrow(x),])
  coda_samples <- coda::as.mcmc.list(lapply(fit_warmup, function(x) coda::as.mcmc(
    x, start = warmup+1, end = nrow(mcmc), thin = thin
  )))

  sum <- summary(coda_samples)
  params <- dimnames(sum$statistics)[[1]]
  # clean up parameter names
  tmp_kappa <- params[which(stringr::str_detect(params, "kappa"))] %>%
    stringr::str_remove(., "t")
  tmp_lambda <- params[which(stringr::str_detect(params, "lambda"))] %>%
    stringr::str_remove(., "t")
  tmp_meff <- params[which(stringr::str_detect(params, "meff"))] %>%
    stringr::str_remove(., "t")
  params[which(stringr::str_detect(params, "kappa"))] <- tmp_kappa
  params[which(stringr::str_detect(params, "lambda"))] <- tmp_lambda
  params[which(stringr::str_detect(params, "meff"))] <- tmp_meff

  tmp_sum <- cbind(sum$statistics, sum$quantiles)

  # get r hat / n_eff
  mat <- matrix(NA, nrow = nrow(tmp_sum), ncol = 3)
  colnames(mat) <- c("Rhat", "ess_bulk", "ess_tail")
  for(i in 1:nrow(tmp_sum)){
    tmp <- sapply(mcmc, function(x) x[,i])
    mat[i,] <- c(rstan::Rhat(tmp), rstan::ess_bulk(tmp), rstan::ess_tail(tmp))
  }

  # get hdi interval
  all_samps <- do.call("rbind", coda_samples)
  hdi_2.5 <- apply(all_samps, 2, HDInterval::hdi)[1, ]
  hdi_97.5 <- apply(all_samps, 2, HDInterval::hdi)[2, ]

  # get median
  median <- apply(all_samps, 2, median)

  # out
  out <- cbind(tmp_sum,
               Median = median,
               "HDI 2.5%" = hdi_2.5,
               "HDI 97.5%" = hdi_97.5,
               mat)

  meff <- out %>%
    data.frame() %>%
    dplyr::filter(rownames(.) == "mefft")

  colnames(meff) <- colnames(out)

  out <- list(mcmc_summary = out,
              model = model,
              meff = meff,
              opts = object$opts)

  class(out) <- c('summary.horseshoeocc', class(out))

  return(out)
}

#' Printing summary of occupancy model fit using the regularized horseshoe prior
#'
#' @details
#' \code{print.summary.horseshoeocc} formats output from \code{summary.horseshoeocc}.
#' The output includes the model used and coefficient summaries for both occupancy
#' and detection. Additionally, the output includes a summary of the effective
#' number of beta coefficients in the model as well the MCMC information.
#'
#' @param fit_sum an object of class \code{summary.horseshoeocc}, typically the
#' result of a call to \code{summary.horseshoeocc}.
#' @param hdi logical; defaults to \code{TRUE}. If \code{FALSE}, the quantile based
#' credibility intervals are used.
#'
#' @export
#'
#' @rdname summary.horseshoeocc
#'
print.summary.horseshoeocc <- function(fit_sum, hdi = TRUE){
  occ_model <- x$model$occ_model
  det_model <- x$model$det_model

  summ <- x$mcmc_summary

  if(hdi == T){
    summ <- summ %>%
      round(3) %>%
      data.frame() %>%
      dplyr::select(-Naive.SE, -Time.series.SE,
                    -Rhat, -ess_bulk, -ess_tail,
                    -X25., -X50., -X75.) %>%
      dplyr::filter(rownames(.) != "deviance",
                    !stringr::str_detect(rownames(.), "kappa"),
                    !stringr::str_detect(rownames(.), "lambda"),
                    rownames(.) != "tau") %>%
      dplyr::select(-X2.5., -X97.5.)

    colnames(summ) <- c("Mean",
                        "SD",
                        "Median",
                        "HDI 2.5%",
                        "HDI 97.5%")
  } else if(hdi == F){
    summ <- summ %>%
      round(3) %>%
      data.frame() %>%
      dplyr::select(-Naive.SE, -Time.series.SE,
                    -Rhat, -ess_bulk, -ess_tail,
                    -X25., -X50., -X75.) %>%
      dplyr::filter(rownames(.) != "deviance",
                    !stringr::str_detect(rownames(.), "kappa"),
                    !stringr::str_detect(rownames(.), "lambda"),
                    rownames(.) != "tau") %>%
      dplyr::select(-HDI.2.5., -HDI.97.5.)

    colnames(summ) <- c("Mean",
                        "SD",
                        "2.5%",
                        "97.5%",
                        "Median")
  }

  beta_intercept <- summ %>%
    dplyr::filter(stringr::str_detect(rownames(.), "beta")) %>%
    dplyr::filter(rownames(.) == "beta0") %>%
    dplyr::select(Mean, SD, Median, dplyr::everything())

  beta_tbl <- summ %>%
    dplyr::filter(stringr::str_detect(rownames(.), "beta")) %>%
    dplyr::filter(!rownames(.) == "beta0") %>%
    dplyr::select(Mean, SD, Median, dplyr::everything())

  beta_tbl <- rbind(beta_intercept,
                    beta_tbl)

  alpha_intercept <- summ %>%
    dplyr::filter(stringr::str_detect(rownames(.), "alpha")) %>%
    dplyr::filter(rownames(.) == "alpha0") %>%
    dplyr::select(Mean, SD, Median, dplyr::everything())

  alpha_tbl <- summ %>%
    dplyr::filter(stringr::str_detect(rownames(.), "alpha")) %>%
    dplyr::filter(!rownames(.) == "alpha0") %>%
    dplyr::select(Mean, SD, Median, dplyr::everything())

  alpha_tbl <- rbind(alpha_intercept,
                     alpha_tbl)

  meff_tbl <- summ %>%
    dplyr::filter(rownames(.) == "mefft") %>%
    dplyr::select(Mean, SD, Median, dplyr::everything())

  # print
  cat("Model for occupancy probability:\n")
  print(occ_model)
  cat("\n")

  cat("Effective number of beta coefficients:\n")
  print(meff_tbl)
  cat("\n")

  cat("Beta coefficient summary:\n")
  print(beta_tbl)
  cat("\n")

  cat("Model for detection probability:\n")
  print(det_model)
  cat("\n")

  cat("Alpha coefficient summary:\n")
  print(alpha_tbl)
  cat("\n")

  cat(paste("Model fit using", x$opts$nchain, "chains of", x$opts$niter, "iterations\n"))
  cat(paste("with", x$opts$nburnin, "samples discarded as burn-in;\n"))
  cat(paste("samples were thinned by", x$opts$thin))
}
