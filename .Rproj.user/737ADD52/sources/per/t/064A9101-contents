#' Credible interval plots for occupancy model parameters fit using the regularized horseshoe prior
#'
#' @description Plot method for objects of class \code{summary.horseshoeocc}.
#'
#' @details
#' You may only select one group of parameters (i.e. \code{"beta"},
#' \code{"alpha"}, \code{"lambda"}, or \code{"kappa"}) to plot at a time.
#'
#' @param fit_sum an object of class \code{summary.horseshoeocc}, typically the
#' result of a call to \code{\link{summary.horseshoeocc}}.
#' @param which the parameters for which to generate credibility interval plots. The value can
#' be a specific parameter (i.e. \code{"beta[1]"}), or you can plot the credibility intervals
#' for all alpha, beta, lambda, or kappa coefficients using \code{"alpha"},
#' \code{"beta"}, \code{"lambda"}, or \code{"kappa"}, respectively; defaults to \code{"beta"}.
#' @param median logical; defaults to \code{TRUE}. If \code{FALSE}, the point
#' corresponds to the mean rather than the median.
#' @param hdi logical; defaults to \code{TRUE}. If \code{FALSE}, the quantile based
#' credibility intervals are used.
#' @param equal optional value to compare to each credibility interval.
#' @param ... Other arguments.
#'
#' @returns Returns a ggplot object that is a credible interval plot of the
#' chosen model parameter(s).
#'
#' @examples
#' # continuing the summary.horseshoeocc example:
#' ## plot HDI credibility intervals for beta
#' plot(fit_summary)
#'
#' ## plot quantile based credibility intervals for alpha
#' plot(fit_summary, "alpha")
#'
#' ## plot the HDI credibility interval for meff
#' plot(fit_summary, "meff")
#'
#' ## plot non-zero beta credible intervals
#' plot(fit_summary, c("beta0", "beta[1]", "beta[2]"), equal = 0)
#'
#' @export
#'
plot.summary.horseshoeocc <- function(fit_sum, which = c("beta"), median = TRUE, hdi = TRUE, equal = NULL, ...){
  x <- fit_sum$mcmc
  if(all(length(which) == 1 & which == "alpha") |
     all(length(which) == 1 & which == "beta") |
     all(length(which) == 1 & which == "lambda") |
     all(length(which) == 1 & which == "kappa") |
     all(length(which) == 1 & which == "meff")){
    if(all(is.numeric(which))){
      plot_df <- x %>%
        as.data.frame() %>%
        dplyr::slice(which) %>%
        dplyr::mutate(var = stringr::str_extract(rownames(.), "[:digit:]+")) %>%
        dplyr::mutate(var = as.numeric(var))

      if(hdi == T){
        title =  paste0('95% HDI credibility intervals for ', rownames(plot_df))
      } else if(hdi == F){
        title =  paste0('95% credibility intervals for ', rownames(plot_df))
      }
    } else if(all(is.character(which))){
      plot_df <- x %>%
        as.data.frame() %>%
        dplyr::filter(stringr::str_detect(rownames(.), which)) %>%
        dplyr::mutate(var = stringr::str_extract(rownames(.), "[:digit:]+")) %>%
        dplyr::mutate(var = as.numeric(var)) %>%
        dplyr::mutate(var = dplyr::case_when(is.na(var) ~ 1,
                                             TRUE ~ var))

      if(hdi == T){
        title =  paste0('95% HDI credibility intervals for ', which)
      } else if(hdi == F){
        title =  paste0('95% credibility intervals for ', which)
      }
    }

  } else if(length(which) == 1){
    if(all(is.numeric(which))){
      plot_df <- x %>%
        as.data.frame() %>%
        dplyr::slice(which) %>%
        dplyr::mutate(var = stringr::str_extract(rownames(.), "[:digit:]+")) %>%
        dplyr::mutate(var = as.numeric(var))

      if(hdi == T){
        title =  paste0('95% HDI credibility intervals for ', paste(rownames(plot_df), collapse = " and "))
      } else if(hdi == F){
        title =  paste0('95% credibility intervals for ', paste(rownames(plot_df), collapse = " and "))
      }
    } else if (all(is.character(which))){
      plot_df <- x %>%
        as.data.frame() %>%
        dplyr::mutate(param = rownames(.)) %>%
        dplyr::filter(param == which) %>%
        dplyr::select(-param) %>%
        dplyr::mutate(var = stringr::str_extract(rownames(.), "[:digit:]+")) %>%
        dplyr::mutate(var = as.numeric(var))

      if(hdi == T){
        title =  paste0('95% HDI credibility intervals for ', paste(rownames(plot_df), collapse = " and "))
      } else if(hdi == F){
        title =  paste0('95% credibility intervals for ', paste(rownames(plot_df), collapse = " and "))
      }
    }


  } else if(length(which) > 1){
    if(all(is.numeric(which))){
      plot_df <- x %>%
        as.data.frame() %>%
        dplyr::slice(which) %>%
        dplyr::mutate(var = stringr::str_extract(rownames(.), "[:digit:]+")) %>%
        dplyr::mutate(var = as.numeric(var))

      if(hdi == T){
        title =  paste0('95% HDI credibility intervals for ', paste(rownames(plot_df), collapse = " and "))
      } else if(hdi == F){
        title =  paste0('95% credibility intervals for ', paste(rownames(plot_df), collapse = " and "))
      }
    } else if(all(is.character(which))){
      plot_df <- x %>%
        as.data.frame() %>%
        dplyr::mutate(param = rownames(.)) %>%
        dplyr::filter(param %in% which) %>%
        dplyr::select(-param) %>%
        dplyr::mutate(var = stringr::str_extract(rownames(.), "[:digit:]+")) %>%
        dplyr::mutate(var = as.numeric(var))

      if(hdi == T){
        title =  paste0('95% HDI credibility intervals for ', paste(which, collapse = " and "))
      } else if(hdi == F){
        title =  paste0('95% HDI credibility intervals for ', paste(which, collapse = " and "))
      }
    }

  }


  if(hdi == T){
    out <- plot_df %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x = `HDI 2.5%`, y = var), shape = "|", size = 3) +
      ggplot2::geom_point(ggplot2::aes(x = `HDI 97.5%`, y = var), shape = "|", size = 3) +
      ggplot2::geom_segment(ggplot2::aes(x = `HDI 2.5%`, y = var, xend = `HDI 97.5%`, yend = var)) +
      # ggplot2::geom_point(ggplot2::aes(y = var, x = Mean)) +
      ggplot2::labs(title = title,
                    x = NULL,
                    y = NULL) +
      ggplot2::theme_bw()
    if(!is.null(equal)){
      out <- out +
        ggplot2::geom_vline(xintercept = equal, color = "grey")
    }


  } else if(hdi == F){
    out <- plot_df %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x = `2.5%`, y = var), shape = "|", size = 3) +
      ggplot2::geom_point(ggplot2::aes(x = `97.5%`, y = var), shape = "|", size = 3) +
      ggplot2::geom_segment(ggplot2::aes(x = `2.5%`, y = var, xend = `97.5%`, yend = var)) +
      # ggplot2::geom_point(ggplot2::aes(y = var, x = Mean)) +
      ggplot2::labs(title = title,
                    x = NULL,
                    y = NULL) +
      ggplot2::theme_bw()

    if(!is.null(equal)){
      out <- out +
        ggplot2::geom_vline(xintercept = equal, color = "grey")
    }

  }

  if(median == T){
    out <- out +
      ggplot2::geom_point(ggplot2::aes(y = var, x = Median))
  } else if(median == F){
    out <- out +
      ggplot2::geom_point(ggplot2::aes(y = var, x = Mean))
  }

  return(out)
}
