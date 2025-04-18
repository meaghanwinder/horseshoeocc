#' Credible interval plots for occupancy model derived parameters fit using the regularized horseshoe prior
#'
#' @description Plot method for objects of class \code{summary.horseshoeocc_derived}.
#'
#' @param sum_fit_derived an object of class \code{summary.horseshoeocc_derived}, typically the
#' result of a call to \code{\link{summary.horseshoeocc_derived}}.
#' @param which the parameters for which to generate trace plots. The value can
#' be either \code{"psi"} or \code{"p"}; defaults to \code{"psi"}.
#' @param median logical; defaults to \code{TRUE}. If \code{FALSE}, the point
#' corresponds to the mean rather than the median.
#' @param hdi logical; defaults to \code{TRUE}. If \code{FALSE}, the quantile based
#' credibility intervals are used.
#' @param equal optional value to compare to each credibility interval.
#' @param ... Other arguments.
#'
#' @returns A ggplot (if \code{"psi"}) or list of ggplot objects (if \code{"p"}).
#' The ggplot objects are credible interval plots of the estimated occupancy or
#' detection probabilities.
#'
#' @examples
#' # continuing the summary.horseshoe_derived example:
#' ## plot HDI credibility intervals for psi
#' plot(derived_summary, "psi")
#'
#' ## list of plots of quantile based credibility intervals for p
#' p_plots <- plot(derived_summary, "p", median = F, hdi = F)
#' p_plots[[1]]
#'
#' @export
#'
plot.summary.horseshoeocc_derived <- function(sum_fit_derived, which = c("psi"), median = TRUE, hdi = TRUE, equal = NULL, ...){
  psi_sum <- sum_fit_derived$psi
  p_sum <- sum_fit_derived$p

  if(which != "psi" & which != "p"){
    stop("The argument 'which' must take on the value of 'psi' or 'p'.")
  }

  if(which == "psi"){
    plot_df <- psi_sum %>%
      as.data.frame() %>%
      dplyr::mutate(var = stringr::str_extract_all(rownames(.), "[:digit:]+"),
                    var = as.numeric(var))
    if(hdi == T){
      title =  paste0('95% HDI credibility intervals for ', which)
    } else if(hdi == F){
      title =  paste0('95% credibility intervals for ', which)
    }
    xlab = "psi"

    if(hdi == T){
      out <- plot_df %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = `HDI 2.5%`, y = var), shape = "|", size = 3) +
        ggplot2::geom_point(ggplot2::aes(x = `HDI 97.5%`, y = var), shape = "|", size = 3) +
        ggplot2::geom_segment(ggplot2::aes(x = `HDI 2.5%`, y = var, xend = `HDI 97.5%`, yend = var)) +
        # ggplot2::geom_point(ggplot2::aes(y = var, x = Mean)) +
        ggplot2::labs(title = title,
                      x = xlab,
                      y = NULL) +
        # ggplot2::scale_y_continuous(breaks = seq(1, nrow(psi_sum))) +
        ggplot2::theme_bw()
      if(!is.null(equal)){
        out <- out +
          ggplot2::geom_vline(xintercept = equal)
      }
    } else if(hdi == F){
      out <- plot_df %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = `2.5%`, y = var), shape = "|", size = 3) +
        ggplot2::geom_point(ggplot2::aes(x = `97.5%`, y = var), shape = "|", size = 3) +
        ggplot2::geom_segment(ggplot2::aes(x = `2.5%`, y = var, xend = `97.5%`, yend = var)) +
        # ggplot2::geom_point(ggplot2::aes(y = var, x = Mean)) +
        ggplot2::labs(title = title,
                      x = xlab,
                      y = NULL) +
        # ggplot2::scale_y_continuous(breaks = seq(1, nrow(psi_sum))) +
        ggplot2::theme_bw()

      if(!is.null(equal)){
        out <- out +
          ggplot2::geom_vline(xintercept = equal)
      }

    }

    if(median == T){
      out <- out +
        ggplot2::geom_point(ggplot2::aes(y = var, x = Median))
    } else if (median == F){
      out <- out +
        ggplot2::geom_point(ggplot2::aes(y = var, x = Mean))
    }
  }


  if(which == "p"){
    plot_df <- p_sum %>%
      as.data.frame() %>%
      dplyr::mutate(site = stringr::str_extract_all(rownames(.), "[:digit:]+(?=\\.)"),
                    samp = stringr::str_extract_all(rownames(.), "(?<=\\.)[:digit:]+"),
                    site = as.numeric(site),
                    samp = as.numeric(samp))
    xlab = "p"
    out <- list()
    if(hdi == T){
      for(i in 1:nrow(psi_sum)){
        title =  paste0('95% HDI credibility intervals for ', which, ' at site ', i)

        out[[i]] <- plot_df %>%
          dplyr::filter(site == i) %>%
          ggplot2::ggplot() +
          ggplot2::geom_point(ggplot2::aes(x = `HDI 2.5%`, y = samp), shape = "|", size = 3) +
          ggplot2::geom_point(ggplot2::aes(x = `HDI 97.5%`, y = samp), shape = "|", size = 3) +
          ggplot2::geom_segment(ggplot2::aes(x = `HDI 2.5%`, y = samp, xend = `HDI 97.5%`, yend = samp)) +
          # ggplot2::geom_point(ggplot2::aes(y = samp, x = Mean)) +
          ggplot2::labs(title = title,
                        x = xlab,
                        y = NULL) +
          ggplot2::xlim(0, 1) +
          ggplot2::theme_bw()
        if(!is.null(equal)){
          out <- out +
            ggplot2::geom_vline(xintercept = equal)
        }

        if(median == T){
          out[[i]] <- out[[i]] +
            ggplot2::geom_point(ggplot2::aes(y = samp, x = Median))
        } else if (median == F){
          out[[i]] <- out[[i]] +
            ggplot2::geom_point(ggplot2::aes(y = samp, x = Mean))
        }
      }
    } else if(hdi == F){
      for(i in 1:nrow(psi_sum)){
        title =  paste0('95% credibility intervals for ', which, ' at site ', i)

        out[[i]] <- plot_df %>%
          dplyr::filter(site == i) %>%
          ggplot2::ggplot() +
          ggplot2::geom_point(ggplot2::aes(x = `2.5%`, y = samp), shape = "|", size = 3) +
          ggplot2::geom_point(ggplot2::aes(x = `97.5%`, y = samp), shape = "|", size = 3) +
          ggplot2::geom_segment(ggplot2::aes(x = `2.5%`, y = samp, xend = `97.5%`, yend = samp)) +
          # ggplot2::geom_point(ggplot2::aes(y = samp, x = Mean)) +
          ggplot2::labs(title = title,
                        x = xlab,
                        y = NULL) +
          ggplot2::xlim(0, 1) +
          ggplot2::theme_bw()
        if(!is.null(equal)){
          out <- out +
            ggplot2::geom_vline(xintercept = equal)
        }

        if(median == T){
          out[[i]] <- out[[i]] +
            ggplot2::geom_point(ggplot2::aes(y = samp, x = Median))
        } else if (median == F){
          out[[i]] <- out[[i]] +
            ggplot2::geom_point(ggplot2::aes(y = samp, x = Mean))
        }
      }
    }
  }

  return(out)
}
