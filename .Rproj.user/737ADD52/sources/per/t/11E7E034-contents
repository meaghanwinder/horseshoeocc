#' Fit occupancy model using the regularized horseshoe prior
#'
#' @description This function fits a the single-season single species occupancy
#' model described by
#' \href{https://doi.org/10.1890/0012-9658(2002)083[2248:ESORWD]2.0.CO;2}{MacKenzie et al. (2002)}
#' using the regularized horseshoe prior described by
#' \href{https://doi.org/10.1214/17-EJS1337SI}{Piironen and Vehtari (2017)}.
#' Instead of specifying a value for the guess for the number of relevant coefficients (\eqn{p_0}),
#' the user should specify a mixture of truncated Normal distributions as a prior
#' distribution for \eqn{p_0}.
#'
#'
#' @usage
#' fit_horseshoeocc(
#'   data = list(site_df = NULL,
#'               samp_df = NULL),
#'   occ_model = ~.,
#'   det_model = ~1,
#'   priors = list(
#'     beta0 = list(mu0 = 0, sigma0 = 2),
#'     alpha0 = list(mu0 = 0, sigma0 = 2),
#'     alpha = list(mu0 = 0, sigma0 = 2),
#'     c = list(df_slab = 4, scale_slab = 2)
#'   ),
#'   mixture_mean1 = NULL,
#'   mixture_var1 = 1,
#'   mixture_weight1 = 0.5,
#'   mixture_mean2 = NULL,
#'   mixture_var2 = 25,
#'   mixture_weight2 = 0.5,
#'   pseudovar = 4,
#'   nchain = 3,
#'   niter = 1e+05,
#'   nburnin = niter/2,
#'   thin = 5,
#'   seed = NULL
#' )
#'
#' @param data an object of class \code{list} containing the following elements:
#' \itemize{
#'   \item \code{site_df} object of class \code{data.frame} containing site-specific
#'  covariates; \code{site_df} should have exactly one row for each site and have
#'  a column named \code{site}.
#'   \item \code{samp_df} object of class \code{data.frame} containing sample-specific
#'  covariates; \code{samp_df} should have exactly one row for each sample and
#'  have columns named \code{site} and \code{sample}.
#' }
#'
#' @param occ_model formula describing the site-level model. The default \code{~.}
#' includes all columns in \code{site_df} except \code{site}.
#'
#' @param det_model formula describing the sample-level model.
#'
#' @param priors an object of class \code{list} containing the following elements:
#'  \itemize{
#'     \item \code{beta0} \verb{ }object of class \code{list} containing the
#'  following elements:
#'     \itemize{
#'       \item \code{mu0} \verb{ }prior mean for the site-level intercept
#'       \item \code{sigma0} \verb{ }prior standard deviation for the site-level intercept
#'     }
#'     \item \code{alpha0} \verb{ }object of class \code{list} containing the
#'  following elements:
#'     \itemize{
#'       \item \code{mu0} \verb{ }prior mean for sample-level intercept
#'       \item \code{sigma0} \verb{ }prior standard deviation matrix for sample-level intercept
#'     }
#'    \item \code{alpha} \verb{ }object of class \code{list} containing the
#'  following elements:
#'    \itemize{
#'      \item \code{mu0} \verb{ }prior mean for sample-level regression coefficients \cr
#'      \item \code{sigma0} \verb{ }prior standard deviation for sample-level regression
#'      coefficients
#'    }
#'    \item \code{c} \verb{ }object of class \code{list} containing the
#'  following elements:
#'    \itemize{
#'      \item \code{df_slab} \verb{ }prior degrees of freedom for the Student-t slab  \cr
#'      \item \code{scale_slab} \verb{ }prior scale of the Student-t slab
#'    }
#' }
#'
#' @param mixture_mean1 the prior mean of the Normal distribution truncated
#' between 1 and D-1 that serves as the first component of the mixture prior.
#'
#' @param mixture_var1 the prior variance of the Normal distribution truncated
#' between 1 and D-1 that serves as the first component of the mixture prior;
#' defaults to 1.
#'
#' @param mixture_weight1 the weight of the first component of the mixture prior;
#' defaults to 0.5.
#'
#' @param mixture_mean2 the prior mean of the Normal distribution truncated
#' between 1 and D-1 that serves as the second component of the mixture prior.
#'
#' @param mixture_var2 the prior variance of the Normal distribution truncated
#' between 1 and D-1 that serves as the second component of the mixture prior;
#' defaults to 25.
#'
#' @param mixture_weight2 the weight of the second component of the mixture prior;
#' defaults to 0.5.
#'
#' @param pseudovar estimated value for the pseudo-variance
#' \href{https://doi.org/10.1214/17-EJS1337SI}{(Piironen and Vehtari 2017)}; the
#' default of 4 corresponds to a proportion of occupied sites of 0.5. Changing
#' the value of \eqn{\bar{z}} does not have a meaningful impact on the results, so
#' we recommend the default value.
#'
#' @param nchain number of Markov chains; defaults to 3.
#'
#' @param niter number of iterations per chain (including burn-in); defaults to 100,000.
#'
#' @param nburnin length of burn-in, i.e. number of iterations to discard at the
#' beginning. The default is \code{niter/2}, discarding the first half of
#' the simulations.
#'
#' @param thin thinning rate. Must be a positive integer. Set \code{nthin} > 1 to
#' save memory and computation time if \code{niter} is large; defaults to 5.
#'
#' @param seed an optional seed.
#'
#' @returns An object of class \code{list} containing the elements:
#' \itemize{
#'   \item \code{mcmc} \verb{ }object of class \code{list} containing the posterior
#'   samples for each parameter. The number of elements corresponds to \code{nchain}.
#'   \item \code{data} \verb{ }object of class \code{list} containing:
#'     \itemize{
#'       \item \code{x} \verb{ }object of class \code{data.frame} containing the
#'       site-level covariates used in the model for occupancy probability.
#'       \item \code{w} \verb{ }object of class \code{data.frame} containing the
#'       sample-level covariates used in the model for detection probability.
#'       \item \code{site_df} \verb{ }inherited from input: \code{data$site_df}.
#'       \item \code{samp_df} \verb{ }inherited from input: \code{data$samp_df}.
#'     }
#'   \item \code{model} \verb{ }object of class \code{list} containing:
#'     \itemize{
#'       \item \code{det_model} \verb{ }formula describing the sample-level model.
#'       \item \code{occ_model} \verb{ }formula describing the site-level model.
#'     }
#'   \item \code{opts} \verb{ }object of class \code{list} containing the MCMC options:
#'   \code{niter}, \code{nchain}, \code{nburnin}, and \code{thin}.
#' }
#'
#' @examples
#' # simulated dataset
#' data(koala_data100)
#'
#' # fit the model to simulated data
#' fit <- fit_horseshoeocc(data = koala_data100,
#'                         det_model = ~weather,
#'                         occ_model = ~.,
#'                         mixture_mean1 = 2,
#'                         mixture_var1 = 1,
#'                         mixture_mean2 = 2,
#'                         mixture_var2 = 25,
#'                         niter = 1000, # small number of iterations to run quickly
#'                         seed = 123)
#'
#' # posterior summaries of fit
#' fit_summary <- summary(fit)
#'
#' # trace plots for non-zero beta coefficients
#' plot(fit, c("beta0", "beta[1]", "beta[2]"))
#'
#' # credibility intervals for non-zero beta coefficients
#' plot(fit_summary, c("beta0", "beta[1]", "beta[2]"), equal = 0)
#'
#'
#' @export

fit_horseshoeocc <- function(data = list(site_df = NULL,
                                         samp_df = NULL),
                             occ_model = ~.,
                             det_model = ~1,
                             priors = list(
                               beta0 = list(mu0 = 0, sigma0 = 2),
                               alpha0 = list(mu0 = 0, sigma0 = 2),
                               alpha = list(mu0 = 0, sigma0 = 2),
                               c = list(df_slab = 4, scale_slab = 2)),
                             mixture_mean1 = NULL,
                             mixture_var1 = 1,
                             mixture_weight1 = 0.5,
                             mixture_mean2 = NULL,
                             mixture_var2 = 25,
                             mixture_weight2 = 0.5,
                             pseudovar = 4,
                             nchain = 3,
                             niter = 1e+05,
                             nburnin = niter/2,
                             thin = 5,
                             seed = NULL
){
  # set optional seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  det_model <- as.formula(det_model)
  occ_model <- as.formula(occ_model)

  # site_df
  if(is.null(data$site_df)){
    stop("Missing data$site_df.")
  }
  site_df <- data$site_df
  if(!("site" %in% names(data$site_df))){
    stop("There should be a column named 'site' in data$site_df.")
  }
  if(nrow(site_df) != length(unique(site_df$site))){
    stop("Each row of site_df should represent a unique site.")
  }
  nsite <- nrow(site_df)


  # samp_df
  if(is.null(data$samp_df)){
    stop("Missing data$samp_df.")
  }
  samp_df <- data$samp_df
  if(!("site" %in% names(data$samp_df))){
    stop("There should be a column named 'site' in data$samp_df")
  }
  if(!("samp" %in% names(data$samp_df))){
    stop("There should be a column named 'samp' in data$samp_df")
  }
  if(nrow(samp_df) != length(unique(interaction(samp_df$site, samp_df$samp)))){
    stop("Each row of samp_df should represent one sample from a site.")
  }


  # x matrix
  if((ncol(site_df %>% dplyr::select(-site))) < 2){
    stop("You must include at least 2 site-level covariates.")
  }
  if(occ_model == ~.){
    x <- site_df %>%
      dplyr::select(-site)
    x <- model.matrix(occ_model, x)[, -1]
  } else if(occ_model != ~.){
    x_vars <- all.vars(occ_model)
    if(sum(colnames(site_df)[which(colnames(site_df) %in% x_vars)] %in% x_vars) != length(x_vars)){
      stop("You have specified one or more variables in occ_model which do not exist in site_df.")
    }
    x <- site_df %>%
      dplyr::select(which(colnames(.) %in% x_vars))
    x <- model.matrix(occ_model, x)[, -1]
  }

  D <- ncol(x)
  if(apply(x, 2, is.numeric) %>% sum() != D){
    stop("All columns of data$site_df$x must be numeric.")
  }
  corx <- cor(x)[upper.tri(cor(x))] %>%
    abs() %>%
    max()
  if(corx > 0.5){
    warning("Your site-level covariates should be uncorrelated.")
  }
  colmeansx <- apply(x, 2, mean) %>% round(2) %>% sum()
  colsdx <- apply(x, 2, sd) %>% round(2) %>% sum()
  if(colmeansx != 0 | colsdx != D){
    message("Site-level covaritates have been scaled to have a mean of 0 and a standard deviation of 1.")
    x <- scale(x)
  }

  if(mixture_weight1 + mixture_weight2 != 1){
    mixture_weight1_new <- mixture_weight1/(mixture_weight1 + mixture_weight2)
    mixture_weight2_new <- mixture_weight2/(mixture_weight1 + mixture_weight2)
    message(paste("Mixture weights were scaled to sum to 1. mixture_weight1 is now", mixture_weight1_new,
                   "and mixture_weight2 is now", mixture_weight2_new))
    mixture_weight1 <- mixture_weight1_new
    mixture_weight2 <- mixture_weight2_new
  }

  # w matrix
  if(det_model == ~1){
    w <- NULL

    # jags model without covariates for detection
    jags_reghs_occ <- function(){
      # likelihood
      for(i in 1:nrow_site_df){
        logit(psi[i]) <- beta0 + inprod(beta[1:D], x[i, 1:D])
        z[i] ~ dbern(psi[i])
      }

      logit(p) <- alpha0

      for(j in 1:nrow_samp_df){
        y[j] ~ dbern(p*z[site_season_ndx[j]])
      }

      ## p0_mixture prior
      m1_p0 ~ dnorm(mean_m1, 1/var_m1);T(1, D-1)
      m2_p0 ~ dnorm(mean_m2, 1/var_m2);T(1, D-1)
      p0 <- weight_m1*m1_p0 + weight_m2*m2_p0

      # compute quantities for hyper prior on tau
      tau0 <- p0/(D-p0) * sqrtpseudovar / (sqrtnsite)
      tau02 <- tau0 * tau0

      # priors
      beta0 ~ dnorm(beta0mu, 1/beta0var) # intercept
      alpha0 ~ dnorm(alpha0mu, 1/alpha0var) # intercept

      tau ~ dt(0, 1/tau02, 1);T(0,) # global shrinkage
      tau2 <- tau*tau
      c2.rep ~ dgamma(df.slab/2, scale.slab^2*df.slab/2)
      c2 <- 1/c2.rep

      for(j in 1:D){
        lambda[j] ~ dt(0, 1, 1);T(0,)
        lambda2[j] <- lambda[j]*lambda[j]
        lambdat[j] <- (c2 * lambda2[j])/(c2 + tau2 * lambda2[j])
        lambda2t[j] <- lambdat[j] * lambdat[j]
        beta[j] ~ dnorm(0, 1/(tau2 * lambda2t[j]))
      }

      # calculate shrinkage factor
      b <- 1/(1 + nrow_site_df * 1/pseudovar * c2)

      for(j in 1:D){
        kappa[j] <- 1/(1 + nrow_site_df * 1/pseudovar * tau2 * 1 * lambda2[j])
        kappat[j] <- (1-b)*kappa[j] + b
      }
      meff <- sum(1 - kappa[1:D])
      mefft <- (1-b)*meff
    }

    # initialize
    inits = function(){
      list(
        # coefficients
        beta0 = 0,
        beta = rep(0, D),
        alpha0 = 0,
        # horseshoe prior inits
        lambda = rep(1, D),
        tau = .05,
        c2.rep = 1,
        z = rep(1, nsite)
      )
    }

    data_list = list(
      nsite = nrow(site_df),
      D = ncol(x),
      x = x,
      nrow_site_df = nrow(site_df),
      nrow_samp_df = nrow(samp_df),
      site_season_ndx = samp_df$site,
      y = samp_df$y,
      df.slab = priors$c$df_slab,
      scale.slab = priors$c$scale_slab,
      mean_m1 = mixture_mean1,
      mean_m2 = mixture_mean2,
      var_m1 = mixture_var1,
      var_m2 = mixture_var2,
      weight_m1 = mixture_weight1,
      weight_m2 = mixture_weight2,
      beta0mu = priors$beta0$mu0,
      beta0var = (priors$beta0$sigma0)^2,
      alpha0mu = priors$alpha0$mu0,
      alpha0var = (priors$alpha0$sigma0)^2,
      sqrtpseudovar = sqrt(pseudovar),
      pseudovar = pseudovar,
      sqrtnsite = sqrt(nrow(site_df))
    )

    params <- c("alpha0", "beta0", "beta", "tau", "lambdat", "kappat", "mefft")

  } else if(det_model != ~1){
    w_vars <- all.vars(det_model)
    if(sum(colnames(samp_df)[which(colnames(samp_df) %in% w_vars)] %in% w_vars) != length(w_vars)){
      stop("You have specified one or more variables in det_model which do not exist in samp_df.")
    }
    w <- samp_df %>%
      dplyr::select(which(colnames(.) %in% w_vars))
    w <- model.matrix(det_model, w)[, -1] %>%
      matrix(nrow = nrow(samp_df))

    # jags model with covariates for detection
    jags_reghs_occ <- function(){
      # likelihood
      for(i in 1:nrow_site_df){
        logit(psi[i]) <- beta0 + inprod(beta[1:D], x[i, 1:D])
        z[i] ~ dbern(psi[i])
      }

      for(j in 1:nrow_samp_df){
        logit(p[j]) <- alpha0 + inprod(alpha[1:ncol_W], w[j, 1:ncol_W])
        y[j] ~ dbern(p[j]*z[site_season_ndx[j]])
      }

      ## p0_mixture prior
      m1_p0 ~ dnorm(mean_m1, 1/var_m1);T(1, D-1)
      m2_p0 ~ dnorm(mean_m2, 1/var_m2);T(1, D-1)
      p0 <- weight_m1*m1_p0 + weight_m2*m2_p0

      # compute quantities for hyper prior on tau
      tau0 <- p0/(D-p0) * sqrtpseudovar / (sqrtnsite)
      tau02 <- tau0 * tau0

      # priors
      beta0 ~ dnorm(beta0mu, 1/beta0var) # intercept
      alpha0 ~ dnorm(alpha0mu, 1/alpha0var) # intercept

      for(i in 1:ncol_W){
        alpha[i] ~ dnorm(alphamu, 1/alphavar)
      }

      tau ~ dt(0, 1/tau02, 1);T(0,) # global shrinkage
      tau2 <- tau*tau
      c2.rep ~ dgamma(df.slab/2, scale.slab^2*df.slab/2)
      c2 <- 1/c2.rep

      for(j in 1:D){
        lambda[j] ~ dt(0, 1, 1);T(0,)
        lambda2[j] <- lambda[j]*lambda[j]
        lambdat[j] <- (c2 * lambda2[j])/(c2 + tau2 * lambda2[j])
        lambda2t[j] <- lambdat[j] * lambdat[j]
        beta[j] ~ dnorm(0, 1/(tau2 * lambda2t[j]))
      }

      # calculate shrinkage factor
      b <- 1/(1 + nrow_site_df * 1/pseudovar * c2)

      for(j in 1:D){
        kappa[j] <- 1/(1 + nrow_site_df * 1/pseudovar * tau2 * 1 * lambda2[j])
        kappat[j] <- (1-b)*kappa[j] + b
      }
      meff <- sum(1 - kappa[1:D])
      mefft <- (1-b)*meff
    }

    # initialize
    inits = function(){
      list(
        # coefficients
        beta0 = 0,
        beta = rep(0, D),
        alpha0 = 0,
        alpha = rep(0, ncol(w)),
        # horseshoe prior inits
        lambda = rep(1, D),
        tau = .05,
        c2.rep = 1,
        z = rep(1, nsite)
      )
    }

    data_list = list(
      nsite = nrow(site_df),
      D = ncol(x),
      x = x,
      w = w,
      ncol_W = ncol(w),
      nrow_site_df = nrow(site_df),
      nrow_samp_df = nrow(samp_df),
      site_season_ndx = samp_df$site,
      y = samp_df$y,
      df.slab = priors$c$df_slab,
      scale.slab = priors$c$scale_slab,
      mean_m1 = mixture_mean1,
      mean_m2 = mixture_mean2,
      var_m1 = mixture_var1,
      var_m2 = mixture_var2,
      weight_m1 = mixture_weight1,
      weight_m2 = mixture_weight2,
      beta0mu = priors$beta0$mu0,
      beta0var = (priors$beta0$sigma0)^2,
      alpha0mu = priors$alpha0$mu0,
      alpha0var = (priors$alpha0$sigma0)^2,
      alphamu = priors$alpha$mu0,
      alphavar = (priors$alpha$sigma0)^2,
      sqrtpseudovar = sqrt(pseudovar),
      pseudovar = pseudovar,
      sqrtnsite = sqrt(nrow(site_df))
    )

    params <- c("alpha0", "alpha", "beta0", "beta", "tau", "lambdat", "kappat", "mefft")

  }



  # y
  if(sum(stringr::str_detect(colnames(samp_df), "y")) != 1){
    stop("ensure you have one column in samp_df called `y`")
  }



  # priors
  if(is.null(mixture_mean1)){
    stop("You must specify a value for mixture_mean1.")
  }
  if(is.null(mixture_mean2)){
    stop("You must specify a value for mixture_mean2.")
  }
  if(mixture_mean1 < 1 | mixture_mean1 > (D-1)){
    stop("You must specify a value between 1 and D-1 for mixture_mean1.")
  }
  if(mixture_mean2 < 1 | mixture_mean2 > (D-1)){
    stop("You must specify a value between 1 and D-1 for mixture_mean2.")
  }



  fit <- R2jags::jags.parallel(
    model.file = jags_reghs_occ,
    n.chains = nchain,
    inits = inits,
    data = data_list,
    n.iter = niter,
    n.burnin = nburnin,
    n.thin = thin,
    parameters.to.save = params,
    jags.seed = 1,
    jags.module = c("glm")
  )

  out <- convert_to_list(fit)
  class(out) <- c('horseshoeocc_mcmc', class(out))

  out <- list(mcmc = out,
              data = list(x = x,
                          w = w,
                          site_df = site_df,
                          samp_df = samp_df),
              model = list(det_model = det_model,
                           occ_model = occ_model),
              opts = list(niter = niter,
                          nchain = nchain,
                          nburnin = nburnin,
                          thin = thin)
  )

  class(out) <- c('horseshoeocc', class(out))
  return(out)
}



