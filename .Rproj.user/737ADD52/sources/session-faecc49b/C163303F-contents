#' Internal function: logit transformation
#'
#' @noRd
#'
logit <- function(p){
  log(p/(1 - p))
}

#' Internal function: inverse logit transformation
#'
#' @noRd
#'
invlogit <- function(x){
  exp(x)/(1 + exp(x))
}

#' Internal function: convert fit object to list
#'
#' @noRd
#'
convert_to_list <- function(fit){
  nchains <- fit$BUGSoutput$n.chains
  nkeep <- fit$BUGSoutput$n.keep
  array <- fit$BUGSoutput$sims.array

  out <- list()
  for(i in 1:(dim(array)[2])) out[[i]] <- array[,i,]
  return(out)
}
