#' Full Bayesian Censored Gaussian Process Regression
#'
#' Models exact, left-censored, right-censored, and interval-censored functional
#' measurements using an exact dense Gaussian process covariance.
#'
#' @useDynLib censoredGP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom posterior as_draws
#' @keywords internal
"_PACKAGE"
