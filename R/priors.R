#' Construct Explicit Hyperparameter Priors
#'
#' Priors are specified in the original response and input units. Normal and
#' Student-t priors on amplitudes, length scales, and noise SD are truncated at
#' zero. Their distribution parameters are fixed hyperparameters, not estimated
#' from the observed responses. Changing priors does not recompile Stan.
#'
#' @param mean,sd Normal location and standard deviation.
#' @return A validated `cgp_prior` object, accepted by [cgp_priors()].
#' @examples
#' normal_prior(0, 1)
#' lognormal_prior(log(60), 0.7)
#' @export
normal_prior <- function(mean = 0, sd = 1) {
  new_prior("normal", c(mean = unname(mean), sd = unname(sd)), positive = "sd")
}

#' @rdname normal_prior
#' @param df,location Student-t degrees of freedom and location.
#' @export
student_t_prior <- function(df = 4, location = 0, scale = 1) {
  new_prior("student_t", c(df = unname(df), location = unname(location), scale = unname(scale)),
    positive = c("df", "scale")
  )
}

#' @rdname normal_prior
#' @param meanlog,sdlog Mean and standard deviation of the natural logarithm.
#' @export
lognormal_prior <- function(meanlog = 0, sdlog = 1) {
  new_prior("lognormal", c(meanlog = unname(meanlog), sdlog = unname(sdlog)), positive = "sdlog")
}

#' @rdname normal_prior
#' @param shape,rate Gamma shape and rate (not scale).
#' @export
gamma_prior <- function(shape, rate) {
  new_prior("gamma", c(shape = unname(shape), rate = unname(rate)), positive = c("shape", "rate"))
}

#' @rdname normal_prior
#' @param scale Inverse-gamma scale, or Student-t scale.
#' @export
inv_gamma_prior <- function(shape, scale) {
  new_prior("inv_gamma", c(shape = unname(shape), scale = unname(scale)), positive = c("shape", "scale"))
}

#' @rdname normal_prior
#' @export
exponential_prior <- function(rate = 1) {
  new_prior("exponential", c(rate = unname(rate)), positive = "rate")
}

new_prior <- function(family, arguments, positive) {
  expected_names <- switch(family,
    normal = c("mean", "sd"),
    student_t = c("df", "location", "scale"),
    lognormal = c("meanlog", "sdlog"),
    gamma = c("shape", "rate"),
    inv_gamma = c("shape", "scale"),
    exponential = "rate"
  )
  if (!is.numeric(arguments) || any(!is.finite(arguments)) ||
    !identical(names(arguments), expected_names)) {
    stop("Prior arguments must be finite numeric scalars.", call. = FALSE)
  }
  if (any(arguments[positive] <= 0)) {
    stop("Prior shape, scale, rate, and dispersion arguments must be positive.", call. = FALSE)
  }
  structure(list(family = family, arguments = arguments), class = "cgp_prior")
}

#' @export
print.cgp_prior <- function(x, ...) {
  cat(format_prior(x), "\n")
  invisible(x)
}

format_prior <- function(prior) {
  arguments <- paste(format(prior$arguments, trim = TRUE, digits = 5), collapse = ", ")
  paste0(prior$family, "(", arguments, ")")
}

#' Set Priors for Named Model Components
#'
#' Every argument is a prior object such as [normal_prior()] or
#' [lognormal_prior()]. Unspecified arguments retain their documented defaults.
#' The defaults assume responses and inputs measured on approximately unit
#' scales; specify substantive priors for other units. No automatic response
#' centering or empirical Bayes estimation is performed.
#'
#' @param intercept Prior for the constant mean when `mean = "constant"`.
#'   Only normal or Student-t distributions are supported for this real parameter.
#' @param amplitude,lengthscale Population GP standard deviation and length scale.
#' @param group_amplitude,group_lengthscale Group deviation GP standard
#'   deviation and length scale. Used only in models with group deviations.
#' @param sigma Measurement noise standard deviation.
#' @return A `cgp_priors` object. [prior_summary()] displays its distributions.
#' @examples
#' priors <- cgp_priors(
#'   intercept = normal_prior(3, 1.5),
#'   amplitude = student_t_prior(4, 0, 2),
#'   lengthscale = lognormal_prior(log(60), 0.7),
#'   sigma = normal_prior(0, 0.5)
#' )
#' prior_summary(priors)
#' @export
cgp_priors <- function(intercept = normal_prior(0, 5),
                       amplitude = student_t_prior(4, 0, 2),
                       lengthscale = lognormal_prior(0, 1),
                       group_amplitude = student_t_prior(4, 0, 1),
                       group_lengthscale = lognormal_prior(0, 1),
                       sigma = student_t_prior(4, 0, 1)) {
  priors <- list(
    intercept = intercept, amplitude = amplitude, lengthscale = lengthscale,
    group_amplitude = group_amplitude, group_lengthscale = group_lengthscale,
    sigma = sigma
  )
  if (!all(vapply(priors, inherits, logical(1), "cgp_prior"))) {
    stop("Every component must use a prior constructor such as `normal_prior()`.", call. = FALSE)
  }
  if (!intercept$family %in% c("normal", "student_t")) {
    stop("The intercept requires a normal or Student-t prior on the real line.", call. = FALSE)
  }
  structure(priors, class = "cgp_priors")
}

#' Inspect Model Priors
#' @param object A [cgp_priors()] object or fitted [cgp()] model.
#' @return A data frame with parameter, distribution, support, units, and whether
#'   the prior is active in the supplied fitted model.
#' @export
prior_summary <- function(object) {
  fitted <- inherits(object, "cgp_fit")
  priors <- if (fitted) object$priors else object
  if (!inherits(priors, "cgp_priors")) stop("Supply model priors or a cgp fit.", call. = FALSE)
  parameters <- names(priors)
  active <- rep(TRUE, length(parameters))
  if (fitted) {
    active[parameters == "intercept"] <- object$prepared$estimate_mean == 1L
    active[grepl("^group_", parameters)] <- object$prepared$has_group == 1L
  }
  data.frame(
    parameter = parameters,
    distribution = vapply(priors, format_prior, character(1)),
    support = ifelse(parameters == "intercept", "real", "positive"),
    units = ifelse(grepl("lengthscale", parameters), "input units", "response units"),
    active = active, row.names = NULL
  )
}

#' @export
print.cgp_priors <- function(x, ...) {
  print(prior_summary(x), row.names = FALSE, ...)
  invisible(x)
}

# Match this family ordering to hyperprior_lpdf() in the Stan program.
encode_priors <- function(priors) {
  families <- c("normal", "student_t", "lognormal", "gamma", "inv_gamma", "exponential")
  arguments <- matrix(0, 6L, 3L)
  for (index in seq_along(priors)) {
    values <- priors[[index]]$arguments
    arguments[index, seq_along(values)] <- values
  }
  list(
    prior_family = match(vapply(priors, `[[`, character(1), "family"), families),
    prior_arguments = arguments
  )
}

prior_quantile <- function(prior, probability, positive = TRUE) {
  arguments <- unname(prior$arguments)
  if (prior$family == "normal") {
    if (!positive) {
      return(stats::qnorm(probability, arguments[1], arguments[2]))
    }
    # Work with log survival probabilities even for severely truncated priors.
    log_tail <- stats::pnorm(0, arguments[1], arguments[2], lower.tail = FALSE, log.p = TRUE)
    return(stats::qnorm(log_tail + log1p(-probability), arguments[1], arguments[2],
      lower.tail = FALSE, log.p = TRUE
    ))
  }
  if (prior$family == "student_t") {
    if (!positive) {
      return(arguments[2] + arguments[3] * stats::qt(probability, arguments[1]))
    }
    log_tail <- stats::pt(-arguments[2] / arguments[3], arguments[1],
      lower.tail = FALSE, log.p = TRUE
    )
    return(arguments[2] + arguments[3] * stats::qt(log_tail + log1p(-probability),
      arguments[1],
      lower.tail = FALSE, log.p = TRUE
    ))
  }
  switch(prior$family,
    lognormal = stats::qlnorm(probability, arguments[1], arguments[2]),
    gamma = stats::qgamma(probability, shape = arguments[1], rate = arguments[2]),
    inv_gamma = 1 / stats::qgamma(1 - probability, shape = arguments[1], rate = arguments[2]),
    exponential = stats::qexp(probability, rate = arguments[1])
  )
}

sample_prior <- function(prior, n, positive = TRUE) {
  prior_quantile(prior, stats::runif(n), positive)
}
