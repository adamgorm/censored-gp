#' Fit a Full Bayesian Censored Gaussian Process
#'
#' The formula `cens(value, type) ~ time | subject` fits a population
#' curve and subject deviations. Omit `| subject` for a single curve. A numeric
#' response without [cens()] denotes entirely uncensored observations.
#' Inputs may be a numeric vector or `cbind(x1, x2, ...)` for an isotropic GP.
#'
#' The latent function is integrated out analytically. Only censored noisy
#' measurements are augmented. All covariance hyperparameters are sampled
#' jointly using Stan NUTS; no plug-in estimates or covariance approximations
#' are used. Both available parameterizations have the same posterior target.
#'
#' @param formula Model formula with per-observation censoring on the left.
#' @param data Data frame or named list containing formula variables.
#' @param priors Named prior bundle created with [cgp_priors()].
#' @param kernel One of `"exp_quad"`, `"exponential"`, `"matern32"`, or
#'   `"matern52"`. Length scales use the original coordinate units.
#' @param mean `"constant"` estimates an intercept with its own prior. A finite
#'   number fixes the prior mean of the population GP; use `0` for the paper's
#'   zero-mean specification.
#' @param group_structure `"sum_to_zero"` preserves the paper's constraint that
#'   group deviations sum to zero, with covariance -1/(G-1) between groups.
#'   `"independent"` uses independent deviations and supports new-group prediction.
#' @param parameterization `"quantile"` uses an exact triangular reparameterization;
#'   `"bounded"` directly samples measurements within their censoring limits.
#' @param chains,parallel_chains Numbers of MCMC chains and concurrent chains.
#' @param iter_warmup,iter_sampling Warmup and retained iterations per chain.
#' @param seed Random seed passed to Stan.
#' @param adapt_delta Target HMC acceptance probability.
#' @param max_treedepth Maximum NUTS tree depth.
#' @param refresh Stan progress frequency; zero suppresses progress.
#' @param diagnostic_warnings Warn about observed sampling or convergence issues.
#' @param ... Additional arguments to `cmdstanr::CmdStanModel$sample()`, including
#'   `init` or `output_dir`.
#' @return A serializable `cgp_fit` with posterior draws, priors, original data,
#'   and sampler diagnostics. Prediction does not depend on external CSV files.
#' @examples
#' \dontrun{
#' d <- data.frame(
#'   time = c(0, 1, 2, 3), value = c(1, 0.5, 2, 1.5),
#'   type = c("none", "left", "right", "left")
#' )
#' fit <- cgp(cens(value, type) ~ time, d,
#'   priors = cgp_priors(lengthscale = lognormal_prior(log(2), 0.7))
#' )
#' predict(fit, data.frame(time = seq(0, 3, length.out = 50)))
#' }
#' @export
cgp <- function(formula, data, priors = cgp_priors(), kernel = "exp_quad",
                mean = "constant", group_structure = c("sum_to_zero", "independent"),
                parameterization = c("quantile", "bounded"), chains = 4L,
                parallel_chains = chains, iter_warmup = 1000L, iter_sampling = 1000L,
                seed = NULL, adapt_delta = 0.9, max_treedepth = 12L, refresh = 100L,
                diagnostic_warnings = TRUE, ...) {
  call <- match.call()
  group_structure <- match.arg(group_structure)
  parameterization <- match.arg(parameterization)
  kernel <- match.arg(kernel, c("exp_quad", "exponential", "matern32", "matern52"))
  chains <- check_count(chains, "chains")
  parallel_chains <- check_count(parallel_chains, "parallel_chains")
  iter_warmup <- check_count(iter_warmup, "iter_warmup", 0L)
  iter_sampling <- check_count(iter_sampling, "iter_sampling")
  max_treedepth <- check_count(max_treedepth, "max_treedepth")
  if (!is.numeric(adapt_delta) || length(adapt_delta) != 1L || !is.finite(adapt_delta) ||
    adapt_delta <= 0 || adapt_delta >= 1) {
    stop("`adapt_delta` must be strictly between zero and one.", call. = FALSE)
  }
  prepared <- prepare_model(formula, data, priors, kernel, mean, group_structure, parameterization)
  model <- compile_cgp_model()
  additional <- list(...)
  # Report diagnostics once, through our interface, after retaining all draws.
  # CmdStanR would otherwise issue its own warnings before diagnostic_warnings
  # can take effect, and then repeat them in diagnostic_summary().
  if (!"diagnostics" %in% names(additional)) additional["diagnostics"] <- list(NULL)
  # Preserve constraint endpoints when posterior imputations are written to CSV.
  if (is.null(additional$sig_figs)) additional$sig_figs <- 18L
  if (is.null(additional$init)) additional$init <- initial_values(prepared, priors, parameterization)
  # CmdStanR expects either a function or one initialization list per chain.
  if (is.list(additional$init) && !is.null(names(additional$init))) {
    additional$init <- rep(list(additional$init), chains)
  }
  elapsed <- system.time({
    stan_fit <- do.call(model$sample, c(list(
      data = prepared$stan_data,
      chains = chains, parallel_chains = min(chains, parallel_chains),
      iter_warmup = iter_warmup, iter_sampling = iter_sampling, seed = seed,
      adapt_delta = adapt_delta, max_treedepth = max_treedepth, refresh = refresh
    ), additional))
  })[["elapsed"]]
  covariance_names <- c(
    "amplitude", "lengthscale",
    if (prepared$has_group) c("group_amplitude", "group_lengthscale"), "sigma"
  )
  raw_names <- paste0("covariance_parameters[", seq_along(covariance_names), "]")
  variables <- c(raw_names, "intercept", if (length(prepared$censor_rows)) "y_cens")
  draws <- stan_fit$draws(variables)
  labels <- dimnames(draws)[[3L]]
  labels[match(raw_names, labels)] <- covariance_names
  dimnames(draws)[[3L]] <- labels
  diagnostics <- stan_fit$diagnostic_summary(quiet = TRUE)
  sampler <- stan_fit$sampler_diagnostics(format = "array")
  # BFMI uses adjacent energy differences within each chain, never across chains.
  energy <- sampler[, , "energy__", drop = FALSE]
  diagnostics$bfmi <- vapply(seq_len(chains), function(chain) {
    values <- as.numeric(energy[, chain, 1L])
    mean(diff(values)^2) / stats::var(values)
  }, numeric(1))
  object <- structure(
    list(
      call = call, formula = formula, data = as.data.frame(data),
      prepared = prepared, priors = priors, kernel = kernel, parameterization = parameterization,
      draws = draws, diagnostics = diagnostics, elapsed = elapsed,
      chains = chains, iterations = iter_sampling, max_treedepth = max_treedepth,
      cmdstan_version = as.character(cmdstanr::cmdstan_version()),
      stan_code = cgp_stancode(),
      cpp_md5 = tools::md5sum(list.files(
        system.file("include", "censoredGP", package = "censoredGP"),
        pattern = "\\.hpp$", full.names = TRUE
      )),
      parameter_names = c(covariance_names, if (prepared$estimate_mean) "intercept")
    ),
    class = "cgp_fit"
  )
  if (diagnostic_warnings) warn_diagnostics(object)
  object
}

warn_diagnostics <- function(object) {
  diagnostics <- object$diagnostics
  # Check augmented measurements as well as hyperparameters: a well-mixed
  # length scale does not establish convergence of every censored value.
  variables <- c(
    object$parameter_names,
    if (length(object$prepared$censor_rows)) "y_cens"
  )
  summary <- posterior::summarise_draws(as_draws.cgp_fit(object, variables))
  issues <- character()
  if (any(diagnostics$num_divergent > 0)) issues <- c(issues, "divergent transitions")
  if (any(diagnostics$num_max_treedepth > 0)) issues <- c(issues, "maximum tree depth reached")
  if (any(diagnostics$bfmi < 0.3, na.rm = TRUE)) issues <- c(issues, "low BFMI")
  if (object$chains > 1L && any(summary$rhat > 1.01, na.rm = TRUE)) {
    issues <- c(issues, "R-hat above 1.01")
  }
  if (length(issues)) {
    warning("Sampling diagnostics: ", paste(issues, collapse = "; "),
      ". Inspect summary(fit) and diagnostics(fit).",
      call. = FALSE
    )
  }
}

#' Extract Posterior Draws
#' @param x A fitted [cgp()] model.
#' @param variables Optional posterior variable names. Defaults to the estimated
#'   hyperparameters. Use `"y_cens"` for augmented measurements.
#' @param ... Reserved for compatibility with `posterior::as_draws()`.
#' @return A `posterior::draws_array`, preserving iteration and chain dimensions.
#' @export
as_draws.cgp_fit <- function(x, variables = NULL, ...) {
  if (is.null(variables)) variables <- x$parameter_names
  posterior::subset_draws(x$draws, variable = variables)
}

#' Inspect MCMC Diagnostics
#' @param object A fitted [cgp()] model.
#' @return One row per chain with divergent transitions, tree-depth hits, and BFMI.
#' @export
diagnostics <- function(object) {
  if (!inherits(object, "cgp_fit")) stop("Supply a cgp fit.", call. = FALSE)
  values <- object$diagnostics
  data.frame(
    chain = seq_along(values$num_divergent),
    num_divergent = values$num_divergent,
    num_max_treedepth = values$num_max_treedepth, bfmi = values$bfmi
  )
}

#' @export
summary.cgp_fit <- function(object, ...) {
  structure(
    list(
      parameters = posterior::summarise_draws(as_draws.cgp_fit(object)),
      diagnostics = diagnostics(object), elapsed = object$elapsed,
      n = nrow(object$prepared$bounds), censoring = table(object$prepared$bounds$type),
      kernel = object$kernel, chains = object$chains, iterations = object$iterations
    ),
    class = "summary.cgp_fit"
  )
}

#' @export
print.summary.cgp_fit <- function(x, ...) {
  cat("Full Bayesian censored Gaussian process (", x$kernel, ")\n", sep = "")
  cat(x$n, "observations;", x$chains, "chains with", x$iterations, "draws each\n")
  print(x$censoring)
  print(as.data.frame(x$parameters), row.names = FALSE, digits = 3)
  cat("\nSampler diagnostics:\n")
  print(x$diagnostics, row.names = FALSE, digits = 3)
  cat("Sampling elapsed:", round(x$elapsed, 1), "seconds\n")
  invisible(x)
}

#' @export
print.cgp_fit <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}
