select_draws <- function(object, ndraws = NULL, draw_ids = NULL) {
  # posterior's matrix class deliberately keeps a matrix when selecting one row.
  # Drop that class so downstream single-draw extraction produces a named vector.
  draws <- unclass(posterior::as_draws_matrix(object$draws))
  total <- nrow(draws)
  if (!is.null(ndraws) && !is.null(draw_ids)) {
    stop("Specify either `ndraws` or `draw_ids`.", call. = FALSE)
  }
  if (!is.null(ndraws)) {
    ndraws <- check_count(ndraws, "ndraws")
    if (ndraws > total) stop("`ndraws` exceeds the number of posterior draws.", call. = FALSE)
    # Deterministic thinning spans all chains and never favors the first chain.
    draw_ids <- unique(as.integer(seq(1, total, length.out = ndraws)))
  }
  if (is.null(draw_ids)) draw_ids <- seq_len(total)
  if (!is.numeric(draw_ids) || !length(draw_ids) || any(!is.finite(draw_ids)) ||
    any(draw_ids != floor(draw_ids)) || any(draw_ids < 1 | draw_ids > total)) {
    stop("`draw_ids` must identify existing posterior draws.", call. = FALSE)
  }
  draws[draw_ids, , drop = FALSE]
}

stationary_kernel <- function(distance2, amplitude, lengthscale, kernel) {
  radius <- sqrt(distance2) / lengthscale
  correlation <- switch(kernel,
    exp_quad = exp(-radius^2 / 2),
    exponential = exp(-radius),
    matern32 = (1 + sqrt(3) * radius) * exp(-sqrt(3) * radius),
    matern52 = (1 + sqrt(5) * radius + 5 * radius^2 / 3) * exp(-sqrt(5) * radius)
  )
  amplitude^2 * correlation
}

prediction_design <- function(object, newdata, component, joint = TRUE) {
  prepared <- object$prepared
  if (is.null(newdata)) newdata <- object$data
  environment <- environment(object$formula)
  coordinates <- evaluate_coordinates(
    prepared$formula_parts$coordinates, newdata,
    environment, ncol(prepared$coordinates)
  )
  n <- nrow(coordinates)
  if (length(component) == 1L) component <- rep(component, n)
  if (length(component) != n || anyNA(component) ||
    any(!component %in% c("response", "population", "deviation"))) {
    stop("`component` must give response, population, or deviation per prediction row.", call. = FALSE)
  }
  group <- rep(".population", n)
  need_group <- component != "population"
  if (!is.null(prepared$formula_parts$group) && any(need_group)) {
    group <- tryCatch(eval(prepared$formula_parts$group, newdata, environment),
      error = function(error) stop("New data must contain the model's grouping variable.", call. = FALSE)
    )
    if (length(group) != n || anyNA(group[need_group])) {
      stop("Supply one group label per response or deviation prediction.", call. = FALSE)
    }
    group <- as.character(group)
    group[!need_group] <- ".population"
    if (prepared$group_structure == "sum_to_zero" &&
      any(!group[need_group] %in% prepared$group_levels)) {
      stop("New groups require a model fitted with group_structure = 'independent'.",
        call. = FALSE
      )
    }
  }
  ordered <- prepared$order
  list(
    coordinates = coordinates, group = group, component = component,
    distance_cross = squared_distance(coordinates, prepared$coordinates[ordered, , drop = FALSE]),
    # Pointwise summaries need only prediction-to-training distances. Avoid
    # allocating an M by M geometry when no joint covariance is requested.
    distance_prior = if (joint) squared_distance(coordinates) else NULL,
    group_cross = group_covariance(
      group, prepared$group[ordered],
      length(prepared$group_levels), prepared$group_structure, prepared$has_group
    ),
    group_prior = if (joint) {
      group_covariance(
        group, group, length(prepared$group_levels),
        prepared$group_structure, prepared$has_group
      )
    } else {
      NULL
    },
    population_weight = as.numeric(component != "deviation"),
    group_weight = as.numeric(component != "population" & prepared$has_group == 1L)
  )
}

# Gaussian conditioning may be exactly singular for duplicate predictions or
# sum-to-zero components. Use an eigendecomposition only when Cholesky fails;
# retain every positive eigenvalue and clip only negative floating-point error.
sample_gaussian <- function(mean, covariance) {
  covariance <- (covariance + t(covariance)) / 2
  factor <- tryCatch(chol(covariance), error = function(error) NULL)
  noise <- stats::rnorm(length(mean))
  if (!is.null(factor)) {
    return(as.vector(mean + crossprod(factor, noise)))
  }
  decomposition <- eigen(covariance, symmetric = TRUE)
  tolerance <- 100 * .Machine$double.eps * nrow(covariance) *
    max(1, max(abs(decomposition$values)))
  if (min(decomposition$values) < -tolerance) {
    stop("Prediction covariance is not positive semidefinite; inspect model conditioning.", call. = FALSE)
  }
  as.vector(mean + decomposition$vectors %*% (sqrt(pmax(decomposition$values, 0)) * noise))
}

posterior_function_draws <- function(object, newdata, component, ndraws, draw_ids,
                                     seed, joint, noisy) {
  if (!inherits(object, "cgp_fit")) stop("Supply a fitted cgp model.", call. = FALSE)
  draws <- select_draws(object, ndraws, draw_ids)
  design <- prediction_design(object, newdata, component, joint = joint)
  n_predictions <- nrow(design$coordinates)
  population_active <- any(design$population_weight != 0)
  group_active <- any(design$group_weight != 0)
  n_censored <- length(object$prepared$censor_rows)
  n_observed <- object$prepared$stan_data$n_observed
  censored_names <- paste0("y_cens[", seq_len(n_censored), "]")
  kernel <- object$kernel
  with_local_seed(seed, {
    result <- matrix(NA_real_, nrow(draws), n_predictions)
    for (index in seq_len(nrow(draws))) {
      draw <- draws[index, ]
      theta <- c(
        draw["amplitude"], draw["lengthscale"],
        if (object$prepared$has_group) draw[c("group_amplitude", "group_lengthscale")] else c(1, 1),
        draw["sigma"]
      )
      augmented <- c(
        object$prepared$stan_data$observed,
        if (n_censored) draw[censored_names] else numeric()
      )
      factor <- cgp_evaluate(
        object$prepared$stan_data$distance2,
        object$prepared$stan_data$group_cov, theta, augmented - draw["intercept"],
        kernel_id(kernel), FALSE
      )$factor
      # Do not evaluate kernels for components absent from the request.
      cross <- matrix(0, n_predictions, length(augmented))
      if (population_active) {
        cross <- cross + design$population_weight *
          stationary_kernel(design$distance_cross, theta[1], theta[2], kernel)
      }
      if (group_active) {
        cross <- cross + design$group_weight * design$group_cross *
          stationary_kernel(design$distance_cross, theta[3], theta[4], kernel)
      }
      prior_variance <- design$population_weight * theta[1]^2 +
        design$group_weight * theta[3]^2
      if (joint) {
        prior <- matrix(0, n_predictions, n_predictions)
        if (population_active) {
          prior <- prior +
            outer(design$population_weight, design$population_weight) *
              stationary_kernel(design$distance_prior, theta[1], theta[2], kernel)
        }
        if (group_active) {
          prior <- prior +
            outer(design$group_weight, design$group_weight) * design$group_prior *
              stationary_kernel(design$distance_prior, theta[3], theta[4], kernel)
        }
        conditional <- cgp_condition(factor, cross, prior_variance, prior, TRUE)
        conditional$mean <- conditional$mean + design$population_weight * draw["intercept"]
        if (noisy) diag(conditional$covariance) <- diag(conditional$covariance) + theta[5]^2
        result[index, ] <- sample_gaussian(conditional$mean, conditional$covariance)
      } else {
        conditional <- cgp_condition(factor, cross, prior_variance, matrix(numeric(), 0, 0), FALSE)
        variance <- as.numeric(conditional$variance)
        tolerance <- 100 * .Machine$double.eps * max(1, max(prior_variance)) * length(augmented)
        if (min(variance) < -tolerance) stop("Negative prediction variance.", call. = FALSE)
        variance <- pmax(variance, 0) + if (noisy) theta[5]^2 else 0
        result[index, ] <- conditional$mean + design$population_weight * draw["intercept"] +
          stats::rnorm(n_predictions, sd = sqrt(variance))
      }
    }
    colnames(result) <- paste0("prediction[", seq_len(n_predictions), "]")
    result
  })
}

#' Draw Joint Posterior Smooth Curves
#'
#' Draws include GP function uncertainty, covariance parameter uncertainty, and
#' uncertainty about censored measurements. Dependence between every requested
#' prediction is preserved within a row. Requests for different components can
#' be combined in one call using a component vector, one label per new-data row.
#'
#' @param object A fitted [cgp()] model.
#' @param newdata Prediction data; defaults to the training data. Population
#'   predictions do not require a grouping column.
#' @param component `"response"` for population plus group deviation,
#'   `"population"` for the common curve, or `"deviation"` for group departures.
#'   May be a vector specifying the component separately for each new-data row.
#' @param ndraws Number of evenly spaced posterior draws to use; defaults to all.
#' @param draw_ids Explicit posterior draw indices, as an alternative to `ndraws`.
#' @param seed Seed for conditional function simulation; preserves R's RNG state.
#' @return Numeric matrix with posterior draws in rows and prediction locations
#'   in columns. [posterior_predict()] additionally includes measurement noise.
#' @export
posterior_epred <- function(object, newdata = NULL, component = "response",
                            ndraws = NULL, draw_ids = NULL, seed = NULL) {
  posterior_function_draws(object, newdata, component, ndraws, draw_ids, seed,
    joint = TRUE, noisy = FALSE
  )
}

#' Draw Joint Posterior Predictive Measurements
#' @inheritParams posterior_epred
#' @return A draws-by-location matrix of future noisy measurements, without
#'   censoring them at detection limits. These differ from imputations of the
#'   original censored measurements; use [impute_censored()] for those.
#' @export
posterior_predict <- function(object, newdata = NULL, ndraws = NULL,
                              draw_ids = NULL, seed = NULL) {
  posterior_function_draws(object, newdata, "response", ndraws, draw_ids, seed,
    joint = TRUE, noisy = TRUE
  )
}

#' Summarize Posterior Predictions
#'
#' For efficiency this method simulates marginal predictions without building
#' their full covariance. The returned intervals are pointwise credible or
#' predictive intervals. Use [posterior_epred()] for coherent joint curves.
#'
#' @param object A fitted [cgp()] model.
#' @param newdata Prediction data; defaults to training locations.
#' @param type `"latent"` for smooth population plus group curves,
#'   `"response"` for future noisy measurements, `"population"`, or `"deviation"`.
#' @param level Probability covered by the pointwise interval.
#' @param ndraws,draw_ids,seed See [posterior_epred()].
#' @param ... Reserved for S3 compatibility.
#' @return A data frame with `mean`, `sd`, `lower`, and `upper` per location.
#' @export
predict.cgp_fit <- function(object, newdata = NULL,
                            type = c("latent", "response", "population", "deviation"),
                            level = 0.95, ndraws = NULL, draw_ids = NULL, seed = NULL, ...) {
  type <- match.arg(type)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be strictly between zero and one.", call. = FALSE)
  }
  component <- if (type %in% c("latent", "response")) "response" else type
  draws <- posterior_function_draws(object, newdata, component, ndraws, draw_ids, seed,
    joint = FALSE, noisy = type == "response"
  )
  probabilities <- c((1 - level) / 2, (1 + level) / 2)
  intervals <- apply(draws, 2L, stats::quantile, probs = probabilities, names = FALSE)
  data.frame(
    mean = colMeans(draws), sd = apply(draws, 2L, stats::sd),
    lower = intervals[1L, ], upper = intervals[2L, ]
  )
}

#' Impute the Original Censored Measurements
#' @inheritParams posterior_epred
#' @return A draws-by-observation matrix in the original data order. Exact
#'   observations remain fixed; censored measurements satisfy their own bounds.
#' @export
impute_censored <- function(object, ndraws = NULL, draw_ids = NULL) {
  draws <- select_draws(object, ndraws, draw_ids)
  result <- matrix(rep(object$prepared$bounds$lower, each = nrow(draws)), nrow = nrow(draws))
  rows <- object$prepared$censor_rows
  if (length(rows)) result[, rows] <- draws[, paste0("y_cens[", seq_along(rows), "]"), drop = FALSE]
  colnames(result) <- paste0("observation[", seq_len(ncol(result)), "]")
  result
}

#' @export
fitted.cgp_fit <- function(object, ...) {
  predict.cgp_fit(object, type = "latent", ...)
}
