# Keep these identifiers synchronized with the C++ stationary kernels.
kernel_id <- function(kernel) {
  match(
    match.arg(kernel, c("exp_quad", "exponential", "matern32", "matern52")),
    c("exp_quad", "exponential", "matern32", "matern52")
  )
}

check_count <- function(value, name, minimum = 1L) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
    value != floor(value) || value < minimum || value > .Machine$integer.max) {
    stop(sprintf("`%s` must be a whole number of at least %s.", name, minimum), call. = FALSE)
  }
  as.integer(value)
}

with_local_seed <- function(seed, expression) {
  if (is.null(seed)) {
    return(force(expression))
  }
  seed <- check_count(seed, "seed", 0L)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(
    {
      if (had_seed) {
        assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  set.seed(seed)
  force(expression)
}

formula_parts <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("Use a formula such as `cens(value, type) ~ time | subject`.", call. = FALSE)
  }
  rhs <- formula[[3L]]
  while (is.call(rhs) && identical(rhs[[1L]], as.name("("))) rhs <- rhs[[2L]]
  grouped <- is.call(rhs) && identical(rhs[[1L]], as.name("|"))
  list(
    response = formula[[2L]], coordinates = if (grouped) rhs[[2L]] else rhs,
    group = if (grouped) rhs[[3L]] else NULL
  )
}

evaluate_coordinates <- function(expression, data, environment, dimension = NULL) {
  coordinates <- eval(expression, data, environment)
  if (!is.numeric(coordinates) || !length(coordinates) || any(!is.finite(coordinates))) {
    stop("GP inputs must be finite numeric values or a numeric `cbind(...)` matrix.", call. = FALSE)
  }
  if (is.null(dim(coordinates))) coordinates <- matrix(coordinates, ncol = 1L)
  if (!is.matrix(coordinates) || (length(dimension) && ncol(coordinates) != dimension)) {
    stop("New GP inputs must have the same number of columns as the training inputs.", call. = FALSE)
  }
  coordinates
}

squared_distance <- function(first, second = first) {
  # Subtract coordinates before squaring. The norm/cross-product identity loses
  # short distances when coordinates have a large common offset (e.g. dates).
  result <- matrix(0, nrow(first), nrow(second))
  for (column in seq_len(ncol(first))) {
    result <- result + outer(first[, column], second[, column], "-")^2
  }
  if (any(!is.finite(result))) {
    stop("Squared GP distances overflow; express inputs in smaller units.", call. = FALSE)
  }
  result
}

group_covariance <- function(first, second, n_groups, structure, active = TRUE) {
  if (!active) {
    return(matrix(0, length(first), length(second)))
  }
  same_group <- outer(first, second, "==")
  if (structure == "independent") {
    return(1 * same_group)
  }
  ifelse(same_group, 1, -1 / (n_groups - 1))
}

prepare_model <- function(formula, data, priors, kernel, mean, group_structure,
                          parameterization) {
  if (!inherits(priors, "cgp_priors")) stop("`priors` must be made with `cgp_priors()`.", call. = FALSE)
  if (!is.data.frame(data) && !is.list(data)) stop("Supply a data frame or named list as `data`.", call. = FALSE)
  parts <- formula_parts(formula)
  environment <- new.env(parent = environment(formula))
  # A formula can be used with censoredGP::cgp() without attaching the package.
  environment$cens <- cens
  coordinates <- evaluate_coordinates(parts$coordinates, data, environment)
  response <- eval(parts$response, data, environment)
  if (is.numeric(response) && !inherits(response, "cgp_response") && is.null(dim(response))) {
    response <- cens(response)
  }
  if (!inherits(response, "cgp_response") || nrow(response) != nrow(coordinates)) {
    stop("The formula response must be numeric or `cens(...)`, with one row per GP input.", call. = FALSE)
  }
  bounds <- as_censoring_data(response)
  group <- if (is.null(parts$group)) {
    rep(".single", nrow(coordinates))
  } else {
    eval(parts$group, data, environment)
  }
  if (!is.atomic(group) || length(group) != nrow(coordinates) || anyNA(group)) {
    stop("The grouping expression must give one non-missing label per observation.", call. = FALSE)
  }
  group <- as.character(group)
  group_levels <- unique(group)
  group_index <- match(group, group_levels)
  has_group <- !is.null(parts$group) &&
    (length(group_levels) > 1L || group_structure == "independent")
  estimate_mean <- identical(mean, "constant")
  if (!estimate_mean && (!is.numeric(mean) || length(mean) != 1L || !is.finite(mean))) {
    stop("`mean` must be 'constant' or a fixed numeric mean, such as 0.", call. = FALSE)
  }
  # Observed rows precede censored rows for the exact triangular quantile map.
  # Keep an inverse permutation so imputations return in the user's row order.
  order <- order(bounds$type != "none", seq_len(nrow(bounds)))
  n_observed <- sum(bounds$type == "none")
  censor_rows <- order[bounds$type[order] != "none"]
  censor_types <- bounds$type[censor_rows]
  geometry <- squared_distance(coordinates[order, , drop = FALSE])
  group_cov <- group_covariance(
    group_index[order], group_index[order],
    length(group_levels), group_structure, has_group
  )
  stan_data <- c(
    list(
      N = nrow(bounds), n_observed = n_observed,
      has_group = as.integer(has_group), estimate_mean = as.integer(estimate_mean),
      fixed_mean = if (estimate_mean) 0 else mean, kernel = kernel_id(kernel),
      quantile = as.integer(parameterization == "quantile"), distance2 = geometry,
      group_cov = group_cov, observed = bounds$lower[bounds$type == "none"],
      censor_kind = match(censor_types, c("left", "right")),
      # Each censored measurement has one finite detection limit. Stan supplies
      # the infinite endpoint from the censoring direction.
      censor_limit = ifelse(censor_types == "left", bounds$upper[censor_rows],
        bounds$lower[censor_rows]
      )
    ),
    encode_priors(priors)
  )
  list(
    stan_data = stan_data, coordinates = coordinates, bounds = bounds,
    group = group, group_levels = group_levels, group_index = group_index,
    group_structure = group_structure, has_group = as.integer(has_group),
    estimate_mean = as.integer(estimate_mean), fixed_mean = stan_data$fixed_mean,
    order = order, censor_rows = censor_rows, formula_parts = parts,
    coordinate_label = paste(deparse(parts$coordinates), collapse = "")
  )
}

initial_values <- function(prepared, priors, parameterization) {
  names <- c(
    "amplitude", "lengthscale",
    if (prepared$has_group) c("group_amplitude", "group_lengthscale"), "sigma"
  )
  covariance_parameters <- vapply(priors[names], prior_quantile, numeric(1), probability = 0.5)
  mean_parameter <- if (prepared$estimate_mean) {
    prior_quantile(priors$intercept, 0.5, positive = FALSE)
  } else {
    numeric()
  }
  bounds <- prepared$bounds[prepared$censor_rows, , drop = FALSE]
  values <- numeric(nrow(bounds))
  values[bounds$type == "left"] <- bounds$upper[bounds$type == "left"] - utils::tail(covariance_parameters, 1L)
  values[bounds$type == "right"] <- bounds$lower[bounds$type == "right"] + utils::tail(covariance_parameters, 1L)
  list(
    covariance_parameters = unname(covariance_parameters), mean_parameter = mean_parameter,
    censored = if (parameterization == "bounded") values else numeric(),
    unit = if (parameterization == "quantile") rep(0.5, nrow(bounds)) else numeric()
  )
}
