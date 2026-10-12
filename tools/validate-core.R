# Development check before installing the package. Run at the repository root.
Sys.setenv(PKG_CPPFLAGS = paste0("-I", shQuote(normalizePath("inst/include"))))
Rcpp::sourceCpp("src/gaussian.cpp", cacheDir = "/tmp/censoredGP-source-cache")
set.seed(914)
worst <- 0
for (kernel in 1:4) {
  for (n_observed in c(0L, 2L)) {
    n_censored <- 5L
    n <- n_observed + n_censored
    coordinates <- matrix(runif(n * 2), ncol = 2)
    distance2 <- as.matrix(dist(coordinates))^2
    group <- rep(1:3, length.out = n)
    group_cov <- ifelse(outer(group, group, "=="), 1, -0.5)
    theta <- c(1.2, 0.6, 0.7, 0.4, 0.25)
    observed <- rnorm(n_observed)
    lower <- c(-Inf, -0.2, -0.5, -Inf, 0.1)
    upper <- c(0.3, Inf, 0.8, 0.1, Inf)
    unit <- seq(0.2, 0.8, length.out = n_censored)
    state <- c(theta, 0.2, unit)
    evaluate <- function(state) {
      cgp_censor_evaluate(
        distance2, group_cov, state[1:5], state[6],
        observed, state[-(1:6)], lower, upper, kernel
      )
    }
    result <- evaluate(state)
    analytical <- c(result$gradient_theta, result$gradient_mean, result$gradient_unit)
    numerical <- numDeriv::grad(function(state) evaluate(state)$log_density, state)
    error <- max(abs(analytical - numerical) / pmax(1, abs(numerical)))
    stopifnot(error < 2e-6)
    worst <- max(worst, error)
    jacobian <- numDeriv::jacobian(function(u) evaluate(c(theta, 0.2, u))$censored_values, unit)
    gaussian <- cgp_evaluate(
      distance2, group_cov, theta,
      c(observed, result$censored_values) - 0.2, kernel, FALSE
    )
    error_identity <- abs(gaussian$log_density + sum(log(abs(diag(jacobian)))) - result$log_density)
    stopifnot(error_identity < 1e-7)
  }
}
for (bounds in list(c(-Inf, -40), c(40, Inf), c(39.9, 40), c(-40, -39.9))) {
  result <- cgp_censor_evaluate(
    matrix(0, 1, 1), matrix(0, 1, 1),
    c(1, 1, 1, 1, 0.01), 0, numeric(), 0.4, bounds[1], bounds[2], 1L
  )
  stopifnot(
    is.finite(result$log_density), all(is.finite(result$gradient_theta)),
    is.finite(result$gradient_mean), result$censored_values > bounds[1],
    result$censored_values < bounds[2]
  )
}
cat("Exact quantile core: PASS. Maximum relative gradient error:", worst, "\n")
