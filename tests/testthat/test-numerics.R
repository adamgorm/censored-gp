test_that("Gaussian and exact quantile gradients match numerical derivatives", {
  skip_if_not_installed("numDeriv")
  set.seed(914)
  for (kernel in 1:4) {
    for (n_observed in c(0L, 2L)) {
      n_censored <- 3L
      n <- n_observed + n_censored
      coordinates <- matrix(runif(n * 2), ncol = 2)
      distance2 <- as.matrix(dist(coordinates))^2
      group <- rep(1:3, length.out = n)
      group_cov <- ifelse(outer(group, group, "=="), 1, -0.5)
      theta <- c(1.2, 0.6, 0.7, 0.4, 0.25)
      observed <- rnorm(n_observed)
      lower <- c(-Inf, -0.2, -0.5)
      upper <- c(0.3, Inf, 0.8)
      state <- c(theta, 0.2, 0.3, 0.5, 0.7)
      evaluate <- function(state) {
        cgp_censor_evaluate(
          distance2, group_cov, state[1:5], state[6],
          observed, state[-(1:6)], lower, upper, kernel
        )
      }
      result <- evaluate(state)
      analytical <- c(result$gradient_theta, result$gradient_mean, result$gradient_unit)
      numerical <- numDeriv::grad(function(state) evaluate(state)$log_density, state)
      expect_lt(max(abs(analytical - numerical) / pmax(1, abs(numerical))), 2e-6)

      # Check the density against an independently parameterized Gaussian plus
      # the numerical Jacobian, rather than only differentiating the same code.
      jacobian <- numDeriv::jacobian(
        function(u) evaluate(c(theta, 0.2, u))$censored_values,
        state[-(1:6)]
      )
      gaussian <- cgp_evaluate(
        distance2, group_cov, theta,
        c(observed, result$censored_values) - 0.2, kernel, TRUE
      )
      expect_equal(gaussian$log_density + sum(log(abs(diag(jacobian)))), result$log_density,
        tolerance = 1e-7
      )
      residual <- c(observed, result$censored_values) - 0.2
      gaussian_function <- function(v) {
        cgp_evaluate(
          distance2, group_cov,
          v[1:5], v[-(1:5)], kernel, FALSE
        )$log_density
      }
      score <- numDeriv::grad(gaussian_function, c(theta, residual))
      expect_equal(c(gaussian$gradient_theta, gaussian$gradient_y), as.numeric(score), tolerance = 1e-5)
    }
  }
})

test_that("extreme tails and narrow intervals have finite derivatives", {
  for (bounds in list(
    c(-Inf, -40), c(40, Inf), c(39.9, 40),
    c(-40, -39.9), c(-1e-7, 1e-7)
  )) {
    result <- cgp_censor_evaluate(
      matrix(0, 1, 1), matrix(0, 1, 1),
      c(1, 1, 1, 1, 0.01), 0, numeric(), 0.4, bounds[1], bounds[2], 1L
    )
    expect_true(is.finite(result$log_density))
    expect_true(all(is.finite(result$gradient_theta)))
    expect_true(is.finite(result$gradient_mean))
    expect_gt(as.numeric(result$censored_values), bounds[1])
    expect_lt(as.numeric(result$censored_values), bounds[2])
  }
})

test_that("input offsets do not erase short GP distances", {
  coordinates <- cbind(c(0, 1e9, 1e9 + 0.25), c(1e9, 1e9 + 0.5, 1e9 + 0.5))
  expect_equal(squared_distance(coordinates), unname(as.matrix(dist(coordinates)))^2,
    tolerance = 1e-15
  )
  expect_equal(squared_distance(coordinates)[2, 3], 0.25^2)
  expect_error(squared_distance(matrix(c(-1e308, 1e308), ncol = 1)), "overflow")
})
