test_that("joint predictions reproduce an independent Gaussian conditional", {
  fit <- gaussian_fixture(n_draws = 4000L)
  locations <- c(0.2, 0.5, 0.9)
  prediction <- posterior_epred(fit, data.frame(time = locations), seed = 172)
  times <- fit$data$time
  covariance <- 1.2^2 * exp(-outer(times, times, "-")^2 / (2 * 0.5^2)) + diag(0.25^2, length(times))
  cross <- 1.2^2 * exp(-outer(locations, times, "-")^2 / (2 * 0.5^2))
  prior <- 1.2^2 * exp(-outer(locations, locations, "-")^2 / (2 * 0.5^2))
  expected_mean <- as.vector(0.2 + cross %*% solve(covariance, fit$data$y - 0.2))
  expected_covariance <- prior - cross %*% solve(covariance, t(cross))
  expect_equal(colMeans(prediction), expected_mean, tolerance = 0.03, ignore_attr = TRUE)
  expect_equal(unname(stats::cov(prediction)), expected_covariance, tolerance = 0.03)
  # Off-diagonal agreement checks that locations were not sampled independently.
  expect_gt(abs(stats::cov(prediction)[1, 2]), 0.001)
})

test_that("joint component predictions preserve the sum-to-zero constraint", {
  fit <- gaussian_fixture(grouped = TRUE, n_draws = 40L)
  newdata <- data.frame(time = 0.4, subject = c("A", "B", "C", "A", "A"))
  component <- c("deviation", "deviation", "deviation", "population", "response")
  draws <- posterior_epred(fit, newdata, component = component, seed = 46)
  expect_lt(max(abs(rowSums(draws[, 1:3]))), 1e-6)
  expect_lt(max(abs(draws[, 5] - draws[, 4] - draws[, 1])), 1e-6)
  expect_error(posterior_epred(fit, data.frame(time = 0.2, subject = "new")), "independent")
  independent <- gaussian_fixture(grouped = TRUE, structure = "independent", n_draws = 5L)
  expect_true(all(is.finite(posterior_epred(independent,
    data.frame(time = c(0.2, 0.7), subject = c("new", "new")),
    seed = 5
  ))))
})

test_that("saved fits predict without CmdStan files and preserve RNG state", {
  fit <- gaussian_fixture(n_draws = 20L)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(fit, path)
  restored <- readRDS(path)
  grid <- data.frame(time = c(0.1, 0.7))
  set.seed(57)
  before <- .Random.seed
  first <- posterior_predict(fit, grid, ndraws = 5, seed = 98)
  expect_identical(.Random.seed, before)
  expect_identical(first, posterior_predict(restored, grid, ndraws = 5, seed = 98))
  expect_equal(dim(posterior_epred(fit, grid, ndraws = 5, seed = 3)), c(5L, 2L))
  summary <- stats::predict(fit, grid, ndraws = 20, seed = 5)
  expect_true(all(summary$lower <= summary$upper))
  expect_equal(nrow(summary), 2L)
  expect_error(posterior_epred(fit, grid, ndraws = 21), "exceeds")
  expect_error(posterior_epred(fit, grid, draw_ids = 0), "existing")
  expect_s3_class(plot(fit, ndraws = 10), "ggplot")
})
