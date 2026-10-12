test_that("observation-specific censoring produces the intended bounds", {
  response <- cens(c(1, 2, 3, 4, 0.5, -1),
    c("none", "left", "right", "interval", "left", "right"),
    upper = c(NA, NA, NA, 5, NA, NA)
  )
  expect_equal(as.numeric(response[, 1]), c(1, -Inf, 3, 4, -Inf, -1))
  expect_equal(as.numeric(response[, 2]), c(1, 2, Inf, 5, 0.5, Inf))
  expect_equal(cens(c(1, 2), c(FALSE, TRUE))[, 2], c(1, 2))
  expect_error(cens(NA), "finite")
  expect_error(cens(1, "interval", 0), "exceed")
  expect_error(cens(1, "interval"), "upper")
  expect_error(cens(c(1, 2), c("left", "unknown")), "every observation")
})

test_that("formula preparation retains mixed mechanisms and row identities", {
  data <- data.frame(
    time = c(3, 1, 5, 2, 4), id = c("B", "A", "B", "A", "B"),
    y = c(0.4, 1, -0.2, 0.5, -0.1),
    kind = c("interval", "none", "right", "left", "left"), high = c(0.7, NA, NA, NA, NA)
  )
  prepared <- prepare_model(
    cens(y, kind, high) ~ time | id, data, cgp_priors(),
    "exp_quad", "constant", "sum_to_zero", "quantile"
  )
  expect_equal(prepared$order, c(2L, 1L, 3L, 4L, 5L))
  expect_equal(prepared$censor_rows, c(1L, 3L, 4L, 5L))
  expect_equal(prepared$stan_data$censor_kind, c(3L, 2L, 1L, 1L))
  expect_equal(prepared$stan_data$censor_upper, c(0.7, 0, 0.5, -0.1))
  expect_equal(prepared$stan_data$censor_lower, c(0.4, -0.2, 0, 0))
  expect_equal(prepared$stan_data$group_cov[1, 2], -1)
  expect_error(prepare_model(
    y ~ time, transform(data, time = NA), cgp_priors(),
    "exp_quad", 0, "sum_to_zero", "quantile"
  ), "finite")
})

test_that("prior constructors validate scalar parameters and use explicit scales", {
  expect_error(normal_prior(sd = 0), "positive")
  expect_error(normal_prior(mean = c(1, 2)), "scalars")
  expect_error(student_t_prior(df = -1), "positive")
  expect_error(gamma_prior(2, NA), "finite")
  expect_error(cgp_priors(sigma = "normal(0,1)"), "constructor")
  expect_error(cgp_priors(intercept = gamma_prior(2, 1)), "real line")
  expect_equal(normal_prior(c(center = 1), 2), normal_prior(1, 2))
  expect_equal(prior_quantile(normal_prior(0, 2), 0.5), stats::qnorm(0.75) * 2)
  expect_equal(prior_quantile(student_t_prior(4, 0, 2), 0.5), stats::qt(0.75, 4) * 2)
  expect_equal(prior_quantile(inv_gamma_prior(3, 2), 0.5), 1 / stats::qgamma(0.5, 3, rate = 2))
  expect_true(prior_quantile(normal_prior(-40, 1), 0.5) > 0)
  priors <- cgp_priors(lengthscale = lognormal_prior(log(60), 0.7))
  expect_equal(priors$lengthscale$arguments, c(meanlog = log(60), sdlog = 0.7))
  encoded <- encode_priors(priors)
  expect_equal(encoded$prior_family, c(1, 2, 3, 2, 3, 2))
  expect_equal(encoded$prior_arguments[3, 1:2], c(log(60), 0.7))
})

test_that("prior prediction is reproducible and does not depend on responses", {
  data <- data.frame(t = seq(0, 1, length.out = 8), y = 0, id = rep(1:2, 4))
  set.seed(45)
  seed_before <- .Random.seed
  first <- prior_predictive(y ~ t | id, data, ndraws = 5, seed = 67)
  expect_identical(.Random.seed, seed_before)
  second <- prior_predictive(y ~ t | id, transform(data, y = 100), ndraws = 5, seed = 67)
  expect_identical(first, second)
  expect_equal(dim(first$latent), c(5L, 8L))
  expect_true(all(is.finite(first$response)))
})

test_that("bundled HIV observations retain the assay threshold", {
  hiv <- hiv_a5055()
  expect_equal(nrow(hiv), 365L)
  expect_equal(nlevels(hiv$subject), 44L)
  expect_equal(sum(hiv$censoring == "left"), 110L)
  expect_true(all(hiv$log_rna[hiv$censoring == "left"] == log10(50)))
})
