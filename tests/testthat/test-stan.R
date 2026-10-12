test_that("Stan samples left and right censoring with both parameterizations", {
  skip_if(Sys.getenv("CENSOREDGP_TEST_STAN") != "true", "Enable explicit Stan integration tests")
  skip_if_not_installed("cmdstanr")
  data <- data.frame(
    time = c(0, 0.1, 0.25, 0.4, 0.6, 0.75, 0.9, 1),
    y = c(0.2, -0.1, 0.1, -0.3, 0.4, 0.2, 0.5, -0.2),
    kind = rep(c("none", "left", "right", "left"), 2),
    subject = rep(c("z", "a"), 4)
  )
  priors <- cgp_priors(
    amplitude = lognormal_prior(-0.5, 0.5),
    lengthscale = inv_gamma_prior(3, 1), group_amplitude = exponential_prior(2),
    group_lengthscale = gamma_prior(2, 4), sigma = student_t_prior(4, 0, 0.5)
  )
  for (parameterization in c("bounded", "quantile")) {
    for (kernel in c("exp_quad", "exponential", "matern32", "matern52")) {
      fit <- cgp(cens(y, kind) ~ time | subject, data,
        priors = priors,
        parameterization = parameterization, kernel = kernel,
        chains = 2, parallel_chains = 2, iter_warmup = 250, iter_sampling = 100,
        adapt_delta = 0.95, seed = 914, refresh = 0, diagnostic_warnings = FALSE
      )
      imputed <- impute_censored(fit)
      expect_equal(imputed[, data$kind == "none"],
        matrix(rep(data$y[data$kind == "none"], each = nrow(imputed)), nrow = nrow(imputed)),
        ignore_attr = TRUE
      )
      for (row in seq_len(nrow(data))) {
        expect_true(all(imputed[, row] >= fit$prepared$bounds$lower[row]))
        expect_true(all(imputed[, row] <= fit$prepared$bounds$upper[row]))
      }
      expect_true(all(is.finite(posterior::as_draws_matrix(fit))))
      expect_true(all(is.finite(posterior_epred(fit, ndraws = 5, seed = 1))))
    }
  }
})

test_that("Stan handles entirely censored and entirely observed single curves", {
  skip_if(Sys.getenv("CENSOREDGP_TEST_STAN") != "true", "Enable explicit Stan integration tests")
  skip_if_not_installed("cmdstanr")
  for (types in list(rep("none", 4), c("left", "right", "right", "left"))) {
    data <- data.frame(
      time = c(0, 0.2, 0.7, 1), y = c(0.2, -0.1, -0.2, 0.3),
      type = types
    )
    fit <- cgp(cens(y, type) ~ time, data,
      mean = 0,
      chains = 2, parallel_chains = 2, iter_warmup = 250, iter_sampling = 100,
      adapt_delta = 0.95, seed = 817, refresh = 0, diagnostic_warnings = FALSE
    )
    expect_equal(fit$prepared$has_group, 0L)
    expect_false("intercept" %in% fit$parameter_names)
    expect_equal(ncol(impute_censored(fit)), 4L)
    expect_true(all(is.finite(posterior_predict(fit, ndraws = 5, seed = 2))))
  }
})
