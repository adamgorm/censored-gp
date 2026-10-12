test_that("stored named-list data preserve matrix formula variables", {
  fit <- gaussian_fixture(grouped = TRUE, n_draws = 20L)
  original <- list(
    time = cbind(x1 = fit$data$time, x2 = c(0, 0, 1, 1, 0, 0)),
    y = fit$data$y, subject = factor(fit$data$subject), constant = 2
  )
  fit$prepared <- prepare_model(
    fit$formula, original, fit$priors, fit$kernel,
    "constant", "sum_to_zero", "quantile"
  )
  fit$data <- model_data_frame(original)
  expect_named(fit$data, names(original))
  expect_identical(fit$data$time, original$time)
  expect_identical(fit$data$subject, original$subject)
  expect_equal(fit$data$constant, rep(2, nrow(fit$data)))
  expect_identical(
    posterior_epred(fit, seed = 13),
    posterior_epred(fit, newdata = original, seed = 13)
  )
  expect_identical(
    posterior_predict(fit, seed = 14),
    posterior_predict(fit, newdata = original, seed = 14)
  )
  expect_identical(predict(fit, seed = 15), predict(fit, original, seed = 15))
  expect_identical(fitted(fit, seed = 16), predict(fit, original, seed = 16))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(fit, path)
  expect_identical(posterior_epred(readRDS(path), seed = 13), posterior_epred(fit, seed = 13))
})

test_that("storing data preserves existing data frames and censored responses", {
  data <- data.frame(time = 1:3, y = c(0.2, 0.3, 0.4))
  data$coordinates <- cbind(x1 = 1:3, x2 = 3:1)
  expect_identical(model_data_frame(data), data)
  response <- cens(data$y, c("none", "left", "right"))
  expect_identical(model_data_frame(list(time = data$time, response = response))$response, response)
  nested <- model_data_frame(list(
    extra = data.frame(index = 1:3, date = as.Date("2026-01-01") + 0:2),
    coordinates = data$coordinates, y = data$y
  ))
  expect_identical(nested$coordinates, data$coordinates)
  expect_identical(nested$extra.date, as.Date("2026-01-01") + 0:2)
})
