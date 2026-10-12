test_that("automatic latent plots preserve transformed grouping labels", {
  for (structure in c("sum_to_zero", "independent")) {
    fit <- gaussian_fixture(grouped = TRUE, structure = structure, n_draws = 20L)
    fit$data$subject <- factor(fit$data$subject)
    fit$formula <- y ~ time | paste0(subject, "-patient")
    environment(fit$formula) <- baseenv()
    fit$prepared <- prepare_model(fit$formula, fit$data, fit$priors,
      fit$kernel, "constant", structure, "quantile"
    )
    grid <- data.frame(
      time = rep(seq(0, 1, length.out = 5L), 3L),
      subject = rep(c("A", "B", "C"), each = 5L)
    )
    automatic <- plot(fit, type = "latent", n = 5L, seed = 14L)
    explicit <- plot(fit, newdata = grid, type = "latent", seed = 14L)
    expect_equal(automatic$data, explicit$data)
    expect_equal(unique(automatic$data$group), c("A-patient", "B-patient", "C-patient"))
    expect_equal(automatic$layers[[3L]]$data$group, fit$prepared$group)

    selected <- plot(fit, type = "latent", groups = "B-patient", n = 5L, seed = 14L)
    explicit_selected <- plot(fit, newdata = grid[grid$subject == "B", ],
      type = "latent", seed = 14L
    )
    expect_equal(selected$data, explicit_selected$data)
    expect_equal(selected$layers[[3L]]$data$group, rep("B-patient", 2L))
  }
})

test_that("new transformed groups require explicit plotting data", {
  fit <- gaussian_fixture(grouped = TRUE, structure = "independent", n_draws = 10L)
  fit$formula <- y ~ time | paste0(subject, "-patient")
  environment(fit$formula) <- baseenv()
  fit$prepared <- prepare_model(fit$formula, fit$data, fit$priors,
    fit$kernel, "constant", "independent", "quantile"
  )
  expect_error(plot(fit, type = "latent", groups = "D-patient", n = 5L),
    "Supply `newdata`.*transformed grouping expression"
  )
  explicit <- plot(fit, newdata = data.frame(time = c(0, 1), subject = "D"),
    type = "latent"
  )
  expect_equal(explicit$data$group, rep("D-patient", 2L))
  expect_equal(nrow(explicit$layers[[3L]]$data), 0L)
})

test_that("automatic latent plots support ordinary and new independent groups", {
  fit <- gaussian_fixture(grouped = TRUE, structure = "independent", n_draws = 10L)
  fit$data$subject <- factor(fit$data$subject)
  automatic <- plot(fit, type = "latent", groups = c("B", "D"), n = 5L, seed = 18L)
  explicit <- plot(fit, newdata = data.frame(
    time = rep(seq(0, 1, length.out = 5L), 2L),
    subject = rep(c("B", "D"), each = 5L)
  ), type = "latent", seed = 18L)
  expect_equal(automatic$data, explicit$data)
  expect_equal(unique(automatic$data$group), c("B", "D"))
  expect_equal(automatic$layers[[3L]]$data$group, rep("B", 2L))

  fit <- gaussian_fixture(grouped = TRUE, n_draws = 10L)
  ordinary <- plot(fit, type = "latent", n = 5L)
  expect_equal(unique(ordinary$data$group), c("A", "B", "C"))
  expect_equal(nrow(ordinary$layers[[3L]]$data), nrow(fit$data))
  expect_error(plot(fit, type = "latent", groups = "D", n = 5L), "independent")
})
