library(censoredGP)
set.seed(451)
d <- data.frame(
  time = runif(12), subject = rep(1:3, 4),
  y = rep(c(0.2, -0.1, 0.3, -0.2), 3),
  type = rep(c("none", "left", "right", "interval"), 3), upper = 0.4
)
fit <- cgp(cens(y, type, upper) ~ time | subject, d,
  chains = 2, parallel_chains = 2, iter_warmup = 150, iter_sampling = 100,
  seed = 711, refresh = 0, diagnostic_warnings = FALSE
)
stopifnot(all(is.finite(posterior::as_draws_matrix(fit))))
imputed <- impute_censored(fit)
stopifnot(
  all(imputed[, d$type == "left"] <= -0.1),
  all(imputed[, d$type == "right"] >= 0.3),
  all(imputed[, d$type == "interval"] >= -0.2),
  all(imputed[, d$type == "interval"] <= 0.4)
)
grid <- data.frame(time = seq(0, 1, length.out = 10), subject = 1)
stopifnot(all(is.finite(posterior_epred(fit, grid, ndraws = 10, seed = 1))))
stopifnot(inherits(plot(fit, ndraws = 10), "ggplot"))
saveRDS(fit, "/tmp/censoredGP-smoke-fit.rds")
cat("Mixed-censoring package smoke test: PASS\n")
