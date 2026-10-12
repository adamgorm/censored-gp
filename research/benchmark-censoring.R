# Run from the repository root: Rscript research/benchmark-censoring.R
# Same priors, data, initialization, and iteration counts for every method.
# Short pilots measure effective sample size per second, not convergence proof.
library(cmdstanr)
p <- tempfile("censored-gp-censoring-benchmark-")
dir.create(p)
file.copy("research/censored-gp.stan", p)
file.copy("research/censor-transform.staninc", p)
model <- cmdstan_model(file.path(p, "censored-gp.stan"),
  stanc_options = list("allow-undefined" = TRUE),
  user_header = normalizePath("research/cpp/gaussian_stan.hpp"), quiet = TRUE
)
make_data <- function(N, G, scheme, seed) {
  set.seed(seed)
  t <- runif(N)
  group <- rep(seq_len(G), length.out = N)
  d2 <- outer(t, t, "-")^2
  B <- if (G > 1) ifelse(outer(group, group, "=="), 1, -1 / (G - 1)) else matrix(0, N, N)
  K <- 1.2^2 * exp(-d2 / (2 * .25^2)) + B * (.6^2 * exp(-d2 / (2 * .15^2))) + diag(.25^2, N)
  y <- as.vector(t(chol(K)) %*% rnorm(N))
  kind <- integer(N)
  lower <- rep(-Inf, N)
  upper <- rep(Inf, N)
  if (scheme != "none") {
    eligible <- sample(1:3, N, replace = TRUE)
    cut <- if (scheme == "heavy") .9 else -.3
    for (i in seq_len(N)) {
      if (eligible[i] == 1 && y[i] < cut + .2 * sin(8 * t[i])) {
        kind[i] <- 1L
        upper[i] <- cut + .2 * sin(8 * t[i])
      }
      if (eligible[i] == 2 && y[i] > -cut + .2 * cos(6 * t[i])) {
        kind[i] <- 2L
        lower[i] <- -cut + .2 * cos(6 * t[i])
      }
      if (eligible[i] == 3 && runif(1) < if (scheme == "heavy") .9 else .5) {
        kind[i] <- 3L
        lower[i] <- floor(y[i] / .4) * .4
        upper[i] <- lower[i] + .4
      }
    }
  }
  ix <- order(kind)
  list(N = N, G = G, No = sum(kind == 0), Nl = sum(kind == 1), Nr = sum(kind == 2), Ni = sum(kind == 3), yo = y[kind == 0], left_upper = upper[kind == 1], right_lower = lower[kind == 2], interval_lower = lower[kind == 3], interval_upper = upper[kind == 3], t = t[ix], group = group[ix])
}
cases <- list(list(name = "single-observed", N = 60L, G = 1L, scheme = "none"), list(name = "single-mixed", N = 80L, G = 1L, scheme = "mixed"), list(name = "grouped-mixed", N = 160L, G = 5L, scheme = "mixed"), list(name = "grouped-heavy", N = 160L, G = 5L, scheme = "heavy"))
results <- list()
for (k in seq_along(cases)) {
  cc <- cases[[k]]
  dat <- make_data(cc$N, cc$G, cc$scheme, 1700 + k)
  cat("DESIGN", cc$name, "observed/left/right/interval", dat$No, dat$Nl, dat$Nr, dat$Ni, "\n")
  for (method in c("bounded", "custom", "quantile")) {
    parameterization <- as.integer(method == "quantile")
    dat$custom_gradient <- as.integer(method == "custom")
    dat$parameterization <- parameterization
    Nc <- dat$Nl + dat$Nr + dat$Ni
    initial <- function() list(alpha_mu = 1.2, rho_mu = .25, alpha_eta = .6, rho_eta = .15, sigma = .25, yc_left = if (parameterization == 0) dat$left_upper - .2 else numeric(), yc_right = if (parameterization == 0) dat$right_lower + .2 else numeric(), yc_interval = if (parameterization == 0) (dat$interval_lower + dat$interval_upper) / 2 else numeric(), u = if (parameterization == 1) rep(.5, Nc) else numeric())
    name <- paste0(cc$name, "-", method)
    elapsed <- system.time(fit <- model$sample(data = dat, init = initial, seed = 1935, chains = 2, parallel_chains = 2, iter_warmup = 250, iter_sampling = 250, adapt_delta = .9, max_treedepth = 10, refresh = 0, output_dir = p, output_basename = name))[["elapsed"]]
    params <- c("alpha_mu", "rho_mu", if (cc$G > 1) c("alpha_eta", "rho_eta"), "sigma")
    summ <- fit$summary(params)
    diag <- fit$diagnostic_summary()
    sd <- fit$sampler_diagnostics(format = "df")
    row <- data.frame(design = cc$name, parameterization = method, N = cc$N, censored = Nc, elapsed = elapsed, min_ess = min(summ$ess_bulk), min_ess_per_second = min(summ$ess_bulk) / elapsed, max_rhat = max(summ$rhat), divergences = sum(diag$num_divergent), treedepth = sum(diag$num_max_treedepth), mean_leapfrog = mean(sd$n_leapfrog__))
    print(row)
    print(summ)
    results[[length(results) + 1]] <- row
    saveRDS(list(row = row, summary = summ, diagnostics = diag, output = fit$output_files(), data = dat), file.path(p, paste0(name, "-results.rds")))
    write.csv(do.call(rbind, results), "research/censoring-benchmarks.csv", row.names = FALSE)
  }
}
cat("GENERAL BENCHMARKS COMPLETE\n")
print(do.call(rbind, results))
