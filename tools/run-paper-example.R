# Reproduce the full Bayesian analysis used in the HIV vignette.
# Install the current package first, then run this script from the repo root.
library(censoredGP)

hiv <- hiv_a5055()
priors <- cgp_priors(
  intercept = normal_prior(3, 1.5),
  amplitude = student_t_prior(4, 0, 2),
  lengthscale = lognormal_prior(log(60), 0.7),
  group_amplitude = normal_prior(0, 1),
  group_lengthscale = lognormal_prior(log(35), 0.7),
  sigma = normal_prior(0, 0.5)
)

fit_file <- Sys.getenv("CENSOREDGP_FIT_FILE", file.path(tempdir(), "hiv-fit.rds"))
# Reusing a saved fit regenerates only the vignette figures and summaries. It
# never silently changes priors or substitutes an empirical Bayes fit.
if (identical(Sys.getenv("CENSOREDGP_REUSE_FIT"), "true")) {
  fit <- readRDS(fit_file)
  stopifnot(identical(fit$priors, priors), identical(fit$data, hiv))
} else {
  fit <- cgp(cens(log_rna, censoring) ~ day | subject, hiv,
    priors = priors, kernel = "exp_quad", mean = "constant",
    group_structure = "sum_to_zero", parameterization = "quantile",
    chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000,
    adapt_delta = 0.95, max_treedepth = 12, seed = 1627, refresh = 100,
    sig_figs = 18
  )
  saveRDS(fit, fit_file)
}
print(fit)

grid <- data.frame(day = seq(min(hiv$day), max(hiv$day), length.out = 80))
individual_grid <- do.call(rbind, lapply(c("1", "9", "41"), function(subject) {
  transform(grid, subject = factor(subject, levels = levels(hiv$subject)))
}))
population_draws <- posterior_epred(fit, grid, component = "population", ndraws = 800, seed = 53)
individual_draws <- posterior_epred(fit, individual_grid, ndraws = 800, seed = 54)
replicated <- posterior_predict(fit, ndraws = 400, seed = 55)
prior <- prior_predictive(log_rna ~ day, transform(grid, log_rna = 0), priors,
  ndraws = 200, seed = 56
)
imputed <- impute_censored(fit)
left <- hiv$censoring == "left"
stopifnot(
  all(imputed[, left, drop = FALSE] <= log10(50)),
  all(imputed[, !left, drop = FALSE] == rep(hiv$log_rna[!left], each = nrow(imputed)))
)

# Record the source actually used for sampling, even if the package has since
# received a formatting or documentation update.
stan_file <- tempfile(fileext = ".stan")
writeLines(fit$stan_code, stan_file)
stan_md5 <- unname(tools::md5sum(stan_file))
unlink(stan_file)
report <- list(
  summary = summary(fit), priors = prior_summary(fit), grid = grid,
  censored_diagnostics = posterior::summarise_draws(posterior::as_draws(fit, variables = "y_cens")),
  individual_grid = individual_grid, population_draws = population_draws,
  individual_draws = individual_draws, prior = prior,
  predictive_censoring_fraction = rowMeans(replicated < log10(50)),
  observed_censoring_fraction = mean(left),
  provenance = list(
    seed = 1627L, chains = 4L, iter_warmup = 1000L,
    iter_sampling = 1000L, adapt_delta = 0.95, kernel = "exp_quad",
    group_structure = "sum_to_zero", mean = "constant",
    cmdstan = fit$cmdstan_version, package = as.character(packageVersion("censoredGP")),
    stan_md5 = stan_md5, stan_source = fit$stan_code,
    cpp_md5 = stats::setNames(unname(fit$cpp_md5), basename(names(fit$cpp_md5))),
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
)
saveRDS(report, "inst/extdata/hiv-paper-analysis.rds", compress = "xz")
cat("Full fitted object:", fit_file, "\nVignette results: inst/extdata/hiv-paper-analysis.rds\n")
