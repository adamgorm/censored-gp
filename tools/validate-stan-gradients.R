# Check the actual compiled Stan target, including all prior, constraint, and
# external C++ adjoint contributions. The core is checked separately against R.
# Run from the repository root after installing the package.
library(censoredGP)
library(cmdstanr)

model <- compile_cgp_model()
output <- Sys.getenv("CENSOREDGP_GRADIENT_DIR", tempfile("censoredGP-stan-gradients-"))
dir.create(output, showWarnings = FALSE)
data <- data.frame(
  time = c(0, 0.1, 0.25, 0.4, 0.6, 0.75, 0.9, 1),
  y = c(0.2, -0.1, 0.1, -0.3, 0.4, 0.2, 0.5, -0.2),
  kind = rep(c("none", "left", "right", "interval"), 2),
  upper = c(NA, NA, NA, 0.1, NA, NA, NA, 0.3),
  subject = rep(c("z", "a"), 4)
)
priors <- cgp_priors(
  intercept = student_t_prior(4, 0.1, 1),
  amplitude = lognormal_prior(-0.5, 0.5),
  lengthscale = inv_gamma_prior(3, 1),
  group_amplitude = exponential_prior(2),
  group_lengthscale = gamma_prior(2, 4),
  sigma = normal_prior(0, 0.5)
)
worst <- 0
for (parameterization in c("bounded", "quantile")) {
  for (kernel in c("exp_quad", "exponential", "matern32", "matern52")) {
    prepared <- censoredGP:::prepare_model(
      cens(y, kind, upper) ~ time | subject, data, priors,
      kernel, "constant", "sum_to_zero", parameterization
    )
    # All free Stan coordinates are unconstrained here. In particular the last
    # six entries are logit quantiles or unconstrained bounded measurements.
    state <- c(
      log(c(1.2, 0.6, 0.7, 0.4, 0.25)), 0.2,
      seq(-0.7, 0.7, length.out = 6)
    )
    step <- 1e-5
    perturbation <- diag(step, length(state))
    baseline <- matrix(rep(state, each = length(state)), nrow = length(state))
    states <- rbind(state, baseline + perturbation, baseline - perturbation)
    data_file <- file.path(output, "data.json")
    params_file <- file.path(output, "params.json")
    results_file <- file.path(output, "results.csv")
    write_stan_json(prepared$stan_data, data_file)
    arguments <- c(
      "log_prob", paste0("unconstrained_params=", shQuote(params_file)),
      "jacobian=1", "data", paste0("file=", shQuote(data_file)),
      "output", paste0("file=", shQuote(results_file)), "sig_figs=18", "refresh=0"
    )
    # Pass one vector at a time. CmdStan 2.37's get_uparams_r() reshapes a batch
    # in row-major order although JSON var_context stores column-major values.
    # A single vector avoids that version-specific batch-input ambiguity.
    result <- do.call(rbind, lapply(seq_len(nrow(states)), function(row) {
      write_stan_json(list(params_r = states[row, ]), params_file)
      status <- system2(model$exe_file(), arguments, stdout = FALSE, stderr = FALSE)
      stopifnot(status == 0)
      as.matrix(read.csv(results_file, comment.char = "#", check.names = FALSE))
    }))
    numerical <- (result[1 + seq_along(state), 1] -
      result[1 + length(state) + seq_along(state), 1]) / (2 * step)
    analytical <- result[1, -1]
    error <- max(abs(analytical - numerical) / pmax(1, abs(numerical)))
    cat(parameterization, kernel, "maximum relative gradient error:", error, "\n")
    if (error >= 2e-6) print(data.frame(analytical, numerical))
    stopifnot(error < 2e-6)
    worst <- max(worst, error)
  }
}
cat("Compiled Stan target: PASS. Maximum relative gradient error:", worst, "\n")
unlink(output, recursive = TRUE)
