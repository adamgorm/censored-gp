# A fixed-parameter Gaussian fixture isolates conditional prediction from MCMC
# error. It is not used by fitting and does not impose a package approximation.
gaussian_fixture <- function(grouped = FALSE, structure = "sum_to_zero", n_draws = 2000L) {
  data <- data.frame(
    time = c(0, 0.15, 0.4, 0.6, 0.8, 1),
    y = c(0.1, 0.4, 0.2, -0.1, 0.3, 0.2), subject = rep(c("A", "B", "C"), 2)
  )
  formula <- if (grouped) y ~ time | subject else y ~ time
  # Use a persistent environment, as a formula created in an ordinary user
  # script would. Do not serialize pkgload's temporary test environment.
  environment(formula) <- baseenv()
  prepared <- prepare_model(formula, data, cgp_priors(), "exp_quad", "constant", structure, "quantile")
  parameters <- c(
    amplitude = 1.2, lengthscale = 0.5,
    if (grouped) c(group_amplitude = 0.6, group_lengthscale = 0.3), sigma = 0.25, intercept = 0.2
  )
  draws <- array(rep(parameters, each = n_draws), c(n_draws, 1L, length(parameters)),
    dimnames = list(iteration = as.character(seq_len(n_draws)), chain = "1", variable = names(parameters))
  )
  structure(list(
    formula = formula, data = data, prepared = prepared,
    kernel = "exp_quad", priors = cgp_priors(), draws = posterior::as_draws_array(draws),
    parameter_names = names(parameters)
  ), class = "cgp_fit")
}
