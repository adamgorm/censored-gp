#' Simulate Joint Prior Predictive Curves and Measurements
#'
#' This function requires no Stan installation and never conditions on response
#' values or censoring events. The formula supplies only locations and grouping.
#' Use the resulting curves to assess whether the priors permit scientifically
#' plausible amplitudes, time scales, and noise levels before fitting.
#'
#' @inheritParams cgp
#' @param ndraws Number of prior draws.
#' @return A list with `latent` smooth-curve draws, `response` noisy-measurement
#'   draws, and `hyperparameters`. Matrices have draws in rows.
#' @examples
#' d <- data.frame(time = seq(0, 1, length.out = 10), y = 0)
#' p <- prior_predictive(y ~ time, d, ndraws = 5, seed = 27)
#' dim(p$latent)
#' @export
prior_predictive <- function(formula, data, priors = cgp_priors(), kernel = "exp_quad",
                             mean = "constant", group_structure = c("sum_to_zero", "independent"),
                             ndraws = 200L, seed = NULL) {
  ndraws <- check_count(ndraws, "ndraws")
  group_structure <- match.arg(group_structure)
  prepared <- prepare_model(formula, data, priors, kernel, mean, group_structure, "quantile")
  distance2 <- squared_distance(prepared$coordinates)
  group_cov <- group_covariance(
    prepared$group, prepared$group, length(prepared$group_levels),
    group_structure, prepared$has_group
  )
  with_local_seed(seed, {
    active_names <- c(
      "intercept", "amplitude", "lengthscale",
      if (prepared$has_group) c("group_amplitude", "group_lengthscale"), "sigma"
    )
    hyperparameters <- as.data.frame(lapply(active_names, function(name) {
      sample_prior(priors[[name]], ndraws, positive = name != "intercept")
    }))
    names(hyperparameters) <- active_names
    if (!prepared$estimate_mean) hyperparameters$intercept <- rep(mean, ndraws)
    latent <- response <- matrix(NA_real_, ndraws, nrow(prepared$coordinates))
    for (draw in seq_len(ndraws)) {
      theta <- hyperparameters[draw, ]
      covariance <- stationary_kernel(distance2, theta$amplitude, theta$lengthscale, kernel)
      if (prepared$has_group) {
        covariance <- covariance + group_cov *
          stationary_kernel(distance2, theta$group_amplitude, theta$group_lengthscale, kernel)
      }
      latent[draw, ] <- sample_gaussian(rep(theta$intercept, nrow(covariance)), covariance)
      response[draw, ] <- latent[draw, ] + stats::rnorm(ncol(latent), sd = theta$sigma)
    }
    list(latent = latent, response = response, hyperparameters = hyperparameters)
  })
}
