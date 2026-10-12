# Run from the repository root: Rscript research/validate-cpp.R
stopifnot(requireNamespace("Rcpp", quietly = TRUE),
          requireNamespace("RcppEigen", quietly = TRUE),
          requireNamespace("numDeriv", quietly = TRUE))
Rcpp::sourceCpp("research/cpp/gaussian_rcpp.cpp", cacheDir = tempdir())

kernel_reference <- function(distance2, amplitude, length, kernel) {
  r <- sqrt(distance2) / length
  amplitude^2 * switch(kernel,
    exp(-r^2/2), exp(-r),
    (1 + sqrt(3)*r)*exp(-sqrt(3)*r),
    (1 + sqrt(5)*r + 5*r^2/3)*exp(-sqrt(5)*r))
}
covariance_reference <- function(distance2, B, theta, kernel) {
  kernel_reference(distance2, theta[1], theta[2], kernel) +
    B * kernel_reference(distance2, theta[3], theta[4], kernel)
}
log_density_reference <- function(theta, y, distance2, B, kernel) {
  K <- covariance_reference(distance2, B, theta, kernel) + diag(theta[5]^2, length(y))
  U <- chol(K)
  z <- forwardsolve(t(U), y)
  -length(y)*log(2*pi)/2 - sum(log(diag(U))) - sum(z*z)/2
}

set.seed(6971)
worst_gradient <- worst_prediction <- worst_density <- 0
for (N in c(1L, 13L, 31L)) for (G in c(1L, 2L, 5L)) for (kernel in 1:4) {
  M <- 9L
  # Independent arbitrary 2D inputs, including an input duplicated at prediction.
  x <- matrix(runif((N + M)*2), ncol = 2)
  x[N + 1, ] <- x[1, ]
  d2 <- as.matrix(dist(x))^2
  group <- sample.int(G, N + M, replace = TRUE)
  B <- if (G == 1) matrix(0, N + M, N + M) else
    ifelse(outer(group, group, "=="), 1, -1/(G - 1))
  ii <- seq_len(N); jj <- N + seq_len(M)
  theta <- c(1.2, .6, .7, .4, .25)
  residual <- rnorm(N)
  fit <- gaussian_evaluate_cpp(d2[ii, ii, drop = FALSE], B[ii, ii, drop = FALSE],
    theta, residual, kernel)
  reference <- function(v) log_density_reference(v[1:5], v[-(1:5)],
    d2[ii, ii, drop = FALSE], B[ii, ii, drop = FALSE], kernel)
  fd <- numDeriv::grad(reference, c(theta, residual))
  analytical <- c(fit$gradient_theta, fit$gradient_y)
  gradient_error <- max(abs(fd - analytical) / pmax(1, abs(fd)))
  density_error <- abs(fit$log_density - reference(c(theta, residual)))
  stopifnot(gradient_error < 2e-6, density_error < 1e-10)

  K <- covariance_reference(d2, B, theta, kernel)
  Kyy <- K[ii, ii, drop = FALSE] + diag(theta[5]^2, N)
  cross <- K[jj, ii, drop = FALSE]
  prior <- K[jj, jj, drop = FALSE]
  reference_mean <- as.vector(cross %*% solve(Kyy, residual))
  reference_covariance <- prior - cross %*% solve(Kyy, t(cross))
  joint <- gaussian_predict_joint_cpp(fit$factor, cross, prior)
  marginal <- gaussian_predict_marginal_cpp(fit$factor, cross, diag(prior))
  prediction_error <- max(abs(joint$mean - reference_mean),
    abs(joint$covariance - reference_covariance),
    abs(marginal$mean - reference_mean),
    abs(marginal$variance - diag(reference_covariance)))
  stopifnot(prediction_error < 1e-10)
  worst_gradient <- max(worst_gradient, gradient_error)
  worst_density <- max(worst_density, density_error)
  worst_prediction <- max(worst_prediction, prediction_error)
}
# Failed scales must be rejected rather than silently repaired with jitter.
stopifnot(inherits(try(gaussian_evaluate_cpp(matrix(0,1,1), matrix(0,1,1),
  c(1,1,1,1,0), 0), silent = TRUE), "try-error"))
cat("PASS: 36 arbitrary-input designs and four kernels\n",
    "Maximum relative gradient error:", worst_gradient, "\n",
    "Maximum log-density error:", worst_density, "\n",
    "Maximum prediction error:", worst_prediction, "\n")
