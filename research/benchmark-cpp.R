# Run from the repository root: Rscript research/benchmark-cpp.R
# Both implementations evaluate the same augmented Gaussian density. This
# measures gradient cost, not censoring integration or posterior mixing.
library(cmdstanr)
output <- tempfile("censored-gp-cpp-benchmark-")
dir.create(output)
header <- normalizePath("research/cpp/gaussian_stan.hpp")
models <- list()
for (method in c("reference", "custom")) {
  file.copy(paste0("research/gaussian-", method, ".stan"), output)
  models[[method]] <- cmdstan_model(file.path(output, paste0("gaussian-", method, ".stan")),
    stanc_options = if (method == "custom") list("O1" = TRUE, "allow-undefined" = TRUE) else list("O1" = TRUE),
    user_header = if (method == "custom") header else NULL, quiet = TRUE
  )
}
set.seed(1851)
results <- list()
for (N in c(80L, 160L, 320L)) {
  for (kernel in 1:4) {
    No <- N %/% 2L
    G <- 5L
    x <- matrix(runif(N * 2), ncol = 2)
    distance2 <- as.matrix(dist(x))^2
    group <- sample.int(G, N, replace = TRUE)
    B <- ifelse(outer(group, group, "=="), 1, -1 / (G - 1))
    observed <- rnorm(No)
    data_file <- file.path(output, "data.json")
    write_stan_json(list(
      N = N, No = No, kernel = kernel, distance2 = distance2,
      group_cov = B, observed = observed
    ), data_file)
    repetitions <- 100L
    base <- c(log(c(1.2, .6, .7, .4, .25)), rnorm(N - No))
    states <- matrix(rep(base, each = repetitions), nrow = repetitions) +
      matrix(rnorm(repetitions * length(base), sd = .01), nrow = repetitions)
    params_file <- file.path(output, "params.json")
    write_stan_json(list(params_r = states), params_file)
    timing <- list()
    values <- list()
    for (method in names(models)) {
      out <- file.path(output, paste0(method, ".csv"))
      args <- c(
        "log_prob", paste0("unconstrained_params=", shQuote(params_file)),
        "data", paste0("file=", shQuote(data_file)), "output",
        paste0("file=", shQuote(out)), "sig_figs=18", "refresh=0"
      )
      # Warm caches, then take the median of three runs including process and I/O.
      run <- function() {
        status <- system2(models[[method]]$exe_file(), args, stdout = FALSE, stderr = FALSE)
        stopifnot(status == 0)
      }
      run()
      timing[[method]] <- median(replicate(3, system.time(run())[["elapsed"]]))
      values[[method]] <- as.matrix(read.csv(out, comment.char = "#", check.names = FALSE))
    }
    relative_error <- max(abs(values$custom - values$reference) / pmax(1, abs(values$reference)))
    stopifnot(relative_error < 1e-8)
    row <- data.frame(
      N = N, kernel = kernel, evaluations = repetitions,
      reference_seconds = timing$reference, custom_seconds = timing$custom,
      speedup = timing$reference / timing$custom, max_relative_error = relative_error
    )
    print(row)
    results[[length(results) + 1L]] <- row
  }
}
write.csv(do.call(rbind, results), "research/cpp-benchmarks.csv", row.names = FALSE)
cat("Artifacts:", output, "\n")
