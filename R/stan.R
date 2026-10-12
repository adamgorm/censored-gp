.compiled_models <- new.env(parent = emptyenv())

#' Inspect the Stan Model
#' @return A character string containing the package's complete Stan program.
#' @export
cgp_stancode <- function() {
  paste(readLines(system.file("stan", "censored_gp.stan", package = "censoredGP")), collapse = "\n")
}

#' Compile and Cache the Bayesian Model
#'
#' Compilation happens automatically on the first call to [cgp()]. Compiled
#' models are cached in a user-writable directory and reused across sessions.
#' Priors, kernels, and censoring bounds are data and do not cause recompilation.
#' The cache key includes source files, CmdStan version, and platform.
#'
#' @param force Recompile even if a cached executable exists.
#' @param quiet Suppress compiler progress messages.
#' @param cache_dir Directory for source and executable caches.
#' @return A `cmdstanr::CmdStanModel` object, invisibly.
#' @export
compile_cgp_model <- function(force = FALSE, quiet = TRUE,
                              cache_dir = tools::R_user_dir("censoredGP", "cache")) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Install cmdstanr from https://stan-dev.r-universe.dev, then run cmdstanr::install_cmdstan().",
      call. = FALSE
    )
  }
  version <- tryCatch(cmdstanr::cmdstan_version(), error = function(error) NULL)
  if (is.null(version) || version < "2.37.0") {
    stop("Configure CmdStan >= 2.37.0 with cmdstanr::install_cmdstan().", call. = FALSE)
  }
  stan_file <- system.file("stan", "censored_gp.stan", package = "censoredGP")
  include_dir <- system.file("include", "censoredGP", package = "censoredGP")
  sources <- c(stan_file, list.files(include_dir, full.names = TRUE, pattern = "\\.hpp$"))
  # Hash file contents, not installation paths, so moving an R library does not
  # needlessly invalidate an otherwise identical standalone executable.
  key_file <- tempfile()
  on.exit(unlink(key_file), add = TRUE)
  writeLines(c(unname(tools::md5sum(sources)), as.character(version), R.version$platform), key_file)
  key <- unname(tools::md5sum(key_file))
  if (!force && exists(key, envir = .compiled_models, inherits = FALSE)) {
    return(invisible(get(key, envir = .compiled_models)))
  }
  model_dir <- file.path(cache_dir, key)
  if (!dir.exists(model_dir) && !dir.create(model_dir, recursive = TRUE)) {
    stop("Could not create the model cache directory: ", model_dir, call. = FALSE)
  }
  target <- file.path(model_dir, "censored_gp.stan")
  if (!file.exists(target)) file.copy(stan_file, target)
  model <- cmdstanr::cmdstan_model(target,
    user_header = file.path(include_dir, "stan_adapter.hpp"),
    stanc_options = list("O1" = TRUE, "allow-undefined" = TRUE),
    force_recompile = force, quiet = quiet
  )
  assign(key, model, envir = .compiled_models)
  invisible(model)
}
