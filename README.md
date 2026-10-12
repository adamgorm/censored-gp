# censoredGP

Full Bayesian Gaussian process regression for value-censored functional and
longitudinal data. Each observation can be exact, left censored, right censored,
or interval censored, with its own limits. Models can describe a single curve
or a population curve with subject deviations.

The package extends the model in Hoffmann, Ekstrøm, Christoffersen, and Jensen,
[*Gaussian Process Regression for Value-Censored Functional and Longitudinal
Data*](https://doi.org/10.1002/sim.70277), Statistics in Medicine (2025).
The original empirical Bayes example is preserved in `data_example.R`.

## Install

Install R >= 4.1, a C++17 compiler, CmdStanR, and CmdStan >= 2.37:

```r
install.packages("cmdstanr", repos = c(
  "https://stan-dev.r-universe.dev", getOption("repos")
))
cmdstanr::install_cmdstan()
install.packages("remotes")
remotes::install_github("aejensen/censored-gp", ref = "codex/repo-updates",
                        build_vignettes = TRUE)
```

The first fit compiles the model into a user cache. Subsequent fits reuse the
executable, including when priors, kernels, groups, or censoring limits change.
The repository owner must agree a redistribution license before a public
package repository release; the upstream repository currently specifies none.

## Observation-specific censoring

```r
library(censoredGP)

d <- data.frame(
  time = c(0, 0.2, 0.5, 1),
  value = c(1.2, 0.4, 1.5, 0.8),
  type = c("none", "left", "right", "interval"),
  upper = c(NA, NA, NA, 1.1)
)

fit <- cgp(cens(value, type, upper) ~ time, d,
  priors = cgp_priors(
    intercept = normal_prior(1, 1),
    amplitude = student_t_prior(4, 0, 1),
    lengthscale = lognormal_prior(log(0.5), 0.7),
    sigma = normal_prior(0, 0.5)
  ), seed = 123
)
```

`value` means the exact measurement for `"none"`, the upper limit for `"left"`,
the lower limit for `"right"`, and the lower endpoint for `"interval"`.
`upper` supplies interval upper endpoints and is ignored for the other types.
Rows may mix all four types, and limits may differ within each type. A logical
censoring column is also accepted: `TRUE` means left censored, `FALSE` exact.
Missing measurements are not censoring; omit unmeasured rows.

Use `cens(value, type, upper) ~ time | subject` for hierarchical curves, a
plain `value ~ time` response for entirely observed data, or
`cens(value, type, upper) ~ cbind(x1, x2) | subject` for multidimensional inputs.
Kernels are `"exp_quad"`, `"exponential"`, `"matern32"`, and `"matern52"`.

## Priors, predictions, and diagnostics

```r
prior_summary(fit)
summary(fit)                       # R-hat and effective sample sizes
diagnostics(fit)                   # divergences, tree-depth hits, BFMI

grid <- data.frame(time = seq(0, 1, length.out = 80))
predict(fit, grid)                 # pointwise smooth-curve summaries
curves <- posterior_epred(fit, grid, seed = 24) # joint smooth-curve draws
future <- posterior_predict(fit, grid, seed = 25) # joint noisy measurements
imputed <- impute_censored(fit)    # original measurements, respecting row bounds
plot(fit)

saveRDS(fit, "fit.rds")            # predictions do not require Stan CSV files
```

Priors are validated R objects, specified by component in the original units.
Normal, Student-t, lognormal, gamma, inverse-gamma, and exponential families are
available. Normal and Student-t priors for positive parameters are truncated at
zero. An estimated intercept has a normal or Student-t prior. Use `mean = 0` to
fix the population GP mean to the paper's specification. Default priors assume
approximately unit-scale responses and inputs; select substantive priors for
your problem and check them with `prior_predictive()` before fitting. No priors
or covariance parameters are estimated from the responses before sampling.

Subject deviations default to the paper's sum-to-zero covariance. Use
`group_structure = "independent"` to model independent deviations and predict
new subjects. `component = "population"` or `"deviation"` in `posterior_epred()`
selects those components; a component vector permits a joint request for both.

## HIV example and numerical implementation

```r
hiv <- hiv_a5055()
vignette("hiv-functional-data", package = "censoredGP")
```

The vignette fits the manuscript's 365 measurements from 44 patients with
explicit priors and full Bayesian inference. Its bundled figures and diagnostics
come from a real four-chain Stan run. It includes the complete fitting code;
set `CENSOREDGP_RUN_STAN=true` when rebuilding it to rerun sampling. Maintainers
can regenerate the bundled results with `Rscript tools/run-paper-example.R`.

The dense GP function values are integrated out analytically. Only censored
noisy measurements are augmented. A shared Eigen C++ core supplies the likelihood,
analytical gradients, and conditional predictions. The default exact sequential
Gaussian quantile transformation includes its full Jacobian; a directly bounded
parameterization is also available. Neither requires an approximate multivariate
normal CDF, low-rank covariance, or an empirical Bayes plug-in estimate.

Joint predictions preserve covariance across locations and components. Pointwise
summaries avoid constructing an unnecessary prediction covariance. Dense GP
factorization still costs O(n³) time and O(n²) memory. Computation does not assume
regular grids, shared censoring limits, or the HIV data layout.

For development:

```sh
Rscript -e 'testthat::test_local()'
CENSOREDGP_TEST_STAN=true Rscript -e 'testthat::test_local()'
Rscript tools/validate-core.R
Rscript tools/validate-stan-gradients.R
R CMD build .
R CMD check --no-manual censoredGP_0.1.0.tar.gz
```

`research/` preserves the earlier optimization experiments and their benchmark
results. The installed implementation is in `R/`, `src/`, and `inst/`.
