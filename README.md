# Gaussian process regression for value-censored functional and longitudinal data

This repository contains code for running value-censored Gaussian process
regression as described in the paper:

A. G. Hoffmann, C. T. Ekstrøm, B. Z. Christoffersen, and A. K. Jensen, *Gaussian
Process Regression for Value-Censored Functional and Longitudinal Data*,
Statistics in Medicine 44, no. 20-22 (2025): e70277,
https://doi.org/10.1002/sim.70277.

- `data_example.R` draws from and plots the posterior distribution in the HIV-1
  RNA data set from section 5 of the paper.
- `functions.stan` contains Stan functions for calculating things needed for
  drawing from the posterior, e.g., conditional means and covariances.
- `A5055data.txt` contains the HIV-1 RNA data set.

To run the code, you must install

- [`cmdstan`](https://mc-stan.org/users/interfaces/cmdstan)
- [`cmdstanr`](https://mc-stan.org/cmdstanr/)

and the CRAN packages

- `tidyverse`
- `mvtnorm`
- `DEoptim`
- `TruncatedNormal`
- `abind`

## R package: censoredGP

Full Bayesian Gaussian process regression for value-censored functional and
longitudinal data. Each observation can be exact, left censored, or right censored,
with its own limit. Models can describe a single curve
or a population curve with subject deviations.

The package implements the framework described in the paper above and was added
after publication. The original empirical Bayes example is preserved in
`data_example.R`.

### Install

Install R >= 4.1, a C++17 compiler, CmdStanR, and CmdStan >= 2.37:

```r
install.packages("cmdstanr", repos = c(
  "https://stan-dev.r-universe.dev", getOption("repos")
))
cmdstanr::install_cmdstan()
install.packages("remotes")
remotes::install_github("adamgorm/censored-gp", build_vignettes = TRUE)
```

The first fit compiles the model into a user cache. Subsequent fits reuse the
executable, including when priors, kernels, groups, or censoring limits change.
The repository owner must agree a redistribution license before a public
package repository release; the upstream repository currently specifies none.

### Observation-specific censoring

```r
library(censoredGP)

d <- data.frame(
  time = c(0, 0.2, 0.5, 1),
  value = c(1.2, 0.4, 1.5, 0.8),
  type = c("none", "left", "right", "left")
)

fit <- cgp(cens(value, type) ~ time, d,
  priors = cgp_priors(
    intercept = normal_prior(1, 1),
    amplitude = student_t_prior(4, 0, 1),
    lengthscale = lognormal_prior(log(0.5), 0.7),
    sigma = normal_prior(0, 0.5)
  ), seed = 123
)
```

`value` means the exact measurement for `"none"`, the upper detection limit for
`"left"`, or the lower detection limit for `"right"`. Each row has one value and
one type; limits can differ between observations. A logical
censoring column is also accepted: `TRUE` means left censored, `FALSE` exact.
Missing measurements are not censoring; omit unmeasured rows. Interval censoring,
which is also treated in the paper, is not currently supported by the package.

Use `cens(value, type) ~ time | subject` for hierarchical curves, a
plain `value ~ time` response for entirely observed data, or
`cens(value, type) ~ cbind(x1, x2) | subject` for multidimensional inputs.
Kernels are `"exp_quad"`, `"exponential"`, `"matern32"`, and `"matern52"`.

### Priors, predictions, and diagnostics

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
zero. By default, a constant intercept is estimated with a normal or Student-t
prior. Use `mean = 0` to fix the population GP prior mean to zero, as in the
original example. Default priors assume
approximately unit-scale responses and inputs; select substantive priors for
your problem and check them with `prior_predictive()` before fitting. No priors
or covariance parameters are estimated from the responses before sampling.

Subject deviations default to the paper's sum-to-zero covariance. Use
`group_structure = "independent"` to model independent deviations and predict
new subjects. `component = "population"` or `"deviation"` in `posterior_epred()`
selects those components; a component vector permits a joint request for both.

### HIV example and numerical implementation

```r
hiv <- hiv_a5055()
vignette("hiv-functional-data", package = "censoredGP")
```

The vignette fits the manuscript's 365 measurements from 44 patients with
explicit priors and full Bayesian inference. It is not an exact numerical
reproduction of the original empirical Bayes analysis. Its bundled figures and
diagnostics come from a real four-chain Stan run. It includes the complete fitting code;
set `CENSOREDGP_RUN_STAN=true` when rebuilding it to rerun sampling.

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
