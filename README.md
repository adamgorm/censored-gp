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

## Optimization prototypes

The `research/` directory contains experimental implementations and reproducible
benchmarks for an exact, general censored-GP package. These are numerical
prototypes, not an installable R package yet. The original example above remains
the implementation accompanying the paper.

The new Stan model integrates out the GP function values and augments only the
censored measurements. Its joint target is the Gaussian density of the observed
and augmented measurements, restricted to the censoring bounds, times the prior
on the covariance parameters. Integrating out the augmented measurements recovers
the censored-data likelihood in the paper. A normalized truncated-normal density
must not replace this joint likelihood: its parameter-dependent normalizer would
remove information about the censoring event.

Two exact parameterizations support mixed left, right, and interval censoring,
observation-specific limits, irregular input locations, and entirely censored
data. One directly bounds the augmented measurements; the other uses a sequential
Gaussian quantile transformation with all required Jacobian terms. The latter
evaluates univariate log probabilities and quantiles, without numerically
integrating a multivariate Gaussian CDF. Both fit full Bayesian models; neither
fixes covariance parameters at empirical Bayes estimates.

A dependency-free C++ header (apart from Eigen) implements dense Gaussian
likelihoods, analytical first derivatives, and prediction using a cached Cholesky
factor. Rcpp and Stan adapters share that core. It supports squared exponential,
exponential, and Matern 3/2 and 5/2 kernels with arbitrary squared input distances
and supplied group covariance. The current censoring model uses the squared
exponential kernel and preserves the paper's sum-to-zero group covariance.
The custom adjoint currently accelerates the directly bounded parameterization;
the quantile parameterization uses Stan's derivatives.

Run these commands from the repository root, with `Rcpp`, `RcppEigen`,
`numDeriv`, `mvtnorm`, and `cmdstanr` installed and CmdStan configured:

```sh
Rscript research/validate-cpp.R
Rscript research/validate-censoring.R
Rscript research/benchmark-cpp.R
Rscript research/benchmark-censoring.R
```

The C++ benchmark compares identical likelihoods and gradients against a
vectorized Stan model compiled with `O1`. The censoring benchmark compares
effective sample size per second for two parameterizations and the custom
adjoint, with identical simulated data, priors, and iteration counts. Compilation
is excluded from timings. Results are short pilots on one machine; they do not
establish a universally fastest method or replace convergence assessment.

Prediction can return marginal means and variances without constructing a full
prediction covariance. Joint predictions retain the full conditional covariance;
independently sampling the marginals would lose dependence between predictions.
A cached factor is valid only for its covariance parameters and augmented data.
No low-rank kernel, independence approximation, or dataset-specific grid
structure is used. General dense GP factorization still takes cubic time and
quadratic memory.
