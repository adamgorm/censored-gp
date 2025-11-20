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
