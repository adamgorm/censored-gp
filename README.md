# Gaussian process regression for value-censored functional and longitudinal data

- `data_example.R` draws from and plots the posterior distribution in the HIV-1
  RNA data set.
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
