functions {
  real cgp_log_density(vector state, data matrix distance2, data matrix group_cov,
                       data vector observed, data vector lower_bound,
                       data vector upper_bound, data int kernel, data int quantile);
  vector cgp_imputed_values(vector state, data matrix distance2, data matrix group_cov,
                            data vector observed, data vector lower_bound,
                            data vector upper_bound, data int kernel);

  // Priors are data, so changing a distribution or its parameters does not
  // require recompiling the model. Positive normal/t priors are truncated at 0.
  real hyperprior_lpdf(real value, int family, vector arguments, int positive) {
    if (family == 1) {
      real result = normal_lpdf(value | arguments[1], arguments[2]);
      if (positive)
        result -= normal_lcdf(0 | -arguments[1], arguments[2]);
      return result;
    }
    if (family == 2) {
      real result = student_t_lpdf(value | arguments[1], arguments[2], arguments[3]);
      if (positive)
        result -= student_t_lccdf(0 | arguments[1], arguments[2], arguments[3]);
      return result;
    }
    if (family == 3)
      return lognormal_lpdf(value | arguments[1], arguments[2]);
    if (family == 4)
      return gamma_lpdf(value | arguments[1], arguments[2]);
    if (family == 5)
      return inv_gamma_lpdf(value | arguments[1], arguments[2]);
    return exponential_lpdf(value | arguments[1]);
  }
}
data {
  int<lower=1> N;
  int<lower=0, upper=N> n_observed;
  int<lower=0, upper=1> has_group;
  int<lower=0, upper=1> estimate_mean;
  real fixed_mean;
  int<lower=1, upper=4> kernel;
  int<lower=0, upper=1> quantile;
  matrix[N, N] distance2;
  matrix[N, N] group_cov;
  vector[n_observed] observed;
  array[N - n_observed] int<lower=1, upper=3> censor_kind;
  vector[N - n_observed] censor_lower;
  vector[N - n_observed] censor_upper;
  array[6] int<lower=1, upper=6> prior_family;
  matrix[6, 3] prior_arguments;
}
transformed data {
  int n_censored = N - n_observed;
  vector[n_censored] lower_bound;
  vector[n_censored] upper_bound;
  if (n_censored > 0)
    for (k in 1 : n_censored) {
      lower_bound[k] = censor_kind[k] == 1 ? negative_infinity() : censor_lower[k];
      upper_bound[k] = censor_kind[k] == 2 ? positive_infinity() : censor_upper[k];
    }
  vector[quantile ? 0 : n_censored] bounded_lower = head(lower_bound,
                                                         quantile ? 0 : n_censored);
  vector[quantile ? 0 : n_censored] bounded_upper = head(upper_bound,
                                                         quantile ? 0 : n_censored);
}
parameters {
  // Group covariance parameters do not exist in an ungrouped model.
  vector<lower=0>[3 + 2 * has_group] covariance_parameters;
  vector[estimate_mean] mean_parameter;
  vector<lower=bounded_lower, upper=bounded_upper>[quantile ? 0 : n_censored] censored;
  vector<lower=0, upper=1>[quantile ? n_censored : 0] unit;
}
transformed parameters {
  vector[5] theta;
  real intercept = estimate_mean ? mean_parameter[1] : fixed_mean;
  theta[1 : 2] = covariance_parameters[1 : 2];
  theta[3] = has_group ? covariance_parameters[3] : 1;
  theta[4] = has_group ? covariance_parameters[4] : 1;
  theta[5] = covariance_parameters[3 + 2 * has_group];
}
model {
  vector[6 + n_censored] state;
  state[1 : 5] = theta;
  state[6] = intercept;
  if (n_censored)
    state[7 : ] = quantile ? unit : censored;
  if (estimate_mean)
    intercept ~ hyperprior(prior_family[1], prior_arguments[1]', 0);
  for (k in 1 : num_elements(covariance_parameters)) {
    int prior_index = has_group || k < 3 ? k + 1 : 6;
    covariance_parameters[k] ~ hyperprior(prior_family[prior_index],
                                          prior_arguments[prior_index]', 1);
  }
  // No truncated-normal normalizer belongs in this augmented likelihood.
  // The quantile implementation instead includes its exact change of variables.
  target += cgp_log_density(state, distance2, group_cov, observed, lower_bound,
                            upper_bound, kernel, quantile);
}
generated quantities {
  vector[n_censored] y_cens;
  if (quantile && n_censored) {
    vector[6 + n_censored] state;
    state[1 : 5] = theta;
    state[6] = intercept;
    state[7 : ] = unit;
    y_cens = cgp_imputed_values(state, distance2, group_cov, observed, lower_bound,
               upper_bound, kernel);
  } else {
    y_cens = censored;
  }
}
