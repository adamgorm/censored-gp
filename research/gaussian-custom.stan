// Same likelihood and state as gaussian-reference.stan, with a custom adjoint.
functions {
  real cgp_gaussian_log_density(vector state, data matrix distance2,
                                data matrix group_cov, data vector observed,
                                data int kernel);
}
data {
  int<lower=1> N;
  int<lower=0, upper=N> No;
  int<lower=1, upper=4> kernel;
  matrix[N, N] distance2;
  matrix[N, N] group_cov;
  vector[No] observed;
}
parameters {
  vector[5 + N - No] q;
}
model {
  vector[5 + N - No] state = append_row(exp(head(q, 5)), tail(q, N - No));
  target += cgp_gaussian_log_density(state, distance2, group_cov, observed, kernel);
}
