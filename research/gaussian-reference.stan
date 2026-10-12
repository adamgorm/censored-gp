// Likelihood/gradient benchmark, not a censored-data sampling model.
data {
  int<lower=1> N;
  int<lower=0,upper=N> No;
  int<lower=1,upper=4> kernel;
  matrix[N,N] distance2;
  matrix[N,N] group_cov;
  vector[No] observed;
}
transformed data {
  matrix[N,N] distance = sqrt(distance2);
}
parameters {
  vector[5 + N - No] q;
}
model {
  vector[5] theta = exp(head(q, 5));
  vector[N] y = append_row(observed, tail(q, N - No));
  matrix[N,N] K_mu;
  matrix[N,N] K_eta;
  if(kernel == 1) {
    K_mu = square(theta[1]) * exp(-0.5 * distance2 / square(theta[2]));
    K_eta = square(theta[3]) * exp(-0.5 * distance2 / square(theta[4]));
  } else {
    matrix[N,N] r_mu = distance / theta[2];
    matrix[N,N] r_eta = distance / theta[4];
    if(kernel == 2) {
      K_mu = square(theta[1]) * exp(-r_mu);
      K_eta = square(theta[3]) * exp(-r_eta);
    } else if(kernel == 3) {
      r_mu *= sqrt(3.0);
      r_eta *= sqrt(3.0);
      K_mu = square(theta[1]) * (1 + r_mu) .* exp(-r_mu);
      K_eta = square(theta[3]) * (1 + r_eta) .* exp(-r_eta);
    } else {
      r_mu *= sqrt(5.0);
      r_eta *= sqrt(5.0);
      K_mu = square(theta[1]) * (1 + r_mu + r_mu .* r_mu / 3) .* exp(-r_mu);
      K_eta = square(theta[3]) * (1 + r_eta + r_eta .* r_eta / 3) .* exp(-r_eta);
    }
  }
  target += multi_normal_cholesky_lpdf(y | rep_vector(0,N),
    cholesky_decompose(add_diag(K_mu + group_cov .* K_eta, square(theta[5]))));
}
