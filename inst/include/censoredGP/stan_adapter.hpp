#ifndef CENSORED_GP_STAN_ADAPTER_HPP
#define CENSORED_GP_STAN_ADAPTER_HPP

#include <stan/model/model_header.hpp>
#include "censoring_core.hpp"

struct CgpStanNormal {
  double log_cdf(double value) const { return stan::math::std_normal_lcdf(value); }
  double log_quantile(double log_probability) const {
    return stan::math::std_normal_log_qf(log_probability);
  }
};

// State order: five covariance parameters, constant mean, then either censored
// measurements or unit-interval quantiles. Only the state requires derivatives.
template <typename V, stan::require_eigen_vector_t<V>* = nullptr>
inline stan::value_type_t<V> cgp_log_density(
    const V& state, const Eigen::MatrixXd& distance2,
    const Eigen::MatrixXd& group_cov, const Eigen::VectorXd& observed,
    const Eigen::VectorXd& lower, const Eigen::VectorXd& upper,
    int kernel, int quantile, std::ostream*) {
  const Eigen::VectorXd values = stan::math::value_of(state);
  const Eigen::Index n_censored = values.size() - 6;
  if (n_censored < 0 || observed.size() + n_censored != distance2.rows())
    throw std::invalid_argument("Invalid augmented state dimensions");
  constexpr bool gradients = stan::is_var<stan::value_type_t<V>>::value;
  Eigen::VectorXd partials(values.size());
  double log_density;
  if (quantile && n_censored) {
    const censored_gp::CensoringEvaluation result(distance2, group_cov,
      values.head(5), values[5], observed, values.tail(n_censored), lower, upper,
      kernel, CgpStanNormal(), gradients);
    log_density = result.log_density;
    if constexpr (gradients) {
      partials.head(5) = result.gradient_theta;
      partials[5] = result.gradient_mean;
      partials.tail(n_censored) = result.gradient_unit;
    }
  } else {
    Eigen::VectorXd residual(distance2.rows());
    residual.head(observed.size()) = observed.array() - values[5];
    residual.tail(n_censored) = values.tail(n_censored).array() - values[5];
    const censored_gp::GaussianEvaluation result(distance2, group_cov,
      values.head(5), residual, kernel, gradients);
    log_density = result.log_density;
    if constexpr (gradients) {
      partials.head(5) = result.gradient_theta;
      partials[5] = -result.gradient_y.sum();
      partials.tail(n_censored) = result.gradient_y.tail(n_censored);
    }
  }
  if constexpr (gradients) {
    std::vector<stan::math::var> operands;
    std::vector<double> scores;
    operands.reserve(values.size());
    scores.reserve(values.size());
    for (Eigen::Index index = 0; index < values.size(); ++index) {
      operands.push_back(state.coeff(index));
      scores.push_back(partials[index]);
    }
    return stan::math::precomputed_gradients(log_density, operands, scores);
  } else {
    return log_density;
  }
}

// Called only in generated quantities, once per saved posterior draw. Computing
// imputed measurements here avoids a second Cholesky factorization per leapfrog.
inline Eigen::VectorXd cgp_imputed_values(
    const Eigen::VectorXd& state, const Eigen::MatrixXd& distance2,
    const Eigen::MatrixXd& group_cov, const Eigen::VectorXd& observed,
    const Eigen::VectorXd& lower, const Eigen::VectorXd& upper,
    int kernel, std::ostream*) {
  const censored_gp::CensoringEvaluation result(distance2, group_cov,
    state.head(5), state[5], observed, state.tail(state.size() - 6), lower, upper,
    kernel, CgpStanNormal(), false);
  return result.censored_values;
}
#endif
