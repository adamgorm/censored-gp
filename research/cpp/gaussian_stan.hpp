#ifndef CENSORED_GP_RESEARCH_GAUSSIAN_STAN_HPP
#define CENSORED_GP_RESEARCH_GAUSSIAN_STAN_HPP
#include <stan/model/model_header.hpp>
#include "gaussian_core.hpp"

// One vector of parameters makes the analytical adjoint easy to audit.
// state = five positive covariance parameters followed by censored values.
template <typename V, stan::require_eigen_vector_t<V>* = nullptr>
inline stan::value_type_t<V> cgp_gaussian_log_density(
    const V& state, const Eigen::MatrixXd& distance2,
    const Eigen::MatrixXd& group_cov, const Eigen::VectorXd& observed,
    int kernel, std::ostream*) {
  const Eigen::VectorXd values = stan::math::value_of(state);
  if (values.size() < 5 || distance2.rows() != observed.size() + values.size() - 5)
    throw std::invalid_argument("Invalid observed/augmented state size");
  Eigen::VectorXd residual(distance2.rows());
  residual.head(observed.size()) = observed;
  residual.tail(values.size() - 5) = values.tail(values.size() - 5);
  constexpr bool gradients = stan::is_var<stan::value_type_t<V>>::value;
  const censored_gp::GaussianEvaluation result(distance2, group_cov,
    values.head(5), residual, kernel, gradients);
  if constexpr (gradients) {
    std::vector<stan::math::var> operands;
    std::vector<double> partials;
    operands.reserve(values.size());
    partials.reserve(values.size());
    for (Eigen::Index i = 0; i < values.size(); ++i) {
      operands.push_back(state.coeff(i));
      partials.push_back(i < 5 ? result.gradient_theta[i]
        : result.gradient_y[observed.size() + i - 5]);
    }
    return stan::math::precomputed_gradients(result.log_density, operands, partials);
  } else {
    return result.log_density;
  }
}
#endif
