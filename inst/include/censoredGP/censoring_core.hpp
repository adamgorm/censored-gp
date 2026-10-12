#ifndef CENSORED_GP_CENSORING_CORE_HPP
#define CENSORED_GP_CENSORING_CORE_HPP

#include "gaussian_core.hpp"
#include <limits>

namespace censored_gp {

// Stable log(exp(a) + exp(b)), also when one argument is negative infinity.
inline double log_add(double a, double b) {
  if (a < b) std::swap(a, b);
  if (!std::isfinite(b)) return a;
  return a + std::log1p(std::exp(b - a));
}

inline double log_subtract(double a, double b) {
  if (!(a > b))
    throw std::domain_error("Censoring interval has zero numerical probability");
  return a + std::log(-std::expm1(b - a));
}

// Contract a symmetric covariance adjoint with the five kernel derivatives.
// Only the lower triangle is visited; off-diagonal entries count twice.
inline Vector covariance_score(const Matrix& distance2, const Matrix& group_cov,
                               const Vector& theta, int kernel,
                               const Matrix& adjoint) {
  Vector score = Vector::Zero(5);
  for (Eigen::Index column = 0; column < distance2.cols(); ++column) {
    for (Eigen::Index row = column; row < distance2.rows(); ++row) {
      const KernelValue population = kernel_value(
        distance2(row, column), theta[0], theta[1], kernel);
      const KernelValue group = kernel_value(
        distance2(row, column), theta[2], theta[3], kernel);
      const double weight = adjoint(row, column) * (row == column ? 1.0 : 2.0);
      score[0] += weight * 2.0 * population.value / theta[0];
      score[1] += weight * population.length_derivative;
      score[2] += weight * 2.0 * group_cov(row, column) * group.value / theta[2];
      score[3] += weight * group_cov(row, column) * group.length_derivative;
    }
  }
  score[4] = 2.0 * theta[4] * adjoint.trace();
  return score;
}

struct CensoringEvaluation {
  double log_density;
  Vector censored_values;
  Vector gradient_theta;
  Vector gradient_unit;
  double gradient_mean = 0.0;

  // NormalFunctions supplies log_cdf(x) and log_quantile(log_p). The R adapter
  // uses Rmath; the Stan adapter uses Stan Math. Neither integrates an MVN CDF.
  template <class NormalFunctions>
  CensoringEvaluation(const Matrix& distance2, const Matrix& group_cov,
                      const Vector& theta, double mean, const Vector& observed,
                      const Vector& unit, const Vector& lower, const Vector& upper,
                      int kernel, const NormalFunctions& normal, bool gradients) {
    const Eigen::Index n_observed = observed.size();
    const Eigen::Index n_censored = unit.size();
    const Eigen::Index n = n_observed + n_censored;
    if (lower.size() != n_censored || upper.size() != n_censored || n < 1)
      throw std::invalid_argument("Invalid censoring dimensions");
    if (!std::isfinite(mean) || !unit.allFinite()
        || (n_censored && (unit.minCoeff() <= 0.0 || unit.maxCoeff() >= 1.0)))
      throw std::domain_error("Quantile parameters must be strictly between zero and one");

    // Reuse the same checked covariance builder and Cholesky implementation as
    // the Gaussian likelihood. This temporary residual is not part of a target.
    const GaussianEvaluation gaussian(distance2, group_cov, theta,
      Vector::Zero(n), kernel, false);
    const Matrix L = gaussian.factor.matrixL();
    Vector z = Vector::Zero(n);
    if (n_observed) {
      const Vector residual = observed.array() - mean;
      z.head(n_observed) = L.topLeftCorner(n_observed, n_observed)
        .triangularView<Eigen::Lower>().solve(residual);
    }
    const double log_two_pi = 1.83787706640934548356;
    log_density = -0.5 * (n_observed * log_two_pi + z.head(n_observed).squaredNorm());
    for (Eigen::Index i = 0; i < n_observed; ++i) log_density -= std::log(L(i, i));

    Vector a(n_censored), b(n_censored), log_mass(n_censored);
    censored_values.resize(n_censored);
    for (Eigen::Index k = 0; k < n_censored; ++k) {
      const Eigen::Index i = n_observed + k;
      if (!(lower[k] < upper[k]))
        throw std::domain_error("Censoring lower bound must be below upper bound");
      const double conditional_mean = i ? L.row(i).head(i).dot(z.head(i)) : 0.0;
      a[k] = (lower[k] - mean - conditional_mean) / L(i, i);
      b[k] = (upper[k] - mean - conditional_mean) / L(i, i);
      const double log_fa = normal.log_cdf(a[k]);
      const double log_fb = normal.log_cdf(b[k]);
      const double log_sa = normal.log_cdf(-a[k]);
      const double log_sb = normal.log_cdf(-b[k]);
      // Reflection avoids subtracting two CDFs rounded to one in the right tail.
      log_mass[k] = a[k] > 0 ? log_subtract(log_sa, log_sb)
                                    : log_subtract(log_fb, log_fa);
      const double log_p = log_add(log_fa, std::log(unit[k]) + log_mass[k]);
      z[i] = log_p < -0.69314718055994530942
        ? normal.log_quantile(log_p)
        : -normal.log_quantile(log_add(log_sb, std::log1p(-unit[k]) + log_mass[k]));
      censored_values[k] = mean + conditional_mean + L(i, i) * z[i];
      log_density += log_mass[k];
    }
    if (!std::isfinite(log_density) || !z.allFinite())
      throw std::domain_error("Censoring transformation is not numerically finite");
    if (!gradients) return;

    // Reverse the triangular transformation explicitly. This keeps O(N^2)
    // matrix operations off the scalar autodiff tape while preserving every
    // dependency of the sequential quantiles and their Jacobian factors.
    Matrix adjoint_L = Matrix::Zero(n, n);
    Vector adjoint_z = Vector::Zero(n);
    adjoint_z.head(n_observed) = -z.head(n_observed);
    gradient_unit = Vector::Zero(n_censored);
    for (Eigen::Index i = 0; i < n_observed; ++i) adjoint_L(i, i) = -1.0 / L(i, i);

    for (Eigen::Index k = n_censored; k-- > 0;) {
      const Eigen::Index i = n_observed + k;
      const double log_phi_z = -0.5 * (z[i] * z[i] + log_two_pi);
      double adjoint_a = 0.0, adjoint_b = 0.0;
      if (std::isfinite(a[k])) {
        const double log_phi_a = -0.5 * (a[k] * a[k] + log_two_pi);
        adjoint_a = -std::exp(log_phi_a - log_mass[k]);
        if (adjoint_z[i] != 0.0)
          adjoint_a += adjoint_z[i] * std::exp(std::log1p(-unit[k]) + log_phi_a - log_phi_z);
      }
      if (std::isfinite(b[k])) {
        const double log_phi_b = -0.5 * (b[k] * b[k] + log_two_pi);
        adjoint_b = std::exp(log_phi_b - log_mass[k]);
        if (adjoint_z[i] != 0.0)
          adjoint_b += adjoint_z[i] * std::exp(std::log(unit[k]) + log_phi_b - log_phi_z);
      }
      if (adjoint_z[i] != 0.0)
        gradient_unit[k] = adjoint_z[i] * std::exp(log_mass[k] - log_phi_z);
      const double adjoint_conditional_mean = -(adjoint_a + adjoint_b) / L(i, i);
      gradient_mean += adjoint_conditional_mean;
      if (std::isfinite(a[k])) adjoint_L(i, i) -= a[k] * adjoint_a / L(i, i);
      if (std::isfinite(b[k])) adjoint_L(i, i) -= b[k] * adjoint_b / L(i, i);
      if (i) {
        adjoint_L.row(i).head(i) += adjoint_conditional_mean * z.head(i).transpose();
        adjoint_z.head(i) += adjoint_conditional_mean * L.row(i).head(i).transpose();
      }
    }
    if (n_observed) {
      const Vector adjoint_residual = L.topLeftCorner(n_observed, n_observed).transpose()
        .triangularView<Eigen::Upper>().solve(adjoint_z.head(n_observed));
      gradient_mean -= adjoint_residual.sum();
      for (Eigen::Index column = 0; column < n_observed; ++column)
        for (Eigen::Index row = column; row < n_observed; ++row)
          adjoint_L(row, column) -= adjoint_residual[row] * z[column];
    }

    // Reverse-mode Cholesky identity: project L' * adjoint_L onto its lower
    // triangle, halve its diagonal, then solve on both sides and symmetrize.
    Matrix projected = L.transpose() * adjoint_L;
    projected.triangularView<Eigen::StrictlyUpper>().setZero();
    projected.diagonal() *= 0.5;
    const Matrix left = L.transpose().triangularView<Eigen::Upper>().solve(projected);
    const Matrix unsymmetric = L.transpose().triangularView<Eigen::Upper>()
      .solve(left.transpose()).transpose();
    const Matrix adjoint_covariance = 0.5 * (unsymmetric + unsymmetric.transpose());
    gradient_theta = covariance_score(distance2, group_cov, theta, kernel, adjoint_covariance);
    if (!gradient_theta.allFinite() || !gradient_unit.allFinite() || !std::isfinite(gradient_mean))
      throw std::domain_error("Censoring derivatives are not numerically finite");
  }
};
}  // namespace censored_gp
#endif
