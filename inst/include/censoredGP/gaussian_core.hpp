#ifndef CENSORED_GP_RESEARCH_GAUSSIAN_CORE_HPP
#define CENSORED_GP_RESEARCH_GAUSSIAN_CORE_HPP

// Exact dense Gaussian calculations. This header has no R or Stan dependency.
// Censored measurements enter as augmented values; their bounds belong in the
// sampler. This is not a numerically integrated Gaussian rectangle probability.
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

namespace censored_gp {
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

// 1: squared exponential; 2: exponential; 3/4: Matern 3/2 and 5/2.
struct KernelValue {
  double value;
  double length_derivative;
};

inline KernelValue kernel_value(double distance2, double amplitude,
                                double length, int kernel) {
  const double a2 = amplitude * amplitude;
  if (kernel == 1) {
    const double r2 = distance2 / (length * length);
    const double value = a2 * std::exp(-0.5 * r2);
    return {value, value * r2 / length};
  }
  double r = std::sqrt(distance2) / length;
  if (kernel == 2) {
    const double value = a2 * std::exp(-r);
    return {value, value * r / length};
  }
  r *= std::sqrt(kernel == 3 ? 3.0 : 5.0);
  const double e = a2 * std::exp(-r);
  if (kernel == 3) return {e * (1.0 + r), e * r * r / length};
  return {e * (1.0 + r + r * r / 3.0),
          e * (r * r + r * r * r) / (3.0 * length)};
}

struct GaussianEvaluation {
  Eigen::LLT<Matrix> factor;
  Vector alpha;
  double log_density;
  Vector gradient_theta;
  Vector gradient_y;

  // theta = population amplitude/length, group amplitude/length, noise SD.
  // group_cov may be zero for a single GP, or any supplied valid group
  // covariance. No repeated inputs or regular grids are assumed.
  GaussianEvaluation(const Matrix& distance2, const Matrix& group_cov,
                     const Vector& theta, const Vector& residual,
                     int kernel, bool gradients = true) {
    const Eigen::Index n = residual.size();
    if (n < 1 || distance2.rows() != n || distance2.cols() != n
        || group_cov.rows() != n || group_cov.cols() != n
        || theta.size() != 5 || kernel < 1 || kernel > 4)
      throw std::invalid_argument("Invalid Gaussian input dimensions or kernel");
    if (!distance2.allFinite() || distance2.minCoeff() < 0
        || !group_cov.allFinite() || !residual.allFinite()
        || !theta.allFinite() || theta.minCoeff() <= 0)
      // Stan treats domain errors as rejected proposals. An invalid_argument
      // here would abort a chain when warmup proposes an overflowing scale.
      throw std::domain_error("Inputs must be finite and scales positive");
    if (!distance2.isApprox(distance2.transpose(), 1e-12)
        || !group_cov.isApprox(group_cov.transpose(), 1e-12))
      throw std::invalid_argument("Distance and group matrices must be symmetric");
    Matrix covariance(n, n);
    for (Eigen::Index j = 0; j < n; ++j) {
      for (Eigen::Index i = j; i < n; ++i) {
        const double v = kernel_value(distance2(i,j), theta[0], theta[1], kernel).value
          + group_cov(i,j) * kernel_value(distance2(i,j), theta[2], theta[3], kernel).value;
        covariance(i,j) = covariance(j,i) = v;
      }
      covariance(j,j) += theta[4] * theta[4];
    }
    if (!covariance.allFinite())
      throw std::domain_error("Covariance entries are not finite");
    factor.compute(covariance);
    if (factor.info() != Eigen::Success)
      throw std::domain_error("Covariance is not positive definite");
    alpha = factor.solve(residual);
    const double log_determinant_half = factor.matrixLLT().diagonal().array().log().sum();
    const double log_two_pi = 1.83787706640934548356;
    log_density = -0.5 * (n * log_two_pi + residual.dot(alpha)) - log_determinant_half;
    if (!gradients) return;

    // W = alpha alpha' - K^-1; solve using the existing Cholesky factor.
    Matrix W = -factor.solve(Matrix::Identity(n, n));
    W.noalias() += alpha * alpha.transpose();
    gradient_y = -alpha;
    gradient_theta = Vector::Zero(5);
    for (Eigen::Index j = 0; j < n; ++j) {
      for (Eigen::Index i = j; i < n; ++i) {
        const KernelValue mu = kernel_value(distance2(i,j), theta[0], theta[1], kernel);
        const KernelValue eta = kernel_value(distance2(i,j), theta[2], theta[3], kernel);
        const double w = W(i,j) * (i == j ? 0.5 : 1.0);
        gradient_theta[0] += w * 2.0 * mu.value / theta[0];
        gradient_theta[1] += w * mu.length_derivative;
        gradient_theta[2] += w * 2.0 * group_cov(i,j) * eta.value / theta[2];
        gradient_theta[3] += w * group_cov(i,j) * eta.length_derivative;
      }
    }
    gradient_theta[4] = theta[4] * W.trace();
    if (!gradient_theta.allFinite() || !gradient_y.allFinite())
      throw std::domain_error("Gaussian gradients are not finite");
  }

  void check_cross(const Matrix& cross) const {
    if (cross.cols() != alpha.size() || !cross.allFinite())
      throw std::invalid_argument("Invalid prediction cross covariance");
  }

  Vector predict_mean(const Matrix& cross) const {
    check_cross(cross);
    return cross * alpha;
  }

  Vector predict_variance(const Matrix& cross, const Vector& prior_variance) const {
    check_cross(cross);
    if (prior_variance.size() != cross.rows() || !prior_variance.allFinite())
      throw std::invalid_argument("Invalid prediction variances");
    const Matrix v = factor.matrixL().solve(cross.transpose());
    return prior_variance - v.array().square().colwise().sum().transpose().matrix();
  }

  Matrix predict_covariance(const Matrix& cross, const Matrix& prior) const {
    check_cross(cross);
    if (prior.rows() != cross.rows() || prior.cols() != cross.rows()
        || !prior.allFinite() || !prior.isApprox(prior.transpose(), 1e-12))
      throw std::invalid_argument("Invalid prediction covariance");
    const Matrix v = factor.matrixL().solve(cross.transpose());
    Matrix result = prior;
    result.noalias() -= v.transpose() * v;
    return result;
  }
};
}  // namespace censored_gp
#endif
