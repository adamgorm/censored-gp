// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppEigen.h>
#include <censoredGP/censoring_core.hpp>

// These wrappers keep R memory management outside the reusable numerical core.
// A factor stores one covariance and one augmented response vector; prediction
// reuses it across all requested components at that posterior draw.

// [[Rcpp::export]]
Rcpp::List cgp_evaluate(const Eigen::MatrixXd& distance2,
                        const Eigen::MatrixXd& group_cov,
                        const Eigen::VectorXd& theta,
                        const Eigen::VectorXd& residual,
                        int kernel, bool gradients) {
  Rcpp::XPtr<censored_gp::GaussianEvaluation> factor(
    new censored_gp::GaussianEvaluation(distance2, group_cov, theta, residual,
                                       kernel, gradients), true);
  return Rcpp::List::create(Rcpp::Named("log_density") = factor->log_density,
    Rcpp::Named("gradient_theta") = factor->gradient_theta,
    Rcpp::Named("gradient_y") = factor->gradient_y,
    Rcpp::Named("factor") = factor);
}

// [[Rcpp::export]]
Rcpp::List cgp_condition(SEXP factor_pointer, const Eigen::MatrixXd& cross,
                         const Eigen::VectorXd& prior_variance,
                         const Eigen::MatrixXd& prior, bool joint) {
  Rcpp::XPtr<censored_gp::GaussianEvaluation> factor(factor_pointer);
  if (joint) return Rcpp::List::create(Rcpp::Named("mean") = factor->predict_mean(cross),
    Rcpp::Named("covariance") = factor->predict_covariance(cross, prior));
  return Rcpp::List::create(Rcpp::Named("mean") = factor->predict_mean(cross),
    Rcpp::Named("variance") = factor->predict_variance(cross, prior_variance));
}

struct RNormalFunctions {
  double log_cdf(double value) const { return R::pnorm(value, 0.0, 1.0, true, true); }
  double log_quantile(double value) const { return R::qnorm(value, 0.0, 1.0, true, true); }
};

// Used by the independent numerical validation suite as well as package tests.
// [[Rcpp::export]]
Rcpp::List cgp_censor_evaluate(const Eigen::MatrixXd& distance2,
                               const Eigen::MatrixXd& group_cov,
                               const Eigen::VectorXd& theta, double mean,
                               const Eigen::VectorXd& observed,
                               const Eigen::VectorXd& unit,
                               const Eigen::VectorXd& lower,
                               const Eigen::VectorXd& upper, int kernel) {
  const censored_gp::CensoringEvaluation result(distance2, group_cov, theta,
    mean, observed, unit, lower, upper, kernel, RNormalFunctions(), true);
  return Rcpp::List::create(Rcpp::Named("log_density") = result.log_density,
    Rcpp::Named("censored_values") = result.censored_values,
    Rcpp::Named("gradient_theta") = result.gradient_theta,
    Rcpp::Named("gradient_mean") = result.gradient_mean,
    Rcpp::Named("gradient_unit") = result.gradient_unit);
}
