// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppEigen.h>
#include "gaussian_core.hpp"

// [[Rcpp::export]]
Rcpp::List gaussian_evaluate_cpp(const Eigen::MatrixXd& distance2,
                                const Eigen::MatrixXd& group_cov,
                                const Eigen::VectorXd& theta,
                                const Eigen::VectorXd& residual,
                                int kernel = 1, bool gradients = true) {
  Rcpp::XPtr<censored_gp::GaussianEvaluation> fit(
    new censored_gp::GaussianEvaluation(distance2, group_cov, theta, residual,
                                        kernel, gradients), true);
  return Rcpp::List::create(Rcpp::Named("log_density") = fit->log_density,
                            Rcpp::Named("gradient_theta") = fit->gradient_theta,
                            Rcpp::Named("gradient_y") = fit->gradient_y,
                            Rcpp::Named("factor") = fit);
}

// [[Rcpp::export]]
Rcpp::List gaussian_predict_marginal_cpp(SEXP factor, const Eigen::MatrixXd& cross,
                                        const Eigen::VectorXd& prior_variance) {
  Rcpp::XPtr<censored_gp::GaussianEvaluation> fit(factor);
  return Rcpp::List::create(Rcpp::Named("mean") = fit->predict_mean(cross),
    Rcpp::Named("variance") = fit->predict_variance(cross, prior_variance));
}

// [[Rcpp::export]]
Rcpp::List gaussian_predict_joint_cpp(SEXP factor, const Eigen::MatrixXd& cross,
                                     const Eigen::MatrixXd& prior) {
  Rcpp::XPtr<censored_gp::GaussianEvaluation> fit(factor);
  return Rcpp::List::create(Rcpp::Named("mean") = fit->predict_mean(cross),
    Rcpp::Named("covariance") = fit->predict_covariance(cross, prior));
}
