#ifndef GAM_HPP
#define GAM_HPP

#include <Eigen/Eigen>

namespace binless {

//class to fit Generalized Additive Models (Normal residuals, performance iteration)
//
// The score is
// | S^{-1}(y - X\beta) |^2 + K^2 \lambda^2 | D\beta |^2 - 2(K-2) \log \lambda + \lambda^2/\sigma^2
// which is minimized with respect to \beta and \lambda. |.| is the L2 norm, and the target
// corresponds to a normal likelihood of y with mean X\beta and diagonal variance matrix S^{-2},
// a normal prior on \beta with zero mean and variance D^{T}D, and a normal prior on \lambda
// with zero mean and variance \sigma^2.
// Equality and inequality constraints on \beta can be incoporated as well, depending on which Library is used
template<typename GAMLibrary>
class GeneralizedAdditiveModel : public GAMLibrary {
public:
  // Start an unconstrained generalized additive model calculation
  // X: design matrix (NxK)
  // D: difference matrix (must have K columns)
  // sigma: prior standard deviation on lambda
  GeneralizedAdditiveModel(const Eigen::SparseMatrix<double>& X,
                           const Eigen::SparseMatrix<double>& D,
                           double sigma) : GAMLibrary(X,D,sigma) {}
  
  // Fit GAM to y and Sm1. Cholesky decomposition is only computed the first time.
  // y: observations (vector of size N)
  // Sm1: inverse standard deviations for each observation (vector of size N)
  // max_iter: maximum number of iterations
  // tol_val: relative tolerance on the final values of X\beta
  void optimize(const Eigen::VectorXd& y, const Eigen::VectorXd& Sm1, unsigned max_iter, double tol_val) {
    GAMLibrary::optimize(y,Sm1,max_iter,tol_val);
  }

  //accessors
  Eigen::VectorXd get_beta() const { return GAMLibrary::get_beta(); }
  Eigen::VectorXd get_mean() const { return GAMLibrary::get_mean(); }
  double get_lambda() const { return GAMLibrary::get_lambda(); }
  bool has_converged() const { return GAMLibrary::has_converged(); }
  
private:
  
};

Eigen::SparseMatrix<double> decreasing_constraint(unsigned K, unsigned Kfree);
Eigen::SparseMatrix<double> second_order_difference_matrix(unsigned K);

}
  
#endif