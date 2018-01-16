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
// Equality and inequality constraints on \beta can be incoporated as well.
class GeneralizedAdditiveModel {
public:
  // Start an unconstrained generalized additive model calculation
  // X: design matrix (NxK)
  // D: difference matrix (must have K columns)
  // sigma: prior standard deviation on lambda
  GeneralizedAdditiveModel(const Eigen::SparseMatrix<double>& X,
                           const Eigen::SparseMatrix<double>& D,
                           double sigma);
  
  // Add n_eq equality constraints on \beta by passing a
  // n_eq x K matrix Ceq such that the constraint is Ceq \beta = 0
  //TODO: cannot use Quadprog++ only for equality constraints
  //void set_equality_constraints(const Eigen::MatrixXd& Ceq) { Ceq_ = Ceq; neq_ = Ceq_.rows(); }
  
  // Add n_in inequality constraints on \beta by passing a
  // n_in x K matrix Cin such that the constraint is Cin \beta >= 0
  void set_inequality_constraints(const Eigen::MatrixXd& Cin) { Cin_ = Cin; nin_ = Cin_.rows(); }
  
  // Fit GAM to y and Sm1. Cholesky decomposition is only computed the first time.
  // y: observations (vector of size N)
  // Sm1: inverse standard deviations for each observation (vector of size N)
  // max_iter: maximum number of iterations
  // tol_val: relative tolerance on the final values of X\beta
  void optimize(const Eigen::VectorXd& y, const Eigen::VectorXd& Sm1, unsigned max_iter, double tol_val);

  //accessors
  Eigen::VectorXd get_beta() const { return beta_; }
  Eigen::VectorXd get_mean() const { return X_ * beta_; }
  double get_lambda() const { return lambda_; }
  bool has_converged() const { return has_converged_; }
  
private:
  
  Eigen::SparseMatrix<double> get_L(const Eigen::SparseMatrix<double>& A);
    
  //data
  unsigned N_, K_;
  Eigen::SparseMatrix<double> X_, D_;
  double sigma_;
  unsigned neq_, nin_;
  Eigen::MatrixXd Ceq_, Cin_;
  
  //params
  Eigen::VectorXd beta_;
  double lambda_;
  bool has_converged_;
  
  //helpers
  bool llt_analyzed_;
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower,
                       Eigen::NaturalOrdering<Eigen::SparseMatrix<double>::StorageIndex> > llt_; //no permutation
  
};

Eigen::SparseMatrix<double> first_order_difference_matrix(unsigned K);
Eigen::SparseMatrix<double> second_order_difference_matrix(unsigned K);

}
  
#endif