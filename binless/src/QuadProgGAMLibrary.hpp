#ifndef QUADPROG_GAM_LIBRARY_HPP
#define QUADPROG_GAM_LIBRARY_HPP

#include <Eigen/Eigen>

namespace binless {

//policy class that implements a generalized additive model using the QuadProg++ library
class QuadProgGAMLibrary {
public:
  // Start an unconstrained generalized additive model calculation
  // X: design matrix (NxK)
  // D: difference matrix (must have K columns)
  // sigma: prior standard deviation on lambda
  QuadProgGAMLibrary(const Eigen::SparseMatrix<double>& X,
                     const Eigen::SparseMatrix<double>& D,
                     double sigma);
  
  // Add n_eq equality constraints on \beta by passing a
  // n_eq x K matrix Ceq such that the constraint is Ceq \beta = 0
  //TODO: cannot use Quadprog++ only for equality constraints
  //void set_equality_constraints(const Eigen::MatrixXd& Ceq) { Ceq_ = Ceq; neq_ = Ceq_.rows(); }
  
  // Add n_in inequality constraints on \beta by passing a
  // n_in x K matrix Cin such that the constraint is Cin \beta >= 0
  void set_inequality_constraints(const Eigen::MatrixXd& Cin) { Cin_ = Cin; nin_ = Cin_.rows(); }
  
  
protected:
  void optimize(const Eigen::VectorXd& y, const Eigen::VectorXd& Sm1, unsigned max_iter, double tol_val);
  
  //accessors
  Eigen::VectorXd get_beta() const { return beta_; }
  Eigen::VectorXd get_mean() const { return X_ * beta_; }
  double get_lambda() const { return lambda_; }
  bool has_converged() const { return has_converged_; }
  
  //to avoid direct destruction by user
  ~QuadProgGAMLibrary() {}
  
private:
  
  Eigen::SparseMatrix<double> get_L(const Eigen::SparseMatrix<double>& A);
    
  //data
  unsigned N_, K_;
  Eigen::SparseMatrix<double> X_, D_;
  double sigmasq_;
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

}
  
#endif