#ifndef ANALYTICAL_GAM_LIBRARY_HPP
#define ANALYTICAL_GAM_LIBRARY_HPP

#include <Eigen/Eigen>

namespace binless {

//policy class that implements a generalized additive model using the normal equations
//does not support inequality constraints but scales to very large datasets
class AnalyticalGAMLibrary {
public:
  // Start an unconstrained generalized additive model calculation
  // X: design matrix (NxK)
  // D: difference matrix (must have K columns)
  // sigma: prior standard deviation on lambda
  AnalyticalGAMLibrary(const Eigen::SparseMatrix<double>& X,
                       const Eigen::SparseMatrix<double>& D,
                       double sigma);
  
  // Add n_eq equality constraints on \beta by passing a
  // n_eq x K matrix Ceq such that the constraint is Ceq \beta = 0
  void set_equality_constraints(const Eigen::SparseMatrix<double, Eigen::RowMajor>& Ceq) { Ceq_ = Ceq; neq_ = Ceq_.rows(); }
  
  //this call must be here, see GAMFitterImpl constructor
  void set_inequality_constraints(const Eigen::SparseMatrix<double, Eigen::RowMajor>& Cin) { Rcpp::stop("Not implemented"); }
  
protected:
  void optimize(const Eigen::VectorXd& y, const Eigen::VectorXd& Sm1, unsigned max_iter, double tol_val);
  
  //accessors
  Eigen::VectorXd get_beta() const { return beta_; }
  Eigen::VectorXd get_mean() const { return X_ * beta_; }
  double get_lambda() const { return lambda_; }
  bool has_converged() const { return has_converged_; }
  
  //to avoid direct destruction by user
  ~AnalyticalGAMLibrary() {}
  
private:
  
  void update_LDLT_A(const Eigen::SparseMatrix<double>& A);
  void update_LDLT_CAC(const Eigen::SparseMatrix<double>& A);
    
  //data
  unsigned N_, K_;
  Eigen::SparseMatrix<double> X_, D_;
  double sigmasqinv_;
  unsigned neq_;
  Eigen::SparseMatrix<double, Eigen::RowMajor> Ceq_;
  
  //params
  Eigen::VectorXd beta_;
  double lambda_;
  bool has_converged_;
  
  //helpers
  bool ldlt_A_analyzed_, ldlt_CAC_analyzed_;
  typedef Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver_t;
  solver_t ldlt_A_, ldlt_CAC_;
  
};

}
  
#endif