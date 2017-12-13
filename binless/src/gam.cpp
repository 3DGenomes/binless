#include "RcppEigen.h"
#include "gam.hpp"
#include "QuadProg++.hpp"

namespace binless {

GeneralizedAdditiveModel::GeneralizedAdditiveModel(const Eigen::VectorXd& y,
                           const Eigen::VectorXd& S,
                           const Eigen::SparseMatrix<double>& X,
                           const Eigen::SparseMatrix<double>& D,
                           double sigma) :
  N_(X.rows()), K_(X.cols()), y_(y), S_(S), X_(X), D_(D), sigma_(sigma), neq_(0), nin_(0),
  Ceq_(Eigen::MatrixXd::Zero(neq_,K_)), Cin_(Eigen::MatrixXd::Zero(nin_,K_)),
  beta_(Eigen::VectorXd::Constant(K_,1)), lambda_(1), has_converged_(false) {}
  
void GeneralizedAdditiveModel::optimize(unsigned max_iter, double tol_val) {
  //build fixed matrices
  const auto XtSm1                         = (S_.asDiagonal().inverse()*X_).transpose();
  const Eigen::SparseMatrix<double> XtSm2X = Eigen::SparseMatrix<double>(K_,K_)
                                             .selfadjointView<Eigen::Lower>().rankUpdate(XtSm1);
  const auto XtSm2y                        = XtSm1 * (S_.asDiagonal().inverse() * y_);
  const Eigen::SparseMatrix<double> K2DtD  = Eigen::SparseMatrix<double>(K_,K_)
                                             .selfadjointView<Eigen::Lower>().rankUpdate(K_*D_.transpose());
  //prepare loop
  Eigen::VectorXd beta_old = beta_;
  unsigned iter = 0;
  has_converged_ = false;
  do {
    //compute A
    const Eigen::SparseMatrix<double> A = XtSm2X + lambda_*lambda_*K2DtD;
    
    Rcpp::Rcout << "ready to go\n";
    //compute L factor (convert to dense for quadprog++)
    Eigen::MatrixXd L = get_L(A);
   
    //compute new beta. In QuadProg++, G = A and g0 = -XtSm2y
    const Eigen::VectorXd g0 = -XtSm2y;
    
    // Find the unconstrained minimizer of the quadratic form 0.5 * x G x + g0 x
    // this is a feasible point in the dual space
    // x = G^-1 * g0
    beta_ = L.triangularView<Eigen::Lower>().solve(g0);
    beta_ = L.triangularView<Eigen::Lower>().adjoint().solve(beta_);
    beta_ = -beta_;
    
    if (nin_>0) {
      double score = 2*QuadProgPP::solve_quadprog_with_guess(L, g0, Ceq_.transpose(), Eigen::VectorXd::Zero(neq_),
                                                             Cin_.transpose(), Eigen::VectorXd::Zero(nin_), beta_);
    }
    
    //compute lambda
    lambda_ = std::sqrt((K_ - 2)/((K_*D_*beta_).squaredNorm()+1));
    
    //compute relative precision on beta and check convergence
    double precision = beta_.maxCoeff()-beta_.minCoeff();
    if ( precision>0 ) precision = (beta_ - beta_old).array().abs().maxCoeff() / precision;
    has_converged_ = precision < tol_val;
    //bookkeeping
    beta_old = beta_;
    iter++;
  } while (iter < max_iter && (!has_converged_));
}

Eigen::MatrixXd GeneralizedAdditiveModel::get_L(const Eigen::SparseMatrix<double>& A) {
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower,
                       Eigen::NaturalOrdering<Eigen::SparseMatrix<double>::StorageIndex> > llt; //no permutation
  Rcpp::Rcout << "info init: " << llt.info() << "\n";
  llt.analyzePattern(A);
  Rcpp::Rcout << "info analyze: " << llt.info() << "\n";
  llt.factorize(A);
  Rcpp::Rcout << "info factorize: " << llt.info() << "\n";
  if (llt.info() != Eigen::Success) Rcpp::stop("decomposition failed");
  return llt.matrixL();
}



Eigen::SparseMatrix<double> first_order_difference_matrix(unsigned K) {
  if (K<=1) return Eigen::SparseMatrix<double>();
  Eigen::MatrixXi D = - Eigen::MatrixXi::Identity(K-1,K);
  D.diagonal(1).setConstant(1);
  return D.sparseView().cast<double>();
}

Eigen::SparseMatrix<double> second_order_difference_matrix(unsigned K) {
  if (K<=2) return Eigen::SparseMatrix<double>();
  Eigen::MatrixXi D = Eigen::MatrixXi::Identity(K-2,K);
  D.diagonal(1).setConstant(-2);
  D.diagonal(2).setConstant(1);
  return D.sparseView().cast<double>();
}


} // namespace binless
