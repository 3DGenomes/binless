#include "RcppEigen.h"
#include "AnalyticalGAMLibrary.hpp"

namespace binless {

AnalyticalGAMLibrary::AnalyticalGAMLibrary(const Eigen::SparseMatrix<double>& X,
                           const Eigen::SparseMatrix<double>& D,
                           double sigma) :
  N_(X.rows()), K_(X.cols()), X_(X), D_(D), sigmasqinv_(1./(sigma*sigma)), neq_(0),
  Ceq_(Eigen::SparseMatrix<double>(neq_,K_)),
  beta_(Eigen::VectorXd::Constant(K_,1)), lambda_(1), has_converged_(false), ldlt_A_analyzed_(false),ldlt_CAC_analyzed_(false) {}
  
void AnalyticalGAMLibrary::update_LDLT_A(const Eigen::SparseMatrix<double>& A) {
  if (!ldlt_A_analyzed_) {
    ldlt_A_.analyzePattern(A); // sparsity pattern will not change if X and D are fixed
    ldlt_A_analyzed_ = true;
  }
  ldlt_A_.factorize(A);
  if (ldlt_A_.info() != Eigen::Success) Rcpp::stop("LDLT(A) decomposition failed");
}

void AnalyticalGAMLibrary::update_LDLT_CAC(const Eigen::SparseMatrix<double>& CAC) {
  if (!ldlt_CAC_analyzed_) {
    ldlt_CAC_.analyzePattern(CAC); // sparsity pattern will not change if X and D are fixed
    ldlt_CAC_analyzed_ = true;
  }
  ldlt_CAC_.factorize(CAC);
  if (ldlt_CAC_.info() != Eigen::Success) Rcpp::stop("LDLT(CAC) decomposition failed");
}

void AnalyticalGAMLibrary::optimize(const Eigen::VectorXd& y, const Eigen::VectorXd& Sm1,
                                        unsigned max_iter, double tol_val) {
  //build fixed matrices
  const auto XtSm1                         = (Sm1.asDiagonal()*X_).transpose();
  const Eigen::SparseMatrix<double> XtSm2X = Eigen::SparseMatrix<double>(K_,K_)
                                             .selfadjointView<Eigen::Lower>().rankUpdate(XtSm1);
  const auto XtSm2y                        = XtSm1 * (Sm1.asDiagonal() * y);
  const Eigen::SparseMatrix<double> K2DtD  = Eigen::SparseMatrix<double>(K_,K_)
                                             .selfadjointView<Eigen::Lower>().rankUpdate(K_*D_.transpose());
  //prepare loop
  Eigen::VectorXd beta_old = beta_;
  unsigned iter = 0;
  has_converged_ = false;
  do {
    //compute A
    auto jitter = Eigen::SparseMatrix<double>(K_,K_);
    jitter.setIdentity();
    jitter *= 1e-10;
    const Eigen::SparseMatrix<double> A = XtSm2X + lambda_*lambda_*K2DtD + jitter;
    
    //compute beta_y = Am1 Xt Sm2 y
    update_LDLT_A(A);
    Eigen::VectorXd beta_y = ldlt_A_.solve(XtSm2y);
    
    //remove constraints if present
    if (neq_>0) {
      //compute decomposition of C Am1 Ct
      auto CAm1Ct = Ceq_ * ldlt_A_.solve(Ceq_.transpose());
      update_LDLT_CAC(CAm1Ct);
      
      //compute constraint term
      Eigen::VectorXd beta_C = ldlt_CAC_.solve(Ceq_*beta_y);
      beta_C = ldlt_A_.solve(Ceq_.transpose()*beta_C);
      
      //subtract from main estimate
      beta_y = beta_y - beta_C;
    }
    
    //compute lambda
    lambda_ = std::sqrt((K_ - 2)/((K_*D_*beta_).squaredNorm()+sigmasqinv_));
    
    //compute relative precision on beta and check convergence
    double precision = beta_.maxCoeff()-beta_.minCoeff();
    if ( precision>0 ) precision = (beta_ - beta_old).array().abs().maxCoeff() / precision;
    has_converged_ = precision < tol_val;
    //bookkeeping
    beta_old = beta_;
    iter++;
  } while (iter < max_iter && (!has_converged_));
}


} // namespace binless
