#include "RcppEigen.h"
#include "spline.hpp"

Eigen::MatrixXd cardinal_cubic_bspline_design(const Eigen::VectorXd& x, double dx, const Eigen::RowVectorXd& knots) {
  const unsigned degree = 3;
  const unsigned N = x.rows();
  const unsigned K = knots.cols();
  Eigen::VectorXi r = Eigen::VectorXi::LinSpaced(K,1,K);
  r(K-1) = 0;
  const Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> R(r);
  const Eigen::ArrayXXd T = knots.replicate(N,1).array();
  const Eigen::ArrayXXd X = x.replicate(1,K).array();
  const Eigen::ArrayXXd P = (X-T)/dx;
  Eigen::ArrayXXd B = ( (T <= X) && (X < T + dx) ).cast<double>();
  for (unsigned k=1; k<=degree; ++k) B = (P * B + (k+1-P) * (B.matrix() * R).array()) / k;
  return B.matrix().leftCols(K-degree-1);
}

Eigen::SparseMatrix<double> generate_spline_base(const Eigen::VectorXd& x, double xmin, double xmax, unsigned K){
  const unsigned degree = 3;
  const double dx = 1.01*(xmax-xmin)/(K-degree);
  const double lower = xmin-dx*0.01;
  const Eigen::RowVectorXd knots = Eigen::RowVectorXd::LinSpaced(K+degree+1,lower-degree*dx,lower+K*dx);
  auto ret = cardinal_cubic_bspline_design(x,dx,knots);
  return ret.sparseView();
}
