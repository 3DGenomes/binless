#include "RcppEigen.h"
#include "gam.hpp"

namespace binless {

Eigen::SparseMatrix<double> decreasing_constraint(unsigned K, unsigned Kfree) {
  if (K <= 1 || Kfree >= K) return Eigen::SparseMatrix<double>();
  //build first order difference matrix
  Eigen::MatrixXi D = - Eigen::MatrixXi::Identity(K-1,K);
  D.diagonal(1).setConstant(1);
  //return only last rows
  return - D.bottomRows(K-Kfree).sparseView().cast<double>();
}

Eigen::SparseMatrix<double> second_order_difference_matrix(unsigned K) {
  if (K<=2) return Eigen::SparseMatrix<double>();
  Eigen::MatrixXi D = Eigen::MatrixXi::Identity(K-2,K);
  D.diagonal(1).setConstant(-2);
  D.diagonal(2).setConstant(1);
  return D.sparseView().cast<double>();
}


} // namespace binless
