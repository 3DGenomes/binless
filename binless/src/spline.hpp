#ifndef SPLINE_HPP
#define SPLINE_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

//generate dense cubic B-spline base along x with basis functions built at the knot locations
//code adapted from the original paper by
//Paul H. C. Eilers and Brian D. Marx
//Flexible Smoothing with B-splines and Penalties
//Statistical Science 1996, Vol. 11, No. 2, pp. 89-121
Eigen::MatrixXd cardinal_cubic_bspline_design(const Eigen::VectorXd& x, double dx, const Eigen::RowVectorXd& knots);

//generate sparse cubic spline base along x with K evenly spaced basis functions
// provided lowest and highest spline knots. Set to min and max of x if unknown.
Eigen::SparseMatrix<double> generate_spline_base(const Eigen::VectorXd& x, double xmin, double xmax, unsigned K);

#endif