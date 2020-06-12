#include "RcppEigen.h"
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::VectorXd cpp_zVz(
  const Eigen::MappedSparseMatrix<double> Z,
  const Eigen::MappedSparseMatrix<double> V
){

  Eigen::SparseMatrix<double> VZ_t = V * Z.adjoint();
  int N = Z.rows();
  Eigen::VectorXd output(N);
  
  for (int i = 0; i < N; i++){
    output(i) = VZ_t.col(i).squaredNorm();
  }
  return output;
  
  // for (int j = 0; j < N; j++){
  //   double norm_j = 0;
  //   for (Eigen::SparseMatrix<double>::InnerIterator i_(VZ_t, j); i_; ++i_){
  //     norm_j += std::pow(i_.value(), 2)
  //   }
  //   output(j) = norm_j;
  // }
}