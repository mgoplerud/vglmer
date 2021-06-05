#include "RcppEigen.h"
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List chol_sparse(
  const Eigen::MappedSparseMatrix<double> X,
  const Eigen::MappedSparseMatrix<double> omega,
  const Eigen::MappedSparseMatrix<double> precision
){
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > Ch(X.adjoint() * omega * X + precision);

  Eigen::SparseMatrix<double> lower_l = Ch.matrixL();
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> invP = Ch.permutationPinv();
  Eigen::VectorXd diag_L = lower_l.diagonal();
  
  return List::create(
    Rcpp::Named("diag_L") = diag_L,
    Rcpp::Named("Pindex") = Ch.permutationP().indices(),
    Rcpp::Named("origL") = lower_l
  );
}

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