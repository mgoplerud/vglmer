#include "RcppEigen.h"
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;


// [[Rcpp::export]]
Eigen::MatrixXd block_diag_product(
    const Eigen::MatrixXd A,
    const Eigen::MatrixXd B,
    const int block_size,
    const int blocks
){
  
  Eigen::MatrixXd out(blocks, block_size * block_size);
  // Loop over each row
  for (int i = 0; i < blocks; i++){
    for (int s = 0; s < block_size; s++){
      Eigen::ArrayXd A_i = A.row(i * block_size + s);
      for (int sprime = 0; sprime < block_size; sprime ++){
        Eigen::ArrayXd B_j = B.row(i * block_size + sprime);
        double out_ij = (A_i * B_j).sum();
        out(i, s + sprime * block_size) = out_ij;
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List invert_rowwise(
    const Eigen::MatrixXd X,
    const Eigen::MatrixXd vec_prior,
    const int dim_X
){
  
  int row_X = X.rows();
  int dim_X_sq = dim_X * dim_X;
  
  Eigen::MatrixXd inv_X(row_X, X.cols());
  Eigen::VectorXd det_inv_X(row_X);
  Eigen::MatrixXd IMatrix(dim_X, dim_X);
  IMatrix.setIdentity();
  // Loop over each row / group
  for (int i = 0; i < row_X; i++){
    Eigen::MatrixXd orig_i = X.row(i);
    orig_i += vec_prior;
    Eigen::Map<Eigen::MatrixXd> X_i(orig_i.data(), dim_X, dim_X);
    Eigen::LLT<Eigen::MatrixXd> llt_of_X(X_i);
    // Get the inverse matrix
    Eigen::MatrixXd inv_Xi = llt_of_X.solve(IMatrix);
    // Get the log-determinant
    Eigen::MatrixXd llt_X_i = llt_of_X.matrixL();
    Eigen::ArrayXd L = llt_X_i.diagonal();
    det_inv_X(i) = - L.log().sum();
    
    Eigen::Map<Eigen::ArrayXd> inv_Xi_flat(inv_Xi.data(), dim_X_sq, 1);
    inv_X.row(i) = inv_Xi_flat;
  }
  return Rcpp::List::create(
    Rcpp::Named("inverse") = inv_X,
    Rcpp::Named("det") = det_inv_X
  );  
}