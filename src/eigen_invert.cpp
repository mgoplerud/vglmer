#include "RcppEigen.h"
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::MatrixXd invert_L(
  const Eigen::MatrixXd L,
  const int size_L,
  const bool direct
){
  
  if (direct){
    Eigen::MatrixXd L_inv = L.inverse().triangularView<Eigen::Lower>();  
    return L_inv;
  }else{
    Eigen::MatrixXd L_inv = Eigen::MatrixXd::Identity(size_L, size_L);
    for (int i = 0; i < size_L; i++) {
      L_inv(i, i) = 1.0 / L(i, i);
      for (int j = i + 1; j < size_L; j++) {
        double sum = 0.0;
        for (int k = i; k < j; k++) {
          sum += L(j, k) * L_inv(k, i);
        }
        L_inv(j, i) = -sum / L(j, j);
      }
    }
    return L_inv;
  }
}

// [[Rcpp::export]]
Rcpp::List invert_rowwise(
    const Eigen::MatrixXd X,
    const Eigen::MatrixXd vec_prior,
    const Eigen::MatrixXd RHS,
    const int dim_X,
    const bool return_chol = false
){
  
  // Initialize Storage Matrices
  
  int row_X = X.rows();
  int dim_X_sq = dim_X * dim_X;
  
  Eigen::MatrixXd solve_X(row_X, dim_X);
  Eigen::MatrixXd inv_X(row_X, X.cols());
  Eigen::VectorXd det_inv_X(row_X);
  Eigen::MatrixXd IMatrix(dim_X, dim_X);
  IMatrix.setIdentity();
  
  if (return_chol){
    
    Eigen::MatrixXd inv_L_out(row_X, X.cols());
    
    for (int i = 0; i < row_X; i++){
      Eigen::MatrixXd orig_i = X.row(i);
      orig_i += vec_prior;
      Eigen::Map<Eigen::MatrixXd> X_i(orig_i.data(), dim_X, dim_X);
      Eigen::LLT<Eigen::MatrixXd> llt_of_X(X_i);
      // Get Cholesky factor
      Eigen::MatrixXd llt_X_i = llt_of_X.matrixL();
      // Get the log determinant
      Eigen::ArrayXd L = llt_X_i.diagonal();
      det_inv_X(i) = - L.log().sum();
      // Get the inverse of L
      Eigen::MatrixXd inv_L = invert_L(llt_X_i, dim_X, false);
      // Get the inverse of X = L L^T
      Eigen::MatrixXd inv_Xi = inv_L.transpose() * inv_L;
      // Save out the flattened L^{-1} and X^{-1}
      Eigen::Map<Eigen::ArrayXd> inv_L_flat(inv_L.data(), dim_X_sq, 1);
      inv_L_out.row(i) = inv_L_flat;
      Eigen::Map<Eigen::ArrayXd> inv_Xi_flat(inv_Xi.data(), dim_X_sq, 1);
      inv_X.row(i) = inv_Xi_flat;
      Eigen::VectorXd RHS_i = RHS.row(i);
      solve_X.row(i) = inv_Xi * RHS_i;
      // solve_X.row(i) = llt_of_X.solve(RHS_i);
    }
    
    return Rcpp::List::create(
      Rcpp::Named("inverse") = inv_X,
      Rcpp::Named("Linv") = inv_L_out,
      Rcpp::Named("det") = det_inv_X,
      Rcpp::Named("mean") = solve_X
    );
    
  }else{
    for (int i = 0; i < row_X; i++){
      // Loop over each row (group)
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
      
      Eigen::VectorXd RHS_i = RHS.row(i);
      solve_X.row(i) = llt_of_X.solve(RHS_i);
    }
    
    return Rcpp::List::create(
      Rcpp::Named("inverse") = inv_X,
      Rcpp::Named("det") = det_inv_X,
      Rcpp::Named("mean") = solve_X
    );
    
  }
}


// [[Rcpp::export]]
Rcpp::List decomp_to_var_rowwise(
    const Eigen::MatrixXd X,
    const int dim_X,
    const bool get_lndet = false
){
  
  int row_X = X.rows();
  int dim_X_sq = dim_X * dim_X;
  double lndet = 0.0;
    
  Eigen::MatrixXd var_rowise(row_X, dim_X_sq);
  // Loop over each row
  for (int i = 0; i < row_X; i++){
    Eigen::VectorXd vec_i = X.row(i);
    Eigen::MatrixXd mat_i = Eigen::Map<Eigen::MatrixXd>(vec_i.data(), dim_X, dim_X);
    Eigen::MatrixXd var_i = mat_i.transpose() * mat_i;
    if (get_lndet){
      Eigen::PartialPivLU<Eigen::MatrixXd> lu(var_i);
      lndet  += lu.matrixLU().diagonal().array().abs().log().sum();
    }
    Eigen::Map<Eigen::ArrayXd> flat_i(var_i.data(), dim_X_sq, 1);
    var_rowise.row(i) = flat_i;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("var") = var_rowise,
    Rcpp::Named("lndet") = lndet
  );
}