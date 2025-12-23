#include "RcppEigen.h"
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::VectorXd cpp_quad_legacy(
    const Eigen::SparseMatrix<double> tZ,
    const Eigen::SparseMatrix<double> varA,
    const Eigen::MatrixXd tP,
    const Eigen::MatrixXd X,
    const Eigen::MatrixXd vi_beta_var
){
  
  int N = tZ.cols();
  Eigen::VectorXd out(N);
  Eigen::MatrixXd P_varA_tP_plus_beta = tP.adjoint() * varA * tP;
  P_varA_tP_plus_beta += vi_beta_var;
  
  // Loop over each column (observation) in tZ //
  for (int k = 0; k < N; ++k){
    // Loop over each non-zero element
    Eigen::VectorXd X_k = X.row(k);
    Eigen::VectorXd tP_x_k = tP * X_k;
    double quad_k = X_k.adjoint() * P_varA_tP_plus_beta * X_k;
    
    for (Eigen::SparseMatrix<double>::InnerIterator it(tZ, k); it; ++it){
      int row_tZ_pos = it.row();
      double row_tZ_value = it.value();
      
      // Faster version of varA.col(row_tZ_pos).dot(tP_x_k);
      // double inter_value = varA.col(row_tZ_pos).dot(tP_x_k);
      double inter_value = 0;
      for (Eigen::SparseMatrix<double>::InnerIterator it3(varA, row_tZ_pos); it3; ++it3){
        inter_value += tP_x_k(it3.row()) * it3.value();
      }
      
      quad_k += -2.0 * inter_value * row_tZ_value;
      double zt_A_z = 0;
      for (Eigen::SparseMatrix<double>::InnerIterator it2(tZ, k); it2; ++it2){
        zt_A_z += varA.coeff(row_tZ_pos,it2.row()) * row_tZ_value * it2.value();
      }
      quad_k += zt_A_z;
    }
    out(k) = quad_k;
  }
  
  return out;
  
}


// [[Rcpp::export]]
Rcpp::List test_f(
    const Eigen::SparseMatrix<double> diag_vi_pg_mean,
    const Eigen::SparseMatrix<double> design_C,
    const Eigen::SparseMatrix<double> Tinv_C,
    const Eigen::VectorXd s,
    const Rcpp::List vi_M_list
){
  
  Eigen::SparseMatrix<double> t_design_C = design_C.adjoint();  
  Eigen::SparseMatrix<double> inter_matrix_C = t_design_C * diag_vi_pg_mean * design_C + Tinv_C;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > chol_C(
      inter_matrix_C
  );
  
  double log_det_C_var = - chol_C.vectorD().array().log().sum();
  
  Eigen::SparseMatrix<double> ident_C(design_C.cols(), design_C.cols());
  ident_C.setIdentity();
  Eigen::SparseMatrix<double> vi_C_var = chol_C.solve(ident_C);
  Eigen::VectorXd C_hat = chol_C.solve(t_design_C * s);
  
  int J = vi_M_list.length();
  Rcpp::List vi_P(J);
  int size_C = design_C.cols();
  
  for (int j = 0; j < J; j++){
    Eigen::SparseMatrix<double> data_M_j = vi_M_list[j];
    if (data_M_j.cols() == 0){
      vi_P[j] = Eigen::SparseMatrix<double>(size_C, 0);
    }else{
      Eigen::SparseMatrix<double> inter_matrix = t_design_C * diag_vi_pg_mean * data_M_j;
      vi_P[j] = chol_C.solve(inter_matrix);
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("vi_P") = vi_P,
    Rcpp::Named("C_hat") = C_hat,
    Rcpp::Named("vi_C_var") = vi_C_var,
    Rcpp::Named("log_det_C_var") = log_det_C_var
  );
}

