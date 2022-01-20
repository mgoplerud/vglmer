
#include "RcppEigen.h"
// [[Rcpp::depends(RcppEigen)]]

// #include "Rcpp/Benchmark/Timer.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cpp_inv_alpha_var(
    const Eigen::MappedSparseMatrix<double> diag_vi_pg_mean,
    const Eigen::MatrixXd P,
    const Eigen::MappedSparseMatrix<double> X,
    const Rcpp::List Tinv,
    const Rcpp::List vi_Z_list,
    const double beta_var_lndet
){
  
  int number_of_RE = Tinv.size();
  Eigen::VectorXd log_det_alpha(number_of_RE);
  
  Rcpp::List out_alpha_var(number_of_RE);

  // Rcpp::Timer time_list;
  // time_list.step("start");
  
  for (int j = 0; j < number_of_RE; j++){
    
    // time_list.step("Prepare A and Meat");
    Eigen::SparseMatrix<double> Z_j = vi_Z_list[j];
    Eigen::SparseMatrix<double> Tinv_j = Tinv[j];
    
    Eigen::SparseMatrix<double> A_j = Z_j.adjoint() * diag_vi_pg_mean * Z_j + Tinv_j;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > Ch_A(A_j);
    Eigen::MatrixXd bread_j = Z_j.adjoint() * diag_vi_pg_mean * X;
    Eigen::MatrixXd Ainv_bread = Ch_A.solve(bread_j);
    Eigen::MatrixXd t_bread = bread_j.adjoint();
    Eigen::MatrixXd meat_woodbury = X.adjoint() * diag_vi_pg_mean * X - t_bread * Ainv_bread;
    
    Eigen::LDLT<Eigen::MatrixXd > chol_wood(meat_woodbury);
      
    // Ainv + Ainv %*% bread_j %*% solve(meat_woodbury) %*% t(bread_j) %*% Ainv
    // Ainv %*% ( A + bread_j %*% solve(meat_woodbury, t(bread_j)) ) Ainv
    
    // time_list.step("Solve");
    Eigen::MatrixXd temp_meat = A_j + bread_j * chol_wood.solve(t_bread);
    // time_list.step("Solve 1");
    Eigen::MatrixXd temp_meat_2 = Ch_A.solve(temp_meat).adjoint();
    // time_list.step("Solve 2");
    Eigen::MatrixXd vi_alpha_j = Ch_A.solve(temp_meat_2);
    
    // time_list.step("VectorD");
    
    Eigen::ArrayXd D_ch_A_lndet = Ch_A.vectorD();
    Eigen::ArrayXd D_ch_wood = chol_wood.vectorD();
    
    // time_list.step("Misc");
    double A_lndet = D_ch_A_lndet.log().sum();
    double meat_woodbury_lndet = D_ch_wood.log().sum();
    
    double log_det_alpha_var_j = meat_woodbury_lndet + beta_var_lndet + A_lndet;
    log_det_alpha_var_j *= -1.0;
    
    // time_list.step("Assign");
    out_alpha_var[j] = vi_alpha_j;
    log_det_alpha(j) = log_det_alpha_var_j;
  }
  
  // time_list.step("Stop");
  double out_det_alpha = log_det_alpha.sum();
  
  return Rcpp::List::create(
    Rcpp::Named("variance") = out_alpha_var,
    Rcpp::Named("logdet") = out_det_alpha
  );
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_quad_collapsed(
  const Eigen::MappedSparseMatrix<double> V, // Var(alpha)
  const Rcpp::List& re_position_list,
  const Rcpp::List& Z_list_raw,
  const Rcpp::List& individual_assignments,
  const Eigen::MatrixXd vi_beta_var,
  const Eigen::MatrixXd P,
  const Eigen::MatrixXd X
){
  
  int N = X.rows();
  int size_P = X.cols();
  int n_REs = re_position_list.size();
  
  Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
  Eigen::MatrixXd weights_cross = Eigen::MatrixXd::Zero(size_P, N);

  Eigen::MatrixXd P_varA_tP_plus_beta = P * V * P.adjoint() + vi_beta_var;
  for (int i = 0; i < N; ++i){
    Eigen::VectorXd X_i = X.row(i);    
    out(i) += X_i.adjoint()  * (P_varA_tP_plus_beta * X_i);
  }
  
  for (int j = 0; j < n_REs; ++j){

    Eigen::MatrixXd Z_j = Z_list_raw[j];
    Eigen::ArrayXi level_j = individual_assignments[j];
    level_j -= 1;
    
    // Used to reconstruct var(alpha_{j,g})
    List re_positions_j = re_position_list[j];
    int g_j = re_positions_j.size();
    
    NumericVector first_position = re_positions_j[0];
    int size_RE_j = first_position.size();
    
    
    // For each group g in random effect j.
    for(int g = 0; g < g_j; ++g) {
      
      // For each random effect alpha_{j,g}, get var(alpha_{j,g})
      
      NumericVector g_prime = re_positions_j[g];
      Eigen::MatrixXd var_alpha_j_g = Eigen::MatrixXd::Zero(size_RE_j,size_RE_j);
      Eigen::MatrixXd cov_jg_P(size_P, size_RE_j);
      
      for (int i = 0; i < size_RE_j; i++){
        int index_i = g_prime[i] - 1;
        cov_jg_P.col(i) = P * V.col(index_i);
        for (int i_prime = 0; i_prime <= i; i_prime++){
          int index_i_prime = g_prime[i_prime] - 1;
          double sum_i = V.coeff(index_i, index_i_prime);
          var_alpha_j_g(i, i_prime) = sum_i;
        }
      }
      
      var_alpha_j_g.triangularView<Eigen::StrictlyUpper>() = var_alpha_j_g.adjoint();

      // Add the outer product if i is in group g
      for (int i = 0; i < N; i++){
        if (level_j(i) == g){
          Eigen::VectorXd Z_ij = Z_j.row(i);
          out(i) += Z_ij.adjoint() * (var_alpha_j_g * Z_ij);
          weights_cross.col(i) += cov_jg_P * Z_ij;
        }
      }
      
    }

  }
  
  Eigen::VectorXd vector_weights = weights_cross.adjoint().array().cwiseProduct(X.array()).rowwise().sum();
  out += -2.0 * vector_weights;
  
  return out;
}

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
  Eigen::MatrixXd P_varA_tP_plus_beta = tP.adjoint() * varA * tP + vi_beta_var;
  
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
