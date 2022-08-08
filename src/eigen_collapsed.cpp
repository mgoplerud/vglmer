
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

// [[Rcpp::export]]
Eigen::VectorXd cpp_var_lp(
  const Eigen::SparseMatrix<double> design_C,
  const Eigen::SparseMatrix<double> vi_C_uncond,
  const Rcpp::List vi_M_var,
  const Rcpp::List vi_M_list,
  const Rcpp::List vi_P,
  const bool sparse_input,
  const Rcpp::LogicalVector skip_vector
){

  int N = design_C.rows();
  int size_C = design_C.cols();
  Eigen::VectorXd joint_quad(N);
  int J = vi_M_var.length();
  
  // Calculate contribution from collapsed set
  // Note that ( (X * V) * Z).sum(axis = 0) gives x_i^T V z_i
  joint_quad = ( (design_C * vi_C_uncond).cwiseProduct(design_C) ) * Eigen::VectorXd::Ones(size_C);
  
  // Loop over each random effect
  for (int j = 0; j < J; j++){
    
    Eigen::SparseMatrix<double> data_Mj = vi_M_list[j];
    Eigen::SparseMatrix<double> vi_P_j = vi_P[j];
    
    bool skip_j = skip_vector[j];
    if (skip_j){
      continue;
    }
    
    Eigen::VectorXd jq1(N);
    Eigen::VectorXd jq2(N);
    
    if (sparse_input){
      Eigen::SparseMatrix<double> var_Mj = vi_M_var[j];
      jq1 = ( (data_Mj * var_Mj).cwiseProduct( data_Mj ) ) * Eigen::VectorXd::Ones(data_Mj.cols());
      jq2 = ( (data_Mj * (vi_P_j * var_Mj).adjoint()).cwiseProduct( design_C )) * Eigen::VectorXd::Ones(size_C);
    }else{
      Eigen::MatrixXd var_Mj = vi_M_var[j];
      jq1 = ( (data_Mj * var_Mj).cwiseProduct( data_Mj ) ) * Eigen::VectorXd::Ones(data_Mj.cols());
      jq2 = ( (data_Mj * (vi_P_j * var_Mj).adjoint()).cwiseProduct( design_C )) * Eigen::VectorXd::Ones(size_C);
    }
    joint_quad += jq1 - 2 * jq2;
  }
  return joint_quad;
}


// [[Rcpp::export]]
Rcpp::List cpp_update_m_var(
  const Eigen::SparseMatrix<double> diag_vi_pg_mean,
  const Eigen::SparseMatrix<double> design_C,
  const Eigen::SparseMatrix<double> Tinv_C,
  const Rcpp::List list_Tinv_M,
  const Rcpp::List vi_M_list,
  const bool any_collapsed_C,
  const double lndet_C
){
  
  int J = vi_M_list.length();
  Eigen::VectorXd running_log_det_M_var(J);
  Rcpp::List vi_M_var(J);
  Rcpp::List misc(J);
  
  Eigen::MatrixXd meat_C = design_C.adjoint() * diag_vi_pg_mean * design_C + Tinv_C;
    
  for (int j = 0; j < J; j++){
    
    Eigen::SparseMatrix<double> data_M_j = vi_M_list[j];
    
    if (data_M_j.cols() == 0){
      vi_M_var[j] = Eigen::MatrixXd(0,0);
      running_log_det_M_var[j] = 0; 
      continue;
    }
    Eigen::SparseMatrix<double> Tinv_Mj = list_Tinv_M[j];
    
    // # Z_j^T Omega Z_j + R - Z_j^T Omega X (X^T Omega X)^{-1} X^T Omega Z_j
    // # A = Z_j^T Omega Z_j + R
    // # U = - Z_j^T Omega X
    // # C = X^T Omega X^{-1}
    // # V = X^T Omega Z_j
    // # A^{-1} - A^{-1} U (C^{-1} + V A^{-1} U) V A^{-1}
    
    Eigen::SparseMatrix<double> prec_M_j = data_M_j.adjoint() * diag_vi_pg_mean * data_M_j + Tinv_Mj;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > chol_prec(prec_M_j);
    Eigen::SparseMatrix<double> ident_j(prec_M_j.cols(), prec_M_j.cols());
    ident_j.setIdentity();
    
    // Eigen::ArrayXd vec_d_chol = ;
    double lndet_Ainv = - chol_prec.vectorD().array().log().sum();
    
    if (any_collapsed_C){
      Eigen::MatrixXd bread_j = data_M_j.adjoint() * diag_vi_pg_mean * design_C;
      Eigen::MatrixXd meat_woodbury = meat_C - bread_j.adjoint() * chol_prec.solve(bread_j);
      Eigen::MatrixXd scaled_bread = chol_prec.solve(bread_j);
      Eigen::LDLT<Eigen::MatrixXd > chol_wood(meat_woodbury);
      Eigen::MatrixXd var_M_j = scaled_bread * chol_wood.solve(scaled_bread.adjoint());  
        
      // Eigen::ArrayXd vec_d_wood = ;
      double lndet_woodbury = chol_wood.vectorD().array().log().sum();
      var_M_j += chol_prec.solve(ident_j);
      vi_M_var[j] = var_M_j;
      running_log_det_M_var[j] = - (lndet_woodbury + lndet_C - lndet_Ainv);

    }else{
      Eigen::SparseMatrix<double> var_M_j = chol_prec.solve(ident_j);
      vi_M_var[j] = var_M_j;
      running_log_det_M_var[j] = lndet_Ainv;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("vi_M_var") = vi_M_var,
    Rcpp::Named("running_log_det_M_var") = running_log_det_M_var
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_update_c_var(
  const Eigen::SparseMatrix<double> diag_vi_pg_mean,
  const Eigen::SparseMatrix<double> design_C,
  const Eigen::SparseMatrix<double> Tinv_C,
  const Eigen::VectorXd s,
  const Rcpp::List vi_M_list
){

  Eigen::SparseMatrix<double> t_design_C = design_C.adjoint();  
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > chol_C(
      t_design_C * diag_vi_pg_mean * design_C + Tinv_C
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
      vi_P[j] = chol_C.solve(t_design_C * diag_vi_pg_mean * data_M_j);
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("vi_P") = vi_P,
    Rcpp::Named("C_hat") = C_hat,
    Rcpp::Named("vi_C_var") = vi_C_var,
    Rcpp::Named("log_det_C_var") = log_det_C_var
  );
}
  