#include "RcppEigen.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

//' Cyclical Calculation of Variance Decomposition
List calculate_alpha_decomp_full_factor(
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::MappedSparseMatrix<double> Z,
    const Eigen::Map<Eigen::MatrixXd> P,
    const Eigen::Map<Eigen::VectorXd> omega,
    const Eigen::Map<Eigen::ArrayXd> d_j,
    const Eigen::Map<Eigen::ArrayXd> g_j,
    const List Tinv,
    const Rcpp::List& re_position_list
){
  
  int n_REs = re_position_list.length();
  
  Eigen::VectorXd sqrt_omega_k = omega.cwiseSqrt();
  Eigen::SparseMatrix<double> scaled_Z = sqrt_omega_k.asDiagonal() * Z;
  Eigen::MatrixXd scaled_X = sqrt_omega_k.asDiagonal() * X;

  Eigen::ArrayXd size_lower_tri = d_j * d_j;
  size_lower_tri += d_j;
  size_lower_tri /= 2;
  
  int total_size = (g_j * size_lower_tri).sum();

  Eigen::VectorXd decomp_var(total_size);
  Eigen::VectorXi decomp_i(total_size);
  Eigen::VectorXi decomp_j(total_size);
  
  int prior_j_offset = 0;
  double running_det = 0;

  for (int j = 0; j < n_REs; ++j){
    
    List re_positions_j = re_position_list[j];
    
    Eigen::MatrixXd Tinv_j = Tinv[j];
    
    int size_RE_j = d_j(j);
    int size_out = size_lower_tri(j);
    int n_groups_j = g_j(j);
    
    // For each group g in random effect j.

    for(int g = 0; g < n_groups_j; ++g) {
      
      NumericVector g_prime = re_positions_j[g];
      
      // Compute the (lower triangular view) precision, i.e. M_{jg}^T O M_{jg}
      Eigen::MatrixXd prec_jg = Eigen::MatrixXd::Zero(size_RE_j,size_RE_j);
      
      for (int i = 0; i < size_RE_j; i++){
        for (int i_prime = 0; i_prime <= i; i_prime++){
          int index_i = g_prime[i] - 1;
          int index_i_prime = g_prime[i_prime] - 1;
          
          // (Z - XP)_i = Z_i - X * P.col(i)
          
          double sum_i = (scaled_Z.col(index_i) - scaled_X * P.col(index_i)).cwiseProduct(scaled_Z.col(index_i_prime) - scaled_X * P.col(index_i_prime)).sum();
          
          prec_jg(i, i_prime) = sum_i;
          
        }
      }
      
      if (size_RE_j == 1){
        double LLT_jg = (prec_jg + Tinv_j).sum();
        decomp_var(g + prior_j_offset) = std::sqrt(LLT_jg);
        decomp_i(g + prior_j_offset) = g_prime[0];
        decomp_j(g + prior_j_offset) = g_prime[0];
        running_det += -1.0 * std::log(LLT_jg);
      }else{
        Eigen::MatrixXd LLT_jg = (prec_jg + Tinv_j).llt().matrixL();
        int counter = 0;
        int offset = g * size_out + prior_j_offset;
        for (int c1 = 0; c1 < size_RE_j; c1++){
          for (int c2 = c1; c2 < size_RE_j; c2++){
            decomp_var(offset + counter) = LLT_jg(c2, c1);
            decomp_i(offset + counter) = g_prime[c2];
            decomp_j(offset + counter) = g_prime[c1];
            if (c1 == c2){
              running_det += -2 * std::log(LLT_jg(c2, c1));
            }
            counter += 1;
          }
        }
      }
    }
    prior_j_offset += n_groups_j * size_out;
  }
  
  return List::create(
    Rcpp::Named("log_det") = running_det,
    Rcpp::Named("decomp_var") = decomp_var,
    Rcpp::Named("decomp_i") = decomp_i,
    Rcpp::Named("decomp_j") = decomp_j
  );
}