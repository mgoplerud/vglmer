#include "RcppEigen.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;


// Q = (X^T X + precision)
// "Ch" finds L L^T = Q.
// "Ch.solve" finds "x" such that (X^T X + precision) "x" = X^T Y
// That is, it finds the ridge solution.
// If z is from i.i.d. Gaussian, then
// L^T \alpha = z
// L^T \alpha = z where \alpha \sim N(0, Q^{-1})

//' Linear Regression by Cholesky
//' 
//' Do linear regression of form (X^T O X + P)^{-1} X^T y where O is omega, P is
//' precision. 
//' 
//' @keywords internal
//' 
//' @param X Design Matrix
//' @param omega Polya-Gamma weights
//' @param prior_precision Prior Precision for Regression
//' @param y Outcome
//' @param save_chol Save cholesky factor
// [[Rcpp::export]]
List LinRegChol(
     const Eigen::MappedSparseMatrix<double> X,
     const Eigen::MappedSparseMatrix<double> omega,
     const Eigen::MappedSparseMatrix<double> prior_precision,
     const Eigen::Map<Eigen::VectorXd> y,
     const bool save_chol = true
  ){
  Eigen::SparseMatrix<double> adj_X = X.adjoint();
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > Ch(adj_X * omega * X + prior_precision);
  Eigen::VectorXd mean = Ch.solve(adj_X * y);
  if (save_chol == false){
    return List::create(
      Rcpp::Named("mean") = mean
    );
  }
  //Extract L to get the transformation of the std.normal
  //Into the correct form of N(0, Q^{-1}) and add the mean.
  Eigen::SparseMatrix<double> lower_l = Ch.matrixL();
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> invP = Ch.permutationPinv();
  Eigen::SparseMatrix<double> invperm_L = invP * lower_l;
  Eigen::VectorXd diag_L = lower_l.diagonal();
  
  // double i_sum = 0;
  // int size_L = lower_l.rows();
  // 
  // for (int i = 0; i < (size_L - 1); ++i){
  //   Eigen::VectorXd unit_i = Eigen::VectorXd::Unit(size_L, i);
  //   i_sum += Ch.solve(unit_i).sum();
  // }
  
  // Eigen::SparseMatrix<double> I(lower_l.rows(),lower_l.rows());
  // I.setIdentity();
  // Eigen::SparseMatrix<double> test_inv = Ch.solve(I);
  // 
  
  // Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> P = Ch.permutationP();
  // Eigen::SparseMatrix<double> permuted_L = P * lower_l;
  // Eigen::VectorXd gibbs = Ch.solve(invperm_L * stdnorm) + mean_eb;
  return List::create(///Rcpp::Named("gibbs") = gibbs,
    Rcpp::Named("mean") = mean,
    Rcpp::Named("diag_L") = diag_L,
    Rcpp::Named("invPindex") =  invP.indices(),
    Rcpp::Named("Pindex") = Ch.permutationP().indices(),
    Rcpp::Named("origL") = lower_l,
    Rcpp::Named("L") = invperm_L
  );
}

// [[Rcpp::export]]
List calculate_expected_outer_alpha(
    const Eigen::MappedSparseMatrix<double> L, // L^T L = Var(alpha)
    const Eigen::Map<Eigen::VectorXd> alpha_mu, // E[alpha]
    const Rcpp::List& re_position_list
){
  //Get the number of REs.
  int n_REs = re_position_list.size();
  
  List outer_alpha(n_REs);
  List var_alpha(n_REs);
  List mu_alpha(n_REs);
  
  Rcpp::NumericVector dim_RE = NumericVector(n_REs);
  //For each random effect:
  for (int j = 0; j < n_REs; ++j){
    List re_positions_j = re_position_list[j];
    int g_j = re_positions_j.size();
    
    NumericVector first_position = re_positions_j[0];
    int size_RE_j = first_position.size();
    
    Eigen::MatrixXd summed_var_alpha = Eigen::MatrixXd::Zero(size_RE_j,size_RE_j);
    Eigen::MatrixXd summed_outer_alpha = Eigen::MatrixXd::Zero(size_RE_j,size_RE_j);
    
    List var_alpha_j(g_j);
    
    // For each group g in random effect j.
    for(int g = 0; g < g_j; ++g) {
      //For each random effect alpha_{j,g}.
      //Adjust for zero indexing.
      NumericVector g_prime = re_positions_j[g];

      Eigen::MatrixXd var_alpha_j_g = Eigen::MatrixXd::Zero(size_RE_j,size_RE_j);
      for (int i = 0; i < size_RE_j; i++){
        for (int i_prime = 0; i_prime <= i; i_prime++){
          int index_i = g_prime[i] - 1;
          int index_i_prime = g_prime[i_prime] - 1;
          
          double sum_i = L.col(index_i).cwiseProduct(L.col(index_i_prime)).sum();
          
          var_alpha_j_g(i, i_prime) = sum_i;
            
          summed_var_alpha(i,i_prime) = summed_var_alpha(i,i_prime) + sum_i;
          summed_outer_alpha(i, i_prime) = summed_outer_alpha(i, i_prime) +
            alpha_mu(index_i) * alpha_mu(index_i_prime);
        }
      }

      var_alpha_j_g.triangularView<Eigen::StrictlyUpper>() = var_alpha_j_g.adjoint();
      var_alpha_j[g] = var_alpha_j_g;
    }
    Eigen::MatrixXd oa_j = summed_var_alpha + summed_outer_alpha;
    oa_j.triangularView<Eigen::StrictlyUpper>() = oa_j.adjoint();
    outer_alpha[j] = oa_j;
    var_alpha[j] = var_alpha_j;
    mu_alpha[j] = summed_outer_alpha;
  }
  
  return List::create(
    Rcpp::Named("outer_alpha") = outer_alpha,
    Rcpp::Named("variance_jg") = var_alpha,
    Rcpp::Named("mu_j") = mu_alpha
  );
}

