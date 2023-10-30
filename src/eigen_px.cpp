#include "RcppEigen.h"
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::MatrixXd vecR_ridge_general(
    const Eigen::MappedSparseMatrix<double> L,    //Decomposition of variance L^T L = VAR(alpha)
    const Rcpp::NumericVector pg_mean,
    const Eigen::Map<Eigen::MatrixXd> Z,
    const Eigen::Map<Eigen::MatrixXi> M,
    const Rcpp::NumericVector mapping_J, // Where to assign the elements to the larger matrix.
    const Rcpp::NumericVector d,
    const Eigen::VectorXi start_z,
    bool diag_only
){
  //const Rcpp::List& mapping_z, // The base data split by observation i and level.
  //const Rcpp::List& mapping_to_re,    //The mapping to RE of each individual observation i.
  //Get preliminary variables.
  int nrow_Z = Z.rows();
  
  Rcpp::NumericVector dsq = d * d;
  int size_vecR = Rcpp::sum(dsq);
  int J = d.size();
  
  Eigen::MatrixXd ridge_R = Eigen::MatrixXd::Zero(size_vecR, size_vecR);
  
  //Loop over each observation.
  for (int i = 0; i < nrow_Z; ++i){
    //Matrix to store the output
    Eigen::MatrixXd store_i = Eigen::MatrixXd::Zero(size_vecR, size_vecR);
    
    // List positions_i = mapping_to_re[i];
    // List l_z_i = mapping_z[i];

    Eigen::VectorXd z_i = Z.row(i);
    Eigen::VectorXi m_i = M.row(i);
    
    //Loop over all pairs of random effects.
    for (int j = 0; j < J; ++j){
      for (int jprime = 0; jprime <= j; ++jprime){
        if (diag_only & (jprime != j)){
          continue;         
        }
        int d_j = d[j];
        int d_jprime = d[jprime];

        Eigen::MatrixXd outer_alpha = Eigen::MatrixXd::Zero(d_j, d_jprime);
        Eigen::MatrixXd outer_z = z_i.segment(start_z[j], d_j) * z_i.segment(start_z[jprime], d_jprime).transpose();
        
        
        int m_ij = m_i[j];
        int m_ijprime = m_i[jprime];
        
        for (int k = 0; k < d_j; ++k){
          for (int kprime = 0; kprime < d_jprime; ++kprime){
            outer_alpha(k,kprime) = L.col(m_ij + k - 1).cwiseProduct(L.col(m_ijprime + kprime - 1)).sum();
          }
        }
        
        Eigen::MatrixXd kron = Eigen::kroneckerProduct(outer_alpha, outer_z);

        store_i.block(mapping_J[j], mapping_J[jprime], dsq[j], dsq[jprime]) = kron;
        if (j != jprime){
          store_i.block(mapping_J[jprime], mapping_J[j], dsq[jprime], dsq[j]) = kron.transpose();
        }
      }
    }
    ridge_R = ridge_R + store_i * pg_mean[i];
  }
  
  
  return(ridge_R);
}

// The BULK of time is here.
// Eigen::VectorXd z_ij = l_z_i[j];
// Eigen::VectorXd z_ijprime = l_z_i[jprime];
// Eigen::MatrixXd outer_z = z_ij * z_ijprime.transpose();
//Decent amount of time is here too.
// Eigen::VectorXi index_j = positions_i[j];
// Eigen::VectorXi index_jprime = positions_i[jprime];
// 
// for (int k = 0; k < d_j; ++k){
//   for (int kprime = 0; kprime < d_jprime; ++kprime){
//     int ik = index_j[k] - 1;
//     int ikprime = index_jprime[kprime] - 1;
// 
//     outer_alpha(k,kprime) = L.col(ik).cwiseProduct(L.col(ikprime)).sum();
// 
//   }
// }


// [[Rcpp::export]]
Eigen::MatrixXd vecR_design(
    const Eigen::Map<Eigen::VectorXd> alpha_mu,
    const Eigen::Map<Eigen::MatrixXd> Z,
    const Eigen::Map<Eigen::MatrixXi> M,
    const Rcpp::NumericVector mapping_J, // Where to assign the elements to the larger matrix.
    const Rcpp::NumericVector d,
    const Eigen::VectorXi start_z
){
  int nrow_Z = Z.rows();
  
  Rcpp::NumericVector dsq = d * d;
  int size_vecR = Rcpp::sum(dsq);
  int J = d.size();
  
  Eigen::MatrixXd design_R(nrow_Z, size_vecR);
  
  // Loop over each observation
  for (int i = 0; i < nrow_Z; ++i){

    Eigen::VectorXd z_i = Z.row(i);
    Eigen::VectorXi m_i = M.row(i);
    
    for (int j = 0; j < J; ++j){
      int d_j = d[j];
      
      design_R.block(i, mapping_J[j], 1, dsq[j]) = Eigen::kroneckerProduct(alpha_mu.segment(m_i[j] - 1, d_j), z_i.segment(start_z[j], d_j)).transpose();
    }    
  }

  return design_R;
}



// [[Rcpp::export]]
Eigen::VectorXd vecR_fast_ridge(
    const Eigen::MappedSparseMatrix<double> X,
    const Eigen::MappedSparseMatrix<double> omega,
    const Eigen::MappedSparseMatrix<double> prior_precision,
    const Eigen::Map<Eigen::VectorXd> y,
    const Eigen::Map<Eigen::VectorXd> adjust_y
){
  Eigen::SparseMatrix<double> adj_X = X.adjoint();
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > Ch(adj_X * omega * X + prior_precision);
  Eigen::VectorXd mean = Ch.solve(adj_X * y + adjust_y);
  return mean;
}


// [[Rcpp::export]]
Eigen::MatrixXd vecR_ridge_new(
    const Eigen::MappedSparseMatrix<double> L,    //Decomposition of variance L^T L = VAR(alpha)
    const Eigen::ArrayXd pg_mean,
    const Rcpp::NumericVector mapping_J, // Where to assign the elements to the larger matrix.
    const Rcpp::NumericVector d,
    const Rcpp::List store_id,
    const Rcpp::List store_re_id,
    const Rcpp::List store_design,
    bool diag_only
){
  
  Rcpp::NumericVector dsq = d * d;
  int size_vecR = Rcpp::sum(dsq);
  int J = d.size();
  
  Eigen::MatrixXd ridge_R = Eigen::MatrixXd::Zero(size_vecR, size_vecR);
  
  for (int j = 0; j < J; ++j){
    
    Rcpp::List id_j = store_id[j];
    Rcpp::List re_id_j = store_re_id[j];
    Eigen::MatrixXd design_j = store_design[j];
    int d_j = d[j];
    for (int jprime = 0; jprime <= j; ++jprime){
      
      if (diag_only & (jprime != j)){
        continue;
      }
      
      Eigen::MatrixXd design_jprime = store_design[jprime];
      int d_jprime = d[jprime];
      Rcpp::List id_jjprime = id_j[jprime];
      Rcpp::List re_id_jjprime = re_id_j[jprime];
      int length_combo = id_jjprime.size();
      
      Eigen::MatrixXd kron = Eigen::MatrixXd::Zero(dsq[j], dsq[jprime]) ;
      
      for (int combo = 0; combo < length_combo; ++combo){
        Rcpp::IntegerVector col_re = re_id_jjprime[combo];
        Rcpp::IntegerVector u_comb = id_jjprime[combo];
        Eigen::MatrixXd outer_Z = Eigen::MatrixXd::Zero(d_j, d_jprime);
        for (int i = 0; i < u_comb.size(); ++i){
        // for (IntegerVector::iterator i = u_comb.begin();  i != u_comb.end(); ++i){
          int pos_i = u_comb[i] - 1;
          Eigen::VectorXd row_j = design_j.row(pos_i);
          Eigen::VectorXd row_jprime = design_jprime.row(pos_i);

          Eigen::MatrixXd tmp_outer = row_j * row_jprime.transpose();
          outer_Z += tmp_outer * pg_mean[pos_i];
        }
        Eigen::MatrixXd outer_alpha = Eigen::MatrixXd::Zero(d_j, d_jprime);
        
        int m_j = col_re[0];
        int m_jprime = col_re[1];
        
        for (int k = 0; k < d_j; ++k){
          for (int kprime = 0; kprime < d_jprime; ++kprime){
            outer_alpha(k,kprime) = L.col(m_j + k - 1).cwiseProduct(L.col(m_jprime + kprime - 1)).sum();
          }
        }
        kron += Eigen::kroneckerProduct(outer_alpha, outer_Z);
      }
      
      ridge_R.block(mapping_J[j], mapping_J[jprime], dsq[j], dsq[jprime]) = kron;
      if (j != jprime){
        ridge_R.block(mapping_J[jprime], mapping_J[j], dsq[jprime], dsq[j]) = kron.transpose();
      }
      
    } 
  }
  return(ridge_R);
}
