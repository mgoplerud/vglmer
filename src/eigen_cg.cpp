#include "RcppEigen.h"
#include "eigen_helpers.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

//' Conjugate Gradient for Large VI Problems
// [[Rcpp::export]]
List cg_custom(
    const Eigen::MappedSparseMatrix<double> X,
    const Eigen::MappedSparseMatrix<double> Z,
    const Eigen::Map<Eigen::MatrixXd> P,
    const Eigen::Map<Eigen::VectorXd> omega,
    const Eigen::MappedSparseMatrix<double> ridge_Z,
    const Eigen::MappedSparseMatrix<double> ridge_X,
    const Eigen::Map<Eigen::VectorXd> s,
    const Eigen::Map<Eigen::VectorXd> offset_ridge_X,
    const Eigen::Map<Eigen::VectorXd> old_alpha,
    const double tol,
    const int it_max = 0,
    const int low_dimension = 5
){
  
  int p_Z = Z.cols();
  
  Eigen::VectorXd new_alpha(p_Z);
  Eigen::VectorXd sqrt_omega_k = omega.cwiseSqrt();
  Eigen::VectorXd adj_s(X.rows());
  
  Eigen::SparseMatrix<double> ridge = ridge_Z + P.adjoint() * ridge_X * P;
    
  adj_s = s.cwiseQuotient(sqrt_omega_k);
  Eigen::SparseMatrix<double> scaled_Z = sparse_diag(sqrt_omega_k) * Z;
  Eigen::SparseMatrix<double> scaled_X = sparse_diag(sqrt_omega_k) * X;
  
  Eigen::VectorXd diag_ridge = ridge.diagonal();
  
  Eigen::VectorXd Z_sqnorm = Eigen::VectorXd::Zero(Z.cols());
  for (int k=0; k < scaled_Z.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(scaled_Z, k); it; ++it){
      int temp_col = it.col();
      double temp_dbl = std::pow(it.value(), 2);
      Z_sqnorm(temp_col) += temp_dbl;
    }
  }

  Eigen::MatrixXd meat_X = (scaled_X.adjoint() * scaled_X) * P;
  Eigen::VectorXd P_sqnorm = meat_X.cwiseProduct(P).colwise().sum();
  
  Eigen::SparseMatrix<double> meat_cross = scaled_X.adjoint() * scaled_Z;
  Eigen::VectorXd cross_norm = Eigen::VectorXd::Zero(p_Z);
  for (int k=0; k < meat_cross.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(meat_cross, k); it; ++it){
      int temp_col = it.col();
      int temp_row = it.row();
      double temp_val = it.value();
      cross_norm(temp_col) +=  temp_val * P(temp_row, temp_col);
    }
  }
  Eigen::VectorXd precond_diag = (Z_sqnorm + P_sqnorm + 2 * cross_norm + diag_ridge).cwiseInverse();

  // Eigen::VectorXd precond_diag_old(p_Z);
  // for (int i = 0; i < p_Z; i++){
  //   Eigen::VectorXd P_i = P.col(i);
  //   Eigen::VectorXd int_vector = scaled_Z.col(i) + scaled_X * P_i;
  //   double int_term = (int_vector).squaredNorm() + diag_ridge(i);
  //   precond_diag_old(i) = 1.0/(int_term);
  // }
  
  Eigen::VectorXd alpha = old_alpha;
  Eigen::VectorXd PA = P * alpha;
  // residual = s - (Z * a - X * P * a)
  Eigen::VectorXd residual = adj_s - scaled_Z * alpha + scaled_X * PA;

  // rhsNorm2 = (Z - X P)^T s = Z^T s - P^T X^T s
  Eigen::VectorXd XTS = scaled_X.adjoint() * s;
  double rhsNorm2 = (scaled_Z.adjoint() * adj_s - P.adjoint() * XTS).squaredNorm();
  // normal_residual = (Z - X P)^T residual = Z^T
  Eigen::VectorXd XTr = scaled_X.adjoint() * residual;
  Eigen::VectorXd normal_residual = scaled_Z.adjoint() * residual - P.adjoint() * XTr  - ridge * alpha + offset_ridge_X;

  Eigen::VectorXd p = normal_residual.cwiseProduct(precond_diag);
  double absNew = normal_residual.dot(p);
    
  double threshold = tol * tol * rhsNorm2;
  double tol_error = 0;
    
  int internal_cg_it;
  if (it_max == 0){
    internal_cg_it = p_Z;
  }else{
    internal_cg_it = it_max;
  }
    
  // Do conjugate gradient with normal equation
  // 
  // (M^T O M + R) beta = M^T y but never form M^T O M or even M
  int it = 0;
  bool convg = false;
    
  while (it < internal_cg_it){
    
    // Xp = nrow(X) by 1
    Eigen::MatrixXd Pp = P * p;
    // tmp = (Z - X P) p
    Eigen::VectorXd tmp = scaled_Z * p - scaled_X * Pp;
    
    double cg_alpha_denom = tmp.squaredNorm() + p.adjoint() * ridge * p;
    double cg_alpha = absNew / cg_alpha_denom;
    
    
    alpha += cg_alpha * p;
    residual -= cg_alpha * tmp;
    
    // XTr = ncol(p) by 1
    Eigen::VectorXd XTr = scaled_X.adjoint() * residual;
    Eigen::VectorXd normal_residual = scaled_Z.adjoint() * residual - P.adjoint() * XTr  - ridge * alpha + offset_ridge_X;

    double residualNorm2 = normal_residual.squaredNorm();
    tol_error = std::sqrt(residualNorm2 / rhsNorm2);
    
    if (residualNorm2 < threshold){
      convg = true;
      break;
    }
    
    Eigen::VectorXd z = normal_residual.cwiseProduct(precond_diag);
    double absOld = absNew;
    absNew = normal_residual.dot(z);
    double step_cg_alpha = absNew / absOld;
    
    p = z + step_cg_alpha * p;
    
    it++;
  }

  return List::create(
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("iter") = it,
    Rcpp::Named("error") = tol_error,
    Rcpp::Named("converged") = convg,
    Rcpp::Named("precond") = precond_diag
  );
 
}