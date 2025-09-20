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

// [[Rcpp::export]]
Eigen::VectorXd cpp_dense_zVz(
    const Eigen::Map<Eigen::MatrixXd> &X,
    const Eigen::Map<Eigen::MatrixXd> &V) {
  
    int n = X.rows();
    int k = V.rows();
    Eigen::MatrixXd obj = X * V.transpose();
    Eigen::VectorXd out = obj.rowwise().squaredNorm();
    return out;
}

// [[Rcpp::export]]
Eigen::VectorXd cpp_zipped_sum(
    const Eigen::VectorXi pos_U,
    const Eigen::VectorXi pos_V,
    const Eigen::Map<Eigen::MatrixXd> Ut,
    const Eigen::Map<Eigen::MatrixXd> Vt
){
  
  int N = pos_U.size();
  Eigen::VectorXd out(N);
  for (int i = 0; i < N; i++){
    out[i] = Ut.col(pos_U(i)).dot(Vt.col(pos_V(i)));
  }
  return out;
}


// [[Rcpp::export]]
List LinRegChol_fe(
    const Eigen::Map<Eigen::MatrixXd> &X,
    const Eigen::MappedSparseMatrix<double> omega,
    const Eigen::Map<Eigen::VectorXd> y,
    const bool save_chol = true
){
  // The LDLt decompsotion gives t(P) L D t(L) P
  Eigen::LDLT<Eigen::MatrixXd> llt_reg(X.adjoint() * omega * X);
  Eigen::VectorXd mean = llt_reg.solve(X.adjoint() * y);
  if (save_chol == false){
    return List::create(
      Rcpp::Named("mean") = mean
    );
  }
  // Adjust to get LL^T without the permutation
  Eigen::VectorXd sqrtD = llt_reg.vectorD().cwiseSqrt();
  Eigen::MatrixXd lower_l = Eigen::MatrixXd(llt_reg.matrixL());
  int k = lower_l.cols();
  for (int j = 0; j < k; j++) {
    lower_l.col(j) *= sqrtD(j);
  }
  // Get the permutation order
  Eigen::VectorXi perm_order = Eigen::VectorXi::LinSpaced(k, 0, k-1);
  perm_order = llt_reg.transpositionsP() * perm_order;
  
  
  return List::create(
    Rcpp::Named("mean") = mean,
    Rcpp::Named("diag_L") = sqrtD,
    Rcpp::Named("Pindex") = perm_order,
    Rcpp::Named("origL") = lower_l
  );
}
