#include "RcppEigen.h"
#include "eigen_helpers.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

Eigen::SparseMatrix<double> sparse_diag(
    Eigen::ArrayXd m
){
  
  int n = m.size();
  
  Eigen::SparseMatrix<double> res(n, n);
  
  std::vector<Eigen::Triplet<double> > tripletList;
  
  tripletList.reserve(n);
  
  for (int i = 0; i < n; i++) {
    
    tripletList.push_back(Eigen::Triplet<double>(i, i, m[i]));
    
  }
  
  res.setFromTriplets(tripletList.begin(), tripletList.end());
  
  return res;
}
