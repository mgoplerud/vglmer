#include "RcppEigen.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;


// [[Rcpp::export]]
IntegerMatrix unique_rows(const IntegerMatrix m) {

  Rcpp::Environment base("package:base");
  Function do_unique = base["unique"];

  return do_unique(m);
}

// [[Rcpp::export]]
Rcpp::List prepare_Z_for_px(
  Rcpp::IntegerMatrix& Mmap //Integer matrix of group IDs
){
  Rcout  << "A";
  int J = Mmap.cols();
  int N = Mmap.rows();
  
  Rcpp::List store_re_id(J);
  Rcpp::List store_id(J);
  
  IntegerVector out(N);
  
  for (int j = 0; j < J; j++){
    Rcout  << "m"; 
    Rcpp::List store_re_id_j(j);
    Rcpp::List store_id_j(j);
    
    for (int jprime = 0; jprime < j; jprime++){
      
      Rcout << "|";
      
      IntegerMatrix sub_M(N, 2);
      sub_M.column( 0 ) = Mmap.column(j);
      sub_M.column( 1 ) = Mmap.column(jprime);
      IntegerMatrix umap = unique_rows(sub_M);

      int n_unique = umap.rows();

      IntegerVector col_0 = umap.column(0);
      IntegerVector col_1 = umap.column(1);

      IntegerVector position_index(N);

      for (int i = 0; i < N; i++){ //Loop over each observation
        IntegerVector M_i = Mmap.row(i);
        for(int g=0; g < n_unique; g++) {
          IntegerVector umap_g = umap.row(g);
          if ( (umap_g(0) == M_i(0)) && (umap_g(1) == M_i(1)) ){
            position_index(i) = g;
            break;
          }
        }
      }
      
      store_re_id_j[jprime] = umap;
      store_id_j[jprime] = position_index;
    }
    
    store_re_id[j] = store_re_id_j;
    store_id[j] = store_id_j;
  }
  
  return store_id;
}
