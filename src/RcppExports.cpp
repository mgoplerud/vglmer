// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// LinRegChol
List LinRegChol(const Eigen::MappedSparseMatrix<double> X, const Eigen::MappedSparseMatrix<double> omega, const Eigen::MappedSparseMatrix<double> prior_precision, const Eigen::Map<Eigen::VectorXd> y, const bool save_chol);
RcppExport SEXP _vglmer_LinRegChol(SEXP XSEXP, SEXP omegaSEXP, SEXP prior_precisionSEXP, SEXP ySEXP, SEXP save_cholSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type prior_precision(prior_precisionSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool >::type save_chol(save_cholSEXP);
    rcpp_result_gen = Rcpp::wrap(LinRegChol(X, omega, prior_precision, y, save_chol));
    return rcpp_result_gen;
END_RCPP
}
// decomp_calculate_expected_outer_alpha
List decomp_calculate_expected_outer_alpha(const Eigen::MappedSparseMatrix<double> L, const Eigen::Map<Eigen::VectorXd> alpha_mu, const Rcpp::List& re_position_list, const Eigen::MatrixXd tP, const Eigen::MatrixXd L_beta, const bool do_adjustment);
RcppExport SEXP _vglmer_decomp_calculate_expected_outer_alpha(SEXP LSEXP, SEXP alpha_muSEXP, SEXP re_position_listSEXP, SEXP tPSEXP, SEXP L_betaSEXP, SEXP do_adjustmentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type L(LSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type alpha_mu(alpha_muSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type re_position_list(re_position_listSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type tP(tPSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type L_beta(L_betaSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_adjustment(do_adjustmentSEXP);
    rcpp_result_gen = Rcpp::wrap(decomp_calculate_expected_outer_alpha(L, alpha_mu, re_position_list, tP, L_beta, do_adjustment));
    return rcpp_result_gen;
END_RCPP
}
// direct_calculate_expected_outer_alpha
List direct_calculate_expected_outer_alpha(const Eigen::MappedSparseMatrix<double> V, const Eigen::Map<Eigen::VectorXd> alpha_mu, const Rcpp::List& re_position_list);
RcppExport SEXP _vglmer_direct_calculate_expected_outer_alpha(SEXP VSEXP, SEXP alpha_muSEXP, SEXP re_position_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type V(VSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type alpha_mu(alpha_muSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type re_position_list(re_position_listSEXP);
    rcpp_result_gen = Rcpp::wrap(direct_calculate_expected_outer_alpha(V, alpha_mu, re_position_list));
    return rcpp_result_gen;
END_RCPP
}
// cg_custom
List cg_custom(const Eigen::MappedSparseMatrix<double> X, const Eigen::MappedSparseMatrix<double> Z, const Eigen::Map<Eigen::MatrixXd> P, const Eigen::Map<Eigen::VectorXd> omega, const Eigen::MappedSparseMatrix<double> ridge_Z, const Eigen::MappedSparseMatrix<double> ridge_X, const Eigen::Map<Eigen::VectorXd> s, const Eigen::Map<Eigen::VectorXd> offset_ridge_X, const Eigen::Map<Eigen::VectorXd> old_alpha, const double tol, const int it_max, const int low_dimension);
RcppExport SEXP _vglmer_cg_custom(SEXP XSEXP, SEXP ZSEXP, SEXP PSEXP, SEXP omegaSEXP, SEXP ridge_ZSEXP, SEXP ridge_XSEXP, SEXP sSEXP, SEXP offset_ridge_XSEXP, SEXP old_alphaSEXP, SEXP tolSEXP, SEXP it_maxSEXP, SEXP low_dimensionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type P(PSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type ridge_Z(ridge_ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type ridge_X(ridge_XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type s(sSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type offset_ridge_X(offset_ridge_XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type old_alpha(old_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< const int >::type low_dimension(low_dimensionSEXP);
    rcpp_result_gen = Rcpp::wrap(cg_custom(X, Z, P, omega, ridge_Z, ridge_X, s, offset_ridge_X, old_alpha, tol, it_max, low_dimension));
    return rcpp_result_gen;
END_RCPP
}
// cpp_inv_alpha_var
Rcpp::List cpp_inv_alpha_var(const Eigen::MappedSparseMatrix<double> diag_vi_pg_mean, const Eigen::MatrixXd P, const Eigen::MappedSparseMatrix<double> X, const Rcpp::List Tinv, const Rcpp::List vi_Z_list, const double beta_var_lndet);
RcppExport SEXP _vglmer_cpp_inv_alpha_var(SEXP diag_vi_pg_meanSEXP, SEXP PSEXP, SEXP XSEXP, SEXP TinvSEXP, SEXP vi_Z_listSEXP, SEXP beta_var_lndetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type diag_vi_pg_mean(diag_vi_pg_meanSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type P(PSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type Tinv(TinvSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type vi_Z_list(vi_Z_listSEXP);
    Rcpp::traits::input_parameter< const double >::type beta_var_lndet(beta_var_lndetSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_inv_alpha_var(diag_vi_pg_mean, P, X, Tinv, vi_Z_list, beta_var_lndet));
    return rcpp_result_gen;
END_RCPP
}
// cpp_quad_collapsed
Eigen::MatrixXd cpp_quad_collapsed(const Eigen::MappedSparseMatrix<double> V, const Rcpp::List& re_position_list, const Rcpp::List& Z_list_raw, const Rcpp::List& individual_assignments, const Eigen::MatrixXd vi_beta_var, const Eigen::MatrixXd P, const Eigen::MatrixXd X);
RcppExport SEXP _vglmer_cpp_quad_collapsed(SEXP VSEXP, SEXP re_position_listSEXP, SEXP Z_list_rawSEXP, SEXP individual_assignmentsSEXP, SEXP vi_beta_varSEXP, SEXP PSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type V(VSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type re_position_list(re_position_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Z_list_raw(Z_list_rawSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type individual_assignments(individual_assignmentsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type vi_beta_var(vi_beta_varSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type P(PSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_quad_collapsed(V, re_position_list, Z_list_raw, individual_assignments, vi_beta_var, P, X));
    return rcpp_result_gen;
END_RCPP
}
// cpp_quad_legacy
Eigen::VectorXd cpp_quad_legacy(const Eigen::SparseMatrix<double> tZ, const Eigen::SparseMatrix<double> varA, const Eigen::MatrixXd tP, const Eigen::MatrixXd X, const Eigen::MatrixXd vi_beta_var);
RcppExport SEXP _vglmer_cpp_quad_legacy(SEXP tZSEXP, SEXP varASEXP, SEXP tPSEXP, SEXP XSEXP, SEXP vi_beta_varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double> >::type tZ(tZSEXP);
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double> >::type varA(varASEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type tP(tPSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type vi_beta_var(vi_beta_varSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_quad_legacy(tZ, varA, tP, X, vi_beta_var));
    return rcpp_result_gen;
END_RCPP
}
// cpp_var_lp
Eigen::VectorXd cpp_var_lp(const Eigen::SparseMatrix<double> design_C, const Eigen::SparseMatrix<double> vi_C_uncond, const Rcpp::List vi_M_var, const Rcpp::List vi_M_list, const Rcpp::List vi_P, const bool sparse_input, const Rcpp::LogicalVector skip_vector);
RcppExport SEXP _vglmer_cpp_var_lp(SEXP design_CSEXP, SEXP vi_C_uncondSEXP, SEXP vi_M_varSEXP, SEXP vi_M_listSEXP, SEXP vi_PSEXP, SEXP sparse_inputSEXP, SEXP skip_vectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double> >::type design_C(design_CSEXP);
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double> >::type vi_C_uncond(vi_C_uncondSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type vi_M_var(vi_M_varSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type vi_M_list(vi_M_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type vi_P(vi_PSEXP);
    Rcpp::traits::input_parameter< const bool >::type sparse_input(sparse_inputSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type skip_vector(skip_vectorSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_var_lp(design_C, vi_C_uncond, vi_M_var, vi_M_list, vi_P, sparse_input, skip_vector));
    return rcpp_result_gen;
END_RCPP
}
// cpp_update_m_var
Rcpp::List cpp_update_m_var(const Eigen::SparseMatrix<double> diag_vi_pg_mean, const Eigen::SparseMatrix<double> design_C, const Eigen::SparseMatrix<double> Tinv_C, const Rcpp::List list_Tinv_M, const Rcpp::List vi_M_list, const bool any_collapsed_C, const double lndet_C);
RcppExport SEXP _vglmer_cpp_update_m_var(SEXP diag_vi_pg_meanSEXP, SEXP design_CSEXP, SEXP Tinv_CSEXP, SEXP list_Tinv_MSEXP, SEXP vi_M_listSEXP, SEXP any_collapsed_CSEXP, SEXP lndet_CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double> >::type diag_vi_pg_mean(diag_vi_pg_meanSEXP);
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double> >::type design_C(design_CSEXP);
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double> >::type Tinv_C(Tinv_CSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type list_Tinv_M(list_Tinv_MSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type vi_M_list(vi_M_listSEXP);
    Rcpp::traits::input_parameter< const bool >::type any_collapsed_C(any_collapsed_CSEXP);
    Rcpp::traits::input_parameter< const double >::type lndet_C(lndet_CSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_update_m_var(diag_vi_pg_mean, design_C, Tinv_C, list_Tinv_M, vi_M_list, any_collapsed_C, lndet_C));
    return rcpp_result_gen;
END_RCPP
}
// test_f
Rcpp::List test_f(const Eigen::SparseMatrix<double> diag_vi_pg_mean, const Eigen::SparseMatrix<double> design_C, const Eigen::SparseMatrix<double> Tinv_C, const Eigen::VectorXd s, const Rcpp::List vi_M_list);
RcppExport SEXP _vglmer_test_f(SEXP diag_vi_pg_meanSEXP, SEXP design_CSEXP, SEXP Tinv_CSEXP, SEXP sSEXP, SEXP vi_M_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double> >::type diag_vi_pg_mean(diag_vi_pg_meanSEXP);
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double> >::type design_C(design_CSEXP);
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double> >::type Tinv_C(Tinv_CSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type s(sSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type vi_M_list(vi_M_listSEXP);
    rcpp_result_gen = Rcpp::wrap(test_f(diag_vi_pg_mean, design_C, Tinv_C, s, vi_M_list));
    return rcpp_result_gen;
END_RCPP
}
// calculate_alpha_decomp_full_factor
List calculate_alpha_decomp_full_factor(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::MappedSparseMatrix<double> Z, const Eigen::Map<Eigen::MatrixXd> P, const Eigen::Map<Eigen::VectorXd> omega, const Eigen::Map<Eigen::ArrayXd> d_j, const Eigen::Map<Eigen::ArrayXd> g_j, const List Tinv, const Rcpp::List& re_position_list);
RcppExport SEXP _vglmer_calculate_alpha_decomp_full_factor(SEXP XSEXP, SEXP ZSEXP, SEXP PSEXP, SEXP omegaSEXP, SEXP d_jSEXP, SEXP g_jSEXP, SEXP TinvSEXP, SEXP re_position_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type P(PSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXd> >::type d_j(d_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXd> >::type g_j(g_jSEXP);
    Rcpp::traits::input_parameter< const List >::type Tinv(TinvSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type re_position_list(re_position_listSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_alpha_decomp_full_factor(X, Z, P, omega, d_j, g_j, Tinv, re_position_list));
    return rcpp_result_gen;
END_RCPP
}
// chol_sparse
Rcpp::List chol_sparse(const Eigen::MappedSparseMatrix<double> X, const Eigen::MappedSparseMatrix<double> omega, const Eigen::MappedSparseMatrix<double> precision);
RcppExport SEXP _vglmer_chol_sparse(SEXP XSEXP, SEXP omegaSEXP, SEXP precisionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type precision(precisionSEXP);
    rcpp_result_gen = Rcpp::wrap(chol_sparse(X, omega, precision));
    return rcpp_result_gen;
END_RCPP
}
// cpp_zVz
Eigen::VectorXd cpp_zVz(const Eigen::MappedSparseMatrix<double> Z, const Eigen::MappedSparseMatrix<double> V);
RcppExport SEXP _vglmer_cpp_zVz(SEXP ZSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_zVz(Z, V));
    return rcpp_result_gen;
END_RCPP
}
// cpp_zAz_nonfact
Eigen::VectorXd cpp_zAz_nonfact(const Eigen::MatrixXd Z, const Eigen::MappedSparseMatrix<double> A);
RcppExport SEXP _vglmer_cpp_zAz_nonfact(SEXP ZSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_zAz_nonfact(Z, A));
    return rcpp_result_gen;
END_RCPP
}
// vecR_ridge_general
Eigen::MatrixXd vecR_ridge_general(const Eigen::MappedSparseMatrix<double> L, const Rcpp::NumericVector pg_mean, const Eigen::Map<Eigen::MatrixXd> Z, const Eigen::Map<Eigen::MatrixXi> M, const Rcpp::NumericVector mapping_J, const Rcpp::NumericVector d, const Eigen::VectorXi start_z, bool diag_only);
RcppExport SEXP _vglmer_vecR_ridge_general(SEXP LSEXP, SEXP pg_meanSEXP, SEXP ZSEXP, SEXP MSEXP, SEXP mapping_JSEXP, SEXP dSEXP, SEXP start_zSEXP, SEXP diag_onlySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type L(LSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type pg_mean(pg_meanSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXi> >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type mapping_J(mapping_JSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type start_z(start_zSEXP);
    Rcpp::traits::input_parameter< bool >::type diag_only(diag_onlySEXP);
    rcpp_result_gen = Rcpp::wrap(vecR_ridge_general(L, pg_mean, Z, M, mapping_J, d, start_z, diag_only));
    return rcpp_result_gen;
END_RCPP
}
// vecR_design
Eigen::MatrixXd vecR_design(const Eigen::Map<Eigen::VectorXd> alpha_mu, const Eigen::Map<Eigen::MatrixXd> Z, const Eigen::Map<Eigen::MatrixXi> M, const Rcpp::NumericVector mapping_J, const Rcpp::NumericVector d, const Eigen::VectorXi start_z);
RcppExport SEXP _vglmer_vecR_design(SEXP alpha_muSEXP, SEXP ZSEXP, SEXP MSEXP, SEXP mapping_JSEXP, SEXP dSEXP, SEXP start_zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type alpha_mu(alpha_muSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXi> >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type mapping_J(mapping_JSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type start_z(start_zSEXP);
    rcpp_result_gen = Rcpp::wrap(vecR_design(alpha_mu, Z, M, mapping_J, d, start_z));
    return rcpp_result_gen;
END_RCPP
}
// vecR_fast_ridge
Eigen::VectorXd vecR_fast_ridge(const Eigen::MappedSparseMatrix<double> X, const Eigen::MappedSparseMatrix<double> omega, const Eigen::MappedSparseMatrix<double> prior_precision, const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::VectorXd> adjust_y);
RcppExport SEXP _vglmer_vecR_fast_ridge(SEXP XSEXP, SEXP omegaSEXP, SEXP prior_precisionSEXP, SEXP ySEXP, SEXP adjust_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type prior_precision(prior_precisionSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type adjust_y(adjust_ySEXP);
    rcpp_result_gen = Rcpp::wrap(vecR_fast_ridge(X, omega, prior_precision, y, adjust_y));
    return rcpp_result_gen;
END_RCPP
}
// vecR_ridge_new
Eigen::MatrixXd vecR_ridge_new(const Eigen::MappedSparseMatrix<double> L, const Eigen::ArrayXd pg_mean, const Rcpp::NumericVector mapping_J, const Rcpp::NumericVector d, const Rcpp::List store_id, const Rcpp::List store_re_id, const Rcpp::List store_design, bool diag_only);
RcppExport SEXP _vglmer_vecR_ridge_new(SEXP LSEXP, SEXP pg_meanSEXP, SEXP mapping_JSEXP, SEXP dSEXP, SEXP store_idSEXP, SEXP store_re_idSEXP, SEXP store_designSEXP, SEXP diag_onlySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type L(LSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXd >::type pg_mean(pg_meanSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type mapping_J(mapping_JSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type store_id(store_idSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type store_re_id(store_re_idSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type store_design(store_designSEXP);
    Rcpp::traits::input_parameter< bool >::type diag_only(diag_onlySEXP);
    rcpp_result_gen = Rcpp::wrap(vecR_ridge_new(L, pg_mean, mapping_J, d, store_id, store_re_id, store_design, diag_only));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_vglmer_LinRegChol", (DL_FUNC) &_vglmer_LinRegChol, 5},
    {"_vglmer_decomp_calculate_expected_outer_alpha", (DL_FUNC) &_vglmer_decomp_calculate_expected_outer_alpha, 6},
    {"_vglmer_direct_calculate_expected_outer_alpha", (DL_FUNC) &_vglmer_direct_calculate_expected_outer_alpha, 3},
    {"_vglmer_cg_custom", (DL_FUNC) &_vglmer_cg_custom, 12},
    {"_vglmer_cpp_inv_alpha_var", (DL_FUNC) &_vglmer_cpp_inv_alpha_var, 6},
    {"_vglmer_cpp_quad_collapsed", (DL_FUNC) &_vglmer_cpp_quad_collapsed, 7},
    {"_vglmer_cpp_quad_legacy", (DL_FUNC) &_vglmer_cpp_quad_legacy, 5},
    {"_vglmer_cpp_var_lp", (DL_FUNC) &_vglmer_cpp_var_lp, 7},
    {"_vglmer_cpp_update_m_var", (DL_FUNC) &_vglmer_cpp_update_m_var, 7},
    {"_vglmer_test_f", (DL_FUNC) &_vglmer_test_f, 5},
    {"_vglmer_calculate_alpha_decomp_full_factor", (DL_FUNC) &_vglmer_calculate_alpha_decomp_full_factor, 8},
    {"_vglmer_chol_sparse", (DL_FUNC) &_vglmer_chol_sparse, 3},
    {"_vglmer_cpp_zVz", (DL_FUNC) &_vglmer_cpp_zVz, 2},
    {"_vglmer_cpp_zAz_nonfact", (DL_FUNC) &_vglmer_cpp_zAz_nonfact, 2},
    {"_vglmer_vecR_ridge_general", (DL_FUNC) &_vglmer_vecR_ridge_general, 8},
    {"_vglmer_vecR_design", (DL_FUNC) &_vglmer_vecR_design, 6},
    {"_vglmer_vecR_fast_ridge", (DL_FUNC) &_vglmer_vecR_fast_ridge, 5},
    {"_vglmer_vecR_ridge_new", (DL_FUNC) &_vglmer_vecR_ridge_new, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_vglmer(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
