
if (FALSE){

vi_alpha_outer <- loop_outer_alpha(vi_alpha_chol,outer_alpha_RE_positions)
vi_alpha_2 <- slow_outer_alpha(vi_alpha_var, outer_alpha_RE_positions)

mapply(vi_alpha_2, vi_alpha_outer, FUN=function(a,b){all.equal(as(a, 'dgCMatrix'),as(b, 'dgCMatrix'))}) 
}
