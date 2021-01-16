
convert_brm <- function(formula = NULL, data = NULL, brms_args = list(),
                        compile = TRUE, prior_method = 'default_iw'){
  message('Compiling original BRM model')
  
  glue <- glue::glue
  brm <- brms::brm
  
  brms_args$formula <- formula
  brms_args$data <- data
  brms_args$chains <- 0
  brm_object <- do.call('brm', brms_args)

  data_brm <- brms::standata(brm_object)
  
  #Get the code and split by whitespace
  code_brm <- brms::stancode(brm_object)
  code_brm <- strsplit(code_brm, split='\\n')[[1]]

  message('Rewriting Code')
  re_id <- names(data_brm)[grepl(names(data_brm), pattern='^M_[0-9]+$')]
  for (v in re_id){
    id <- gsub(v, pattern='^M_', replacement = '')
    size_v <- data_brm[[v]]
    if (size_v == 1){
      if (prior_method == 'default_iw'){
        data_brm[[paste0('prior_S_', id)]] <- (1 + 1)/2
        data_brm[[paste0('nu_', id)]] <- 1/2
      }
      insert_data <- c(glue("  // ----------------- [Manual Prior for IW {id}]"),
                       glue("  real nu_{id};"),
                       glue("  real prior_S_{id};"),
                       "  // -----------------")
      pos_data <- grep(code_brm, pattern = glue("vector\\[N\\] Z_{id}_1;"))
      if (length(pos_data) > 1){stop('---')}
      code_brm <- c(code_brm[1:(pos_data-1)], insert_data, code_brm[(pos_data):length(code_brm)])  
      
      pos_parameter <- grep(code_brm, pattern=glue("vector<lower=0>\\[M_{id}\\] sd_{id};  // group-level standard deviations"))
      insert_param <- c(glue("  vector<lower=0>[M_{id}] var_{id}; // group-level variance"))
      code_brm <- c(code_brm[1:(pos_parameter - 1)], insert_param,
                    code_brm[(pos_parameter + 1):length(code_brm)])      
      
      pos_model <- grep(code_brm, pattern=glue("vector\\[N_{id}\\] r_{id}_1;  // actual group-level effects"))
      
      code_brm <- c(code_brm[1:pos_model],
                    glue("  vector[M_{id}] sd_{id} = sqrt(var_{id});"),
                    code_brm[(pos_model+1):length(code_brm)])
      
      pos_model <- grep(code_brm, pattern=glue("target \\+= student_t_lpdf\\(sd_{id} \\| 3, 0, 2.5\\)|- 1 \\* student_t_lccdf\\(0 \\| 3, 0, 2.5\\);"))
      
      insert_model <- c(glue("  target += inv_gamma_lpdf(var_{id}[1] | nu_{id}, prior_S_{id});"))
      code_brm <- c(code_brm[1:(pos_model[1]-1)], 
                    insert_model, 
                    code_brm[(pos_model[2] + 1):length(code_brm)])
      

    }else{
      if (prior_method == 'default_iw'){
        data_brm[[paste0('prior_S_', id)]] <- diag(size_v)
        data_brm[[paste0('nu_', id)]] <- size_v + 1
      }
      
      insert_data <- c(glue("  // ----------------- [Manual Prior for IW {id}]"),
                       glue("  cov_matrix[M_1] prior_S_{id};"),
                       glue("  real nu_{id};"),
                       "  // -----------------")
      pos_data <- grep(code_brm, pattern = glue("int<lower=1> NC_{id};"))
      if (length(pos_data) > 1){stop('---')}
      code_brm <- c(code_brm[1:(pos_data-1)], insert_data, code_brm[(pos_data + 1):length(code_brm)])  
      
      pos_parameter <- grep(code_brm, pattern=glue("cholesky_factor_corr\\[M_{id}\\] L_{id};|vector<lower=0>\\[M_{id}\\] sd_{id};"))
      if (diff(pos_parameter) != 2){stop('Misaligned')}
      
      insert_param <- c(glue("  cov_matrix[M_{id}] S_{id};  // covariance matrix"))
      code_brm <- c(code_brm[1:pos_parameter[1]-1], 
                    code_brm[pos_parameter[1] + 1],
                    insert_param,
                    code_brm[pos_parameter[2]:length(code_brm)]
      )
      
      
      pos_transform <- grep(code_brm, pattern=glue("r_{id} = \\(diag_pre_multiply\\(sd_{id}, L_{id}\\) \\* z_{id}\\)';"))
      insert_transform <- c(glue("  cholesky_factor_cov[M_{id}] chol_{id} = cholesky_decompose(S_{id});"),
                            glue("  r_{id} = (chol_{id} * z_{id})';"))
      code_brm <- c(
        code_brm[1:(pos_transform-1)],
        insert_transform,
        code_brm[(pos_transform+1):length(code_brm)]
      )
      pat_model <- c(glue("target \\+= student_t_lpdf\\(sd_{id} \\| 3, 0, 2.5\\)"),
                     glue("- 2 \\* student_t_lccdf\\(0 \\| 3, 0, 2.5\\);"),
                     glue("target \\+= lkj_corr_cholesky_lpdf\\(L_{id} | 1);"))
      pos_model <- grep(code_brm, pattern=paste(pat_model, collapse = '|'))
      insert_model <- c(glue("  target += inv_wishart_lpdf(S_{id} | nu_{id}, prior_S_{id});"))
      code_brm <- c(code_brm[1:(pos_model[1]-1)], code_brm[pos_model[1] + 2],
                    insert_model, code_brm[(pos_model[3]+1):length(code_brm)])
      
      pos_gen <- grep(code_brm, pattern=glue("cor_{id}\\[choose\\(k - 1, 2\\) \\+ j\\] = Cor_{id}\\[j, k\\];"))
      if (length(pos_gen) == 0){stop('posgen missing (cor)')}
      code_brm <- code_brm[-(pos_gen + -3:2)]
      pos_gen <- grep(code_brm, pattern = glue("corr_matrix\\[M_{id}\\] Cor_{id} = multiply_lower_tri_self_transpose\\(L_{id}\\);"))
      if (length(pos_gen) == 0){stop('posgen missing (corr)')}
      code_brm <- code_brm[-c(pos_gen + -1:1)]
      
    }
  }  
  

  pos_chol <- grep(code_brm, pattern=glue("cholesky_factor_cov\\[M_[0-9]+\\] chol_[0-9]+ = cholesky_decompose\\(S_[0-9]+\\);"))
  if (length(pos_chol) > 0){
    select_chol <- code_brm[pos_chol]
    code_brm <- code_brm[-pos_chol]
    code_brm <- c(code_brm[1:(pos_chol[1]-1)], select_chol,
      code_brm[(pos_chol[1] ):length(code_brm)])
  }
  
  
  
  code_brm <- paste(code_brm, collapse='\n')
  
  
  data_brm <- data_brm[!(names(data_brm) %in% paste0('NC', gsub(re_id, pattern='^M', replacement = '')))]
  
  if (compile){
    message('Compiling Again')
    compile_adjusted <- rstan::stan_model(model_code = code_brm, save_dso = TRUE)
    #Check convergence
    convg <- suppressMessages(rstan::sampling(compile_adjusted, data = data_brm, chains = 0))
  }else{
    compile_adjusted <- NULL
  }
  
  output <- list(model = compile_adjusted, data = data_brm,
       code = code_brm, orig_brm = brm_object)
  return(output)
}
