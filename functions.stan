functions {
  
  /* Censoring */
  
  tuple(vector, matrix, vector, matrix, matrix)
  prep_multi_cond_post(vector yo, vector yc,
                       array[] real to, array[] real tc,
                       array[] real tpred,
                       array[] int to_is, array[] int tc_is,
                       array[] int Jo, array[] int Jc,
                       int Jpred,
                       real magnitude_mu, real length_scale_mu,
                       real magnitude_eta, real length_scale_eta,
                       real sigma)
  {
    int n = num_elements(Jc);
    int Npred = n * Jpred;
    int Nc = sum(Jc);
    int No = sum(Jo);

    matrix[Jpred, Jpred] K_eta_pred =
      gp_exp_quad_cov(tpred, magnitude_eta, length_scale_eta);
    matrix[Npred, Npred] cov_eta =
      block_mat_AB(K_eta_pred, - K_eta_pred / (n-1), n);

    matrix[Npred, Nc] cov_eta_c;
    int row_start, row_end, col_start, col_end;
    row_end = 0;
    for (i in 1:n) {
      row_start = row_end + 1;
      row_end = row_end + Jpred;
      col_end = 0;
      for (j in 1:n) {
        if (Jc[j] == 0)
          continue;
        col_start = col_end + 1;
        col_end = col_end + Jc[j];
        if (i == j) {
          cov_eta_c[row_start:row_end, col_start:col_end] =
            gp_exp_quad_cov(tpred, tc[col_start:col_end],
                            magnitude_eta, length_scale_eta);
        } else {
          cov_eta_c[row_start:row_end, col_start:col_end] =
            - gp_exp_quad_cov(tpred, tc[col_start:col_end],
                              magnitude_eta, length_scale_eta) / (n-1);
        }
      }
    }

    matrix[Npred, No] cov_eta_o;
    row_end = 0;
    for (i in 1:n) {
      row_start = row_end + 1;
      row_end = row_end + Jpred;
      col_end = 0;
      for (j in 1:n) {
        if (Jo[j] == 0)
          continue;
        col_start = col_end + 1;
        col_end = col_end + Jo[j];
        if (i == j) {
          cov_eta_o[row_start:row_end, col_start:col_end] =
            gp_exp_quad_cov(tpred, to[col_start:col_end],
                            magnitude_eta, length_scale_eta);
        } else {
          cov_eta_o[row_start:row_end, col_start:col_end] =
            - gp_exp_quad_cov(tpred, to[col_start:col_end],
                              magnitude_eta, length_scale_eta) / (n-1);
        }
      }
    }

    matrix[Nc, Nc] cov_c;
    row_end = 0;
    for (i in 1:n) {
      if (Jc[i] == 0)
        continue;
      row_start = row_end + 1;
      row_end = row_end + Jc[i];
      col_end = 0;
      for (j in 1:n) {
        if (Jc[j] == 0)
          continue;
        col_start = col_end + 1;
        col_end = col_end + Jc[j];
        if (i == j) {
          cov_c[row_start:row_end, col_start:col_end] =
            add_diag(gp_exp_quad_cov(tc[row_start:row_end], tc[col_start:col_end],
                                     magnitude_mu, length_scale_mu) +
                     gp_exp_quad_cov(tc[row_start:row_end], tc[col_start:col_end],
                                     magnitude_eta, length_scale_eta),
                     sigma^2);
        } else {
          cov_c[row_start:row_end, col_start:col_end] =
            gp_exp_quad_cov(tc[row_start:row_end], tc[col_start:col_end],
                            magnitude_mu, length_scale_mu) -
            gp_exp_quad_cov(tc[row_start:row_end], tc[col_start:col_end],
                            magnitude_eta, length_scale_eta) / (n-1);
        }
      }
    }

    matrix[Nc, No] cov_c_o;
    row_end = 0;
    for (i in 1:n) {
      if (Jc[i] == 0)
        continue;
      row_start = row_end + 1;
      row_end = row_end + Jc[i];
      col_end = 0;
      array[Jc[i]] real tc_slice = tc[row_start:row_end];
      for (j in 1:n) {
        if (Jo[j] == 0)
          continue;
        col_start = col_end + 1;
        col_end = col_end + Jo[j];
        array[Jo[j]] real to_slice = to[col_start:col_end];
        if (i == j) {
          cov_c_o[row_start:row_end, col_start:col_end] =
            gp_exp_quad_cov(tc_slice, to_slice,
                            magnitude_mu, length_scale_mu) +
            gp_exp_quad_cov(tc_slice, to_slice,
                            magnitude_eta, length_scale_eta);
        } else {
          cov_c_o[row_start:row_end, col_start:col_end] =
            gp_exp_quad_cov(tc_slice, to_slice,
                            magnitude_mu, length_scale_mu) -
            gp_exp_quad_cov(tc_slice, to_slice,
                            magnitude_eta, length_scale_eta) / (n-1);
        }
      }
    }

    matrix[No, No] cov_o;
    row_end = 0;
    for (i in 1:n) {
      if (Jo[i] == 0)
        continue;
      row_start = row_end + 1;
      row_end = row_end + Jo[i];
      col_end = 0;
      for (j in 1:n) {
        if (Jo[j] == 0)
          continue;
        col_start = col_end + 1;
        col_end = col_end + Jo[j];
        if (i == j) {
          cov_o[row_start:row_end, col_start:col_end] =
            add_diag(gp_exp_quad_cov(to[row_start:row_end], to[col_start:col_end],
                                     magnitude_mu, length_scale_mu) +
                     gp_exp_quad_cov(to[row_start:row_end], to[col_start:col_end],
                                     magnitude_eta, length_scale_eta),
                     sigma^2);
        } else {
          cov_o[row_start:row_end, col_start:col_end] =
            gp_exp_quad_cov(to[row_start:row_end], to[col_start:col_end],
                            magnitude_mu, length_scale_mu) -
            gp_exp_quad_cov(to[row_start:row_end], to[col_start:col_end],
                            magnitude_eta, length_scale_eta) / (n-1);
        }
      }
    }

    matrix[Npred, Npred] cov_eta_cond_o =
      cov_eta - cov_eta_o * (cov_o \ cov_eta_o');

    matrix[Nc, Nc] cov_c_cond_o =
      cov_c - cov_c_o * (cov_o \ cov_c_o');

    matrix[Npred, Nc] cov_eta_c_cond_o =
      cov_eta_c - cov_eta_o * (cov_o \ cov_c_o');

    int Npred_except_1_group = Npred - Jpred;

    vector[Npred] mean_eta_cond_o = cov_eta_o * (cov_o \ yo);
    vector[Nc] mean_c_cond_o = cov_c_o * (cov_o \ yo);

    matrix[Jpred, Jpred] K_mu_pred =
      gp_exp_quad_cov(tpred, magnitude_mu, length_scale_mu);
    matrix[Jpred, No] K_mu_pred_o =
      gp_exp_quad_cov(tpred, to, magnitude_mu, length_scale_mu);
    matrix[Jpred, Nc] K_mu_pred_c =
      gp_exp_quad_cov(tpred, tc, magnitude_mu, length_scale_mu);

    matrix[Jpred, Jpred] cov_mu_cond_o =
      K_mu_pred - K_mu_pred_o * (cov_o \ K_mu_pred_o');
    matrix[Jpred, Nc] cov_mu_c_cond_o =
      K_mu_pred_c - K_mu_pred_o * (cov_o \ cov_c_o');

    vector[Jpred] mean_mu_cond_o = K_mu_pred_o * (cov_o \ yo);

    vector[Npred] mean_mueta_cond_o =
      append_row(mean_mu_cond_o, mean_eta_cond_o[:Npred_except_1_group]);

    matrix[Npred, Npred] cov_mueta_cond_o;
    cov_mueta_cond_o[:Jpred, :Jpred] = cov_mu_cond_o;
    cov_mueta_cond_o[(Jpred+1):, (Jpred+1):] =
      cov_eta_cond_o[:Npred_except_1_group, :Npred_except_1_group];
    cov_mueta_cond_o[:Jpred, (Jpred+1):] =
      - K_mu_pred_o * (cov_o \ cov_eta_o[:Npred_except_1_group, ]');
    cov_mueta_cond_o[(Jpred+1):, :Jpred] = cov_mueta_cond_o[:Jpred, (Jpred+1):]';

    matrix[Npred, Nc] cov_mueta_c_cond_o;
    cov_mueta_c_cond_o[:Jpred, ] = K_mu_pred_c - K_mu_pred_o * (cov_o \ cov_c_o');
    cov_mueta_c_cond_o[(Jpred+1):, ] = cov_eta_c[:Npred_except_1_group, ] -
      cov_eta_o[:Npred_except_1_group, ] * (cov_o \ cov_c_o');

    /* Term for multiplying on the truncated normal variable P */
    matrix[Npred, Nc] p_factor_mueta = cov_mueta_c_cond_o / cov_c_cond_o;

    matrix[Npred, Npred] cov_q_mueta = cov_mueta_cond_o -
      cov_mueta_c_cond_o * (cov_c_cond_o \ cov_mueta_c_cond_o');

    return (mean_c_cond_o,
            cov_c_cond_o,
            mean_mueta_cond_o,
            p_factor_mueta,
            cov_q_mueta);
  }

  tuple(real, vector, matrix)
  prep_multi_cens_log_lik(vector yo,
                          array[] real to, array[] real tc,
                          array[] int to_is, array[] int tc_is,
                          array[] int Jo, array[] int Jc,
                          real magnitude_mu, real length_scale_mu,
                          real magnitude_eta, real length_scale_eta,
                          real sigma)
  {
    int n = num_elements(Jc); /* number of groups */
    int No = sum(Jo);
    int Nc = sum(Jc);

    int row_start, row_end, col_start, col_end;

    matrix[Nc, Nc] cov_c;
    row_end = 0;
    for (i in 1:n) {
      if (Jc[i] == 0)
        continue;
      row_start = row_end + 1;
      row_end = row_end + Jc[i];
      col_end = 0;
      for (j in 1:n) {
        if (Jc[j] == 0)
          continue;
        col_start = col_end + 1;
        col_end = col_end + Jc[j];
        if (i == j) {
          cov_c[row_start:row_end, col_start:col_end] =
            add_diag(gp_exp_quad_cov(tc[row_start:row_end], tc[col_start:col_end],
                                     magnitude_mu, length_scale_mu) +
                     gp_exp_quad_cov(tc[row_start:row_end], tc[col_start:col_end],
                                     magnitude_eta, length_scale_eta),
                     sigma^2);
        } else {
          cov_c[row_start:row_end, col_start:col_end] =
            gp_exp_quad_cov(tc[row_start:row_end], tc[col_start:col_end],
                            magnitude_mu, length_scale_mu) -
            gp_exp_quad_cov(tc[row_start:row_end], tc[col_start:col_end],
                            magnitude_eta, length_scale_eta) / (n-1);
        }
      }
    }

    matrix[Nc, No] cov_c_o;
    row_end = 0;
    for (i in 1:n) {
      if (Jc[i] == 0)
        continue;
      row_start = row_end + 1;
      row_end = row_end + Jc[i];
      col_end = 0;
      array[Jc[i]] real tc_slice = tc[row_start:row_end];
      for (j in 1:n) {
        if (Jo[j] == 0)
          continue;
        col_start = col_end + 1;
        col_end = col_end + Jo[j];
        array[Jo[j]] real to_slice = to[col_start:col_end];
        if (i == j) {
          cov_c_o[row_start:row_end, col_start:col_end] =
            gp_exp_quad_cov(tc_slice, to_slice,
                            magnitude_mu, length_scale_mu) +
            gp_exp_quad_cov(tc_slice, to_slice,
                            magnitude_eta, length_scale_eta);
        } else {
          cov_c_o[row_start:row_end, col_start:col_end] =
            gp_exp_quad_cov(tc_slice, to_slice,
                            magnitude_mu, length_scale_mu) -
            gp_exp_quad_cov(tc_slice, to_slice,
                            magnitude_eta, length_scale_eta) / (n-1);
        }
      }
    }

    matrix[No, No] cov_o;
    row_end = 0;
    for (i in 1:n) {
      if (Jo[i] == 0)
        continue;
      row_start = row_end + 1;
      row_end = row_end + Jo[i];
      col_end = 0;
      for (j in 1:n) {
        if (Jo[j] == 0)
          continue;
        col_start = col_end + 1;
        col_end = col_end + Jo[j];
        if (i == j) {
          cov_o[row_start:row_end, col_start:col_end] =
            add_diag(gp_exp_quad_cov(to[row_start:row_end], to[col_start:col_end],
                                     magnitude_mu, length_scale_mu) +
                     gp_exp_quad_cov(to[row_start:row_end], to[col_start:col_end],
                                     magnitude_eta, length_scale_eta),
                     sigma^2);
        } else {
          cov_o[row_start:row_end, col_start:col_end] =
            gp_exp_quad_cov(to[row_start:row_end], to[col_start:col_end],
                            magnitude_mu, length_scale_mu) -
            gp_exp_quad_cov(to[row_start:row_end], to[col_start:col_end],
                            magnitude_eta, length_scale_eta) / (n-1);
        }
      }
    }

    real lpdf_o = multi_normal_lpdf(yo | rep_vector(0, No), cov_o);
    vector[Nc] mean_c_cond_o = cov_c_o * (cov_o \ yo);
    matrix[Nc, Nc] cov_c_cond_o = cov_c - cov_c_o * (cov_o \ cov_c_o');

    return (lpdf_o,
            mean_c_cond_o,
            cov_c_cond_o);
  }
                          

  /* yo1o2 must contain (yo1, yo2) */
  tuple(real, vector, matrix)
  prep_multi_cens_log_lik_base(vector yo1o2,
                               array[] real to1, array[] real to2, array[] real tc,
                               array[] int to2_is, array[] int tc_is,
                               int Jo1, array[] int Jo2, array[] int Jc,
                               real magnitude_mu, real length_scale_mu,
                               real magnitude_eta, real length_scale_eta,
                               real sigma)
  {
    int n = num_elements(Jc); /* number of groups */
    int No1 = n * Jo1;
    int No2 = sum(Jo2);
    int No = No1 + No2;
    int Nc = sum(Jc);
    matrix[Jo1, Jo1] K_mu_o1 =
      gp_exp_quad_cov(to1, magnitude_mu, length_scale_mu);
    matrix[Jo1, Jo1] K_eta_o1 =
      gp_exp_quad_cov(to1, magnitude_eta, length_scale_eta);
    matrix[Jo1, Jo1] diag_cov_o1 = add_diag(K_mu_o1 + K_eta_o1, sigma^2);
    matrix[Jo1, Jo1] off_diag_cov_o1 = K_mu_o1 - K_eta_o1 / (n-1.0);
    matrix[No, No] cov_o;
    cov_o[:No1, :No1] = block_mat_AB(diag_cov_o1, off_diag_cov_o1, n);
    cov_o[:No1, (No1+1):] =
      block_rep(gp_exp_quad_cov(to1, to2, magnitude_mu, length_scale_mu) -
                gp_exp_quad_cov(to1, to2, magnitude_eta, length_scale_eta) / (n-1.0),
                n, 1);
    int row_start, row_end, col_start, col_end;
    /* Modify diagonal entries */
    row_end = 0;
    col_end = 0;    
    for (i in 1:n) {
      if (Jo2[i] == 0)
        continue;
      row_start = row_end + 1;
      row_end = row_end + Jo1;
      col_start = col_end + 1;
      col_end = col_end + Jo2[i];
      cov_o[row_start:row_end, (col_start+No1):(col_end+No1)] +=
        n/(n-1.0) * gp_exp_quad_cov(to1, to2[col_start:col_end],
                                    magnitude_eta, length_scale_eta);
    }

    cov_o[(No1+1):, :No1] = cov_o[:No1, (No1+1):]';
    cov_o[(No1+1):, (No1+1):] =
      irregular_cov_mat_B(No2, n, n, Jo2, to2, to2_is,
                          magnitude_mu, length_scale_mu,
                          magnitude_eta, length_scale_eta,
                          sigma);

    matrix[No, No] L_o = cholesky_decompose(cov_o);

    matrix[Nc, Nc] cov_c =
      irregular_cov_mat_B(Nc, n, n, Jc, tc, tc_is,
                          magnitude_mu, length_scale_mu,
                          magnitude_eta, length_scale_eta,
                          sigma);

    matrix[No, Nc] cov_oc;

    cov_oc[:No1, ] =
      block_rep(gp_exp_quad_cov(to1, tc, magnitude_mu, length_scale_mu) -
                gp_exp_quad_cov(to1, tc, magnitude_eta, length_scale_eta) / (n-1.0),
                n, 1);
    /* Modify diagonal entries */
    row_end = 0;
    col_end = 0;
    for (i in 1:n) {
      if (Jc[i] == 0)
        continue;
      row_start = row_end + 1;
      row_end = row_end + Jo1;
      col_start = col_end + 1;
      col_end = col_end + Jc[i];
      cov_oc[row_start:row_end, col_start:col_end] +=
        n/(n-1.0) * gp_exp_quad_cov(to1, tc[col_start:col_end],
                                    magnitude_eta, length_scale_eta);
    }

    /* cov(Yo2, Yc) */
    row_end = No1;
    for (i in 1:n) {
      if (Jo2[i] == 0)
        continue;
      col_end = 0;
      row_start = row_end + 1;
      row_end = row_end + Jo2[i];
      array[Jo2[i]] real slice_o2 = t_b_slice(i, to2, to2_is);
      for (j in 1:n) {
        if (Jc[j] == 0)
          continue;
        array[Jc[j]] real slice_c = t_b_slice(j, tc, tc_is);
        col_start = col_end + 1;
        col_end = col_end + Jc[j];
        if (i == j) {
          cov_oc[row_start:row_end, col_start:col_end] =
            gp_exp_quad_cov(slice_o2, slice_c, magnitude_mu, length_scale_mu) +
            gp_exp_quad_cov(slice_o2, slice_c, magnitude_eta, length_scale_eta);
        } else {
          cov_oc[row_start:row_end, col_start:col_end] =
            gp_exp_quad_cov(slice_o2, slice_c, magnitude_mu, length_scale_mu) -
            gp_exp_quad_cov(slice_o2, slice_c, magnitude_eta, length_scale_eta) / (n-1.0);
        }
      }
    }
    
    matrix[Nc, No] covco_covoinv = cholesky_left_divide_mat(L_o, cov_oc)';

    return (multi_normal_cholesky_lpdf(yo1o2 | rep_vector(0, No), L_o),
            covco_covoinv * yo1o2,
            cov_c - covco_covoinv * cov_oc);
  }

  
  

  /* Returns the log-lik of yo and the mean and cov for the CDF of the other
     term */
  tuple(real, vector, matrix)
  prep_cens_log_lik(vector yo, array[] real to, array[] real tc,
                    real magnitude, real length_scale, real sigma)
  {
    /* Note: these n correspond to J in the multi setup */
    int no = num_elements(to);
    int nc = num_elements(tc);
    matrix[no, no] cov_o = add_diag(gp_exp_quad_cov(to, magnitude, length_scale), sigma^2);
    matrix[no, no] L_o = cholesky_decompose(cov_o);
    matrix[nc, nc] cov_c = add_diag(gp_exp_quad_cov(tc, magnitude, length_scale), sigma^2);
    matrix[nc, no] cov_co = gp_exp_quad_cov(tc, to, magnitude, length_scale);
    matrix[nc, no] covco_covoinv = cholesky_left_divide_mat(L_o, cov_co')';

    return (multi_normal_cholesky_lpdf(yo | rep_vector(0, no), L_o),
            covco_covoinv * yo,
            cov_c - covco_covoinv * cov_co');
  }

  tuple(vector, matrix, vector, matrix, matrix)
  prep_cond_post(vector yo, vector yc,
                 array[] real to,
                 array[] real tc,
                 array[] real t_pred,
                 real magnitude,
                 real length_scale,
                 real sigma)
  {
    /* Note: these n correspond to J in the multi setup */
    int no = num_elements(to);
    int nc = num_elements(tc);
    int n_pred = num_elements(t_pred);
    matrix[n_pred, n_pred] K_pred =
      gp_exp_quad_cov(t_pred, magnitude, length_scale);
    matrix[nc, nc] cov_c =
      add_diag(gp_exp_quad_cov(tc, magnitude, length_scale), sigma^2);
    matrix[n_pred, no] K_pred_o =
      gp_exp_quad_cov(t_pred, to, magnitude, length_scale);
    matrix[n_pred, nc] K_pred_c =
      gp_exp_quad_cov(t_pred, tc, magnitude, length_scale);    
    matrix[nc, no] K_c_o =
      gp_exp_quad_cov(tc, to, magnitude, length_scale);
    matrix[no, no] cov_o =
      add_diag(gp_exp_quad_cov(to, magnitude, length_scale), sigma^2);
    matrix[no, no] L_o = cholesky_decompose(cov_o);
    vector[no] covoinv_yo = cholesky_left_divide_vec(L_o, yo);
    vector[n_pred] mu_f_cond_o = K_pred_o * covoinv_yo;
    vector[nc] mu_c_cond_o = K_c_o * covoinv_yo;
    
    /* One could use Cholesky for divisions below and save/reuse results. */
    matrix[n_pred, n_pred] cov_f_cond_o = K_pred - K_pred_o * (cov_o \ K_pred_o');
    matrix[nc, nc] cov_c_cond_o = cov_c - K_c_o * (cov_o \ K_c_o');
    matrix[n_pred, nc] cov_fyc_cond_o =
      K_pred_c - K_pred_o * (cov_o \ K_c_o');
    matrix[n_pred, nc] covfyccondo_covccondoinv = cov_fyc_cond_o / cov_c_cond_o;
    matrix[n_pred, n_pred] cov_q =
      cov_f_cond_o - covfyccondo_covccondoinv * cov_fyc_cond_o';

    return (mu_f_cond_o,
            covfyccondo_covccondoinv,
            mu_c_cond_o,
            cov_c_cond_o,
            cov_q);
  }

  /* Utility functions */

  /* returns last eta to make them sum to 0 */
  vector
  get_last_eta(matrix mu_eta_pred)
  {
    int n_group = cols(mu_eta_pred) - 1;
    return mu_eta_pred[, 2:n_group] * rep_vector(-1, n_group - 1);
  }

  /* Return block matrix with repetitions of same diagonal block and
     off diagonal block */
  matrix
  block_mat_AB(matrix diagonal_block, matrix off_diagonal_block, int n_blocks)
  {
    int block_rows = rows(diagonal_block);
    int block_cols = cols(diagonal_block);
    matrix[n_blocks * block_rows, n_blocks * block_cols] block_mat;
    int col_start, col_end, row_start, row_end;
    for (i_col in 1:n_blocks) {
      col_start = (i_col-1) * block_cols + 1;
      col_end = col_start + block_cols - 1;
      for (i_row in 1:n_blocks) {
        row_start = (i_row-1) * block_rows + 1;
        row_end = row_start + block_rows - 1;
        if (i_row == i_col)
          block_mat[row_start:row_end, col_start:col_end] = diagonal_block;
        else          
          block_mat[row_start:row_end, col_start:col_end] = off_diagonal_block;
      }
    }

    return block_mat;
  }


  /* Return block matrix with n_rows x n_cols blocks of the matrix block */
  matrix
  block_rep(matrix block, int n_rows, int n_cols)
  {
    int block_cols = cols(block);
    int block_rows = rows(block);
    int col_start, col_end, row_start, row_end;
    matrix[n_rows * block_rows, n_cols * block_cols] result;
    for (i_row in 1:n_rows) {
      row_start = (i_row-1) * block_rows + 1;
      row_end = row_start + block_rows - 1;
      for (i_col in 1:n_cols) {
        col_start = (i_col-1) * block_cols + 1;
        col_end = col_start + block_cols - 1;
        result[row_start:row_end, col_start:col_end] = block;
      }
    }

    return result;
  }

  matrix
  irregular_cov_mat_B(int N_b, int n_b, int n,
                      array[] int J_b, array[] real t_b, array[] int t_b_is,
                      real magnitude_mu, real length_scale_mu,
                      real magnitude_eta, real length_scale_eta,
                      real sigma)
  {
    matrix[N_b, N_b] B;
    int col_start = 0;
    int col_end = 0;
    int row_start = 0;
    int row_end = 0;
    for (i in 1:n_b) {
      col_start = col_end + 1;
      col_end = col_start + J_b[i] - 1;
      row_start = 0;
      row_end = 0;
      array[J_b[i]] real t_b_slice_i = t_b_slice(i, t_b, t_b_is);
      for (j in 1:n_b) {
        row_start = row_end + 1;
        row_end = row_start + J_b[j] - 1;
        if (i == j) {
          B[col_start:col_end, row_start:row_end] = add_diag(
            gp_exp_quad_cov(t_b_slice_i, magnitude_mu, length_scale_mu) +
              gp_exp_quad_cov(t_b_slice_i, magnitude_eta, length_scale_eta),
            max([square(sigma), 1e-8]));
        } else {
          array[J_b[j]] real t_b_slice_j = t_b_slice(j, t_b, t_b_is);
          B[col_start:col_end, row_start:row_end] =
            gp_exp_quad_cov(t_b_slice_i,
                            t_b_slice_j,
                            magnitude_mu, length_scale_mu) -
            gp_exp_quad_cov(t_b_slice_i,
                            t_b_slice_j,
                            magnitude_eta, length_scale_eta) / (n - 1.0);
        }
      }
    }

    return B;
  }
  
  int
  t_b_slice_lwr(int i, array[] int t_b_is)
  {
    return t_b_is[i]+1;
  }

  int
  t_b_slice_upr(int i, array[] int t_b_is)
  {
    return t_b_is[i+1];
  }

  /* Slice i is of size J_b[i] */
  array[] real
  t_b_slice(int i, array[] real t_b, array[] int t_b_is)
  {
    return t_b[t_b_slice_lwr(i, t_b_is):t_b_slice_upr(i, t_b_is)];
  }

  matrix
  cholesky_left_divide_mat(matrix L, matrix A)
  {
    return mdivide_right_tri_low(
      mdivide_left_tri_low(L, A)',
      L
    )';
  }

  vector
  cholesky_left_divide_vec(matrix L, vector v)
  {
    return mdivide_right_tri_low(
      mdivide_left_tri_low(L, v)',
      L
    )';
  }
}

