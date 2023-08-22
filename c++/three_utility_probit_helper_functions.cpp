#include <RcppArmadillo.h>
#include <cmath>
#include <RcppDist.h>
#include <mvtnorm.h>
#include <R_ext/Rdynload.h>
//Code from RcppTN: https://github.com/olmjo/RcppTN/blob/master/src/rtn1.cpp
#include "rtn1.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::depends(RcppArmadillo, RcppDist, mvtnorm)]]

const double pi2 = pow(datum::pi,2);
const double TWOPI = 6.283185307179586;

mat create_ar_1_m(
    double term_length, double rho, 
    double tau) {
  
  mat ar_1_kernel(term_length, term_length);
  for (int i = 0; i < term_length; i++) {
    for (int j = i; j < term_length; j++) {
      ar_1_kernel(i, j) = tau / (1 - pow(rho, 2)) * pow(rho, j - i);
      ar_1_kernel(j, i) = ar_1_kernel(i, j);
    }
  }
  return(ar_1_kernel);
}

mat create_ar_1_m_chol(double term_length, double rho, double tau) {
  mat ar_1_kernel_chol(term_length, term_length);
  if (rho == 0 || term_length == 1) {
    return(ar_1_kernel_chol.eye(term_length, term_length) * sqrt(tau / (1 - pow(rho, 2))));
  }
  for (int i = 0; i < term_length; i++) {
    for (int j = i; j < term_length; j++) {
      ar_1_kernel_chol(i, j) = sqrt(tau) * pow(rho, j - i);
    }
  }
  for (int j = 0; j < term_length; j++) {
    ar_1_kernel_chol(0, j) = ar_1_kernel_chol(0, j) / sqrt(1 - pow(rho, 2));
  }
  return(ar_1_kernel_chol);
}

mat create_ar_1_m_inverse(double term_length, double rho, double tau) {
  mat inv_m(term_length, term_length);
  if (rho == 0 || term_length == 1) {
    return(inv_m.eye(term_length, term_length) * (1 - pow(rho, 2)) / tau);
  }
  inv_m(0,0) = 1;
  inv_m(term_length - 1, term_length - 1) = 1;
  inv_m(0,1) = -rho;
  inv_m(term_length - 1, term_length - 2) = -rho;
  if (term_length == 2) {
    return(inv_m / tau);
  }
  for (int i = 1; i < term_length - 1; i++) {
    inv_m(i, i - 1) = -rho;
    inv_m(i, i) = 1 + pow(rho, 2);
    inv_m(i, i + 1) = -rho;
  }
  return(inv_m / tau);
}

vec simulate_draw_from_ar_1_m_chol(
    double term_length, double rho, double tau, 
    double mean) {
  
  mat chol_m = create_ar_1_m_chol(term_length, rho, tau);
  vec draw(term_length, fill::randn);
  return(chol_m.t() * draw + mean);
}

double sample_three_utility_probit_beta(
    rowvec y_star_m_1, rowvec y_star_m_3, 
    rowvec alpha_v_1, rowvec alpha_v_2,
    rowvec delta_v_1, rowvec delta_v_2,
    double beta_mean, double beta_s) { 
  
  y_star_m_1 = y_star_m_1 - alpha_v_1 % delta_v_1;
  y_star_m_3 = y_star_m_3 - alpha_v_2 % delta_v_2;
  double post_var = 1.0 / pow(beta_s, 2) + 
    dot(alpha_v_1, alpha_v_1) + dot(alpha_v_2, alpha_v_2);
  double post_mean = beta_mean / pow(beta_s, 2) - 
    dot(alpha_v_1, y_star_m_1) - dot(alpha_v_2, y_star_m_3);
  
  return(randn() / sqrt(post_var) + post_mean / post_var);
}

vec sample_three_utility_probit_beta_gp(
    rowvec y_star_m_1, rowvec y_star_m_3, 
    rowvec alpha_v_1, rowvec alpha_v_2,
    rowvec delta_v_1, rowvec delta_v_2, 
    uvec case_year, double rho) {
  
  int years_served = max(case_year) - min(case_year) + 1;
  mat ar_1_m_inv = create_ar_1_m_inverse(years_served, rho, 1 - rho * rho);
  y_star_m_1 = y_star_m_1 - alpha_v_1 % delta_v_1;
  y_star_m_3 = y_star_m_3 - alpha_v_2 % delta_v_2;

  vec post_mean(years_served, fill::zeros);
  for (int i = 0; i < case_year.n_elem; i++) {
    ar_1_m_inv(case_year(i), case_year(i)) += 
      alpha_v_1(i) * alpha_v_1(i) + alpha_v_2(i) * alpha_v_2(i);
    post_mean(case_year(i)) -=
      alpha_v_1(i) * y_star_m_1(i) + alpha_v_2(i) * y_star_m_3(i);
  }
  post_mean = solve(ar_1_m_inv, post_mean);
  return(rmvnorm(1, post_mean, ar_1_m_inv.i()).t());
}

double logit(double p) {
  return(log(p) - log(1 - p));
}

double inv_logit(double z) {
  return(1.0 / (1.0 + exp(-z)));
}

double sample_rho_pos_logit_gibbs(
    double rho, vec ideal_pos_1_m, 
    uvec judge_start_ind, uvec judge_end_ind,
    double rho_mean, 
    double rho_sigma, double rho_sd) {
  
  double next_rho = inv_logit(logit(rho) + rho_sd * randn());
  double next_log_ll = 
    d_truncnorm(next_rho, rho_mean, rho_sigma, 0, 1, 1) +
                  log(next_rho) + log(1 - next_rho);
  double prev_log_ll = 
    d_truncnorm(rho, rho_mean, rho_sigma, 0, 1, 1) +
                  log(rho) + log(1 - rho);
  for (int i = 0; i < judge_start_ind.n_elem; i++) {
    rowvec pos_v = ideal_pos_1_m(span(judge_start_ind(i),
                                      judge_end_ind(i))).t();
    
    mat prev_ar_1_m = create_ar_1_m(pos_v.n_elem, rho, 1 - rho * rho);
    prev_log_ll += as_scalar(dmvnorm(pos_v, zeros(pos_v.n_elem), prev_ar_1_m, true));
    
    mat next_ar_1_m = create_ar_1_m(pos_v.n_elem, next_rho, 1 - next_rho * next_rho);
    next_log_ll += as_scalar(dmvnorm(pos_v, zeros(pos_v.n_elem), next_ar_1_m, true));
  }
  if (log(randu()) < next_log_ll - prev_log_ll) {
    return(next_rho);
  }
  return(rho);
}

vec sample_three_utility_probit_matched_alpha(
    vec y_star_m_1, vec y_star_m_3,  
    vec beta_v, vec delta_v,
    vec alpha_mean_v, mat alpha_cov_s,
    vec delta_mean_v, mat delta_cov_s) {
  
  vec beta_diff_v_1 = beta_v - delta_v(0);
  vec beta_diff_v_2 = beta_v - delta_v(1);
    
  mat post_cov = alpha_cov_s.i();
  post_cov(0, 0) += dot(beta_diff_v_1, beta_diff_v_1);
  post_cov(1, 1) += dot(beta_diff_v_2, beta_diff_v_2);
  
  vec post_mean = solve(alpha_cov_s, alpha_mean_v);
  post_mean(0) -= dot(beta_diff_v_1, y_star_m_1);
  post_mean(1) -= dot(beta_diff_v_2, y_star_m_3);
  post_mean = solve(post_cov, post_mean);
    
  double sample_order_up_prob = 
    R::pnorm(0, post_mean(0), sqrt(1.0 / post_cov(0,0)), false, true) +
    R::pnorm(0, post_mean(1), sqrt(1.0 / post_cov(1,1)), true, true) +
    as_scalar(dmvnorm(delta_v.t(), delta_mean_v, delta_cov_s, true));
  double sample_order_down_prob = 
    R::pnorm(0, post_mean(0), sqrt(1.0 / post_cov(0,0)), true, true) +
    R::pnorm(0, post_mean(1), sqrt(1.0 / post_cov(1,1)), false, true) +
    as_scalar(dmvnorm(delta_v.t(), -delta_mean_v, delta_cov_s, true));
  
  double log_sample_prob = sample_order_up_prob - 
    (max(sample_order_up_prob, sample_order_down_prob) +
    log(1 + exp(min(sample_order_up_prob, sample_order_down_prob) - 
                  max(sample_order_up_prob, sample_order_down_prob))));
  double match_var = (log(randu()) < log_sample_prob) * 2 - 1;
    
  vec out_v(3);
  if (match_var == 1) {
    out_v(0) = rtn1(post_mean(0), 1.0 / sqrt(post_cov(0, 0)), 
                    0, datum::inf);
    out_v(1) = rtn1(post_mean(1), 1.0 / sqrt(post_cov(1, 1)), 
                    -datum::inf, 0);
  } else {
    out_v(0) = rtn1(post_mean(0), 1.0 / sqrt(post_cov(0, 0)), 
                    -datum::inf, 0);
    out_v(1) = rtn1(post_mean(1), 1.0 / sqrt(post_cov(1, 1)), 
                    0, datum::inf);
  }
  out_v(2) = match_var;
  
  return(out_v);
}

vec sample_three_utility_probit_matched_delta(
    vec y_star_m_1, vec y_star_m_3, 
    vec alpha_v, vec beta_v, double match_var,
    vec delta_mean_v, mat delta_cov_s) {
  
  y_star_m_1 += alpha_v(0) * beta_v;
  y_star_m_3 += alpha_v(1) * beta_v;
  
  mat post_cov = beta_v.n_elem * 
    diagmat(alpha_v) * diagmat(alpha_v) + 
    delta_cov_s.i();
  vec post_mean = match_var * solve(delta_cov_s, delta_mean_v);
  post_mean(0) += accu(alpha_v(0) * y_star_m_1);
  post_mean(1) += accu(alpha_v(1) * y_star_m_3);
  return(rmvnorm(1, solve(post_cov, post_mean),
                 post_cov.i()).t());
}

vec sample_y_star_m_na(double mean_m_1, double mean_m_2) {
  vec out_v(3, fill::randn);
  out_v(0) -= mean_m_1;
  out_v(2) -= mean_m_2;
  return(out_v);
}

vec sample_y_star_m_yea(vec y_star_yea, double mean_m_1, double mean_m_2) {

  y_star_yea(0) = 
    rtn1(-mean_m_1, 1, -datum::inf, y_star_yea(1));
  y_star_yea(1) = 
    rtn1(0, 1, max(y_star_yea(0), y_star_yea(2)), datum::inf);
  y_star_yea(2) = 
    rtn1(-mean_m_2, 1, -datum::inf, y_star_yea(1));
  return(y_star_yea);
}

vec sample_y_star_m_no(vec y_star_no, double mean_m_1, double mean_m_2) {
  
  if (y_star_no(2) < y_star_no(1)) {
    y_star_no(0) = 
      rtn1(-mean_m_1, 1, y_star_no(1), datum::inf);
  } else {
    y_star_no(0) = randn() - mean_m_1;
  }
  
  y_star_no(1) = 
    rtn1(0, 1, -datum::inf, max(y_star_no(0), y_star_no(2)));
  
  if (y_star_no(0) < y_star_no(1)) {
    y_star_no(2) = 
      rtn1(-mean_m_2, 1, y_star_no(1), datum::inf);
  } else {
    y_star_no(2) = randn() - mean_m_2;  
  }
  return(y_star_no);
}

vec sample_y_star_m(vec y_star_vec, double vote, double alpha_1, double alpha_2,
                    double leg_pos, double delta_1, double delta_2) {
  
  vec out_vec(3);
  double mean_m_1 = alpha_1 * (leg_pos - delta_1);
  double mean_m_2 = alpha_2 * (leg_pos - delta_2);
  if (vote == 1) {
    out_vec = sample_y_star_m_yea(y_star_vec, mean_m_1, mean_m_2);
  } else {
    out_vec = sample_y_star_m_no(y_star_vec, mean_m_1, mean_m_2);
  }
  return(out_vec);
}

// [[Rcpp::export]]
List sample_three_utility_probit(
  mat vote_m, mat all_param_draws, mat y_star_m_1, mat y_star_m_2, mat y_star_m_3,
  int leg_start_ind, int alpha_v_1_start_ind, int alpha_v_2_start_ind, 
  int delta_v_1_start_ind, int delta_v_2_start_ind, 
  double leg_mean, double leg_sd, vec alpha_mean_v, mat alpha_cov_s,
  vec delta_mean_v, mat delta_cov_s, int num_iter, int start_iter, 
  int keep_iter, int pos_ind, int neg_ind, bool sample_beta) {
  
  
  vec current_param_val_v = all_param_draws.row(0).t();
  for (int i = 0; i < num_iter; i++) {
    if (i % 100 == 0) {
      Rcout << i << "\n";
    }
    
    for (int j = 0; j < vote_m.n_rows; j++) {
      for (int k = 0; k < vote_m.n_cols; k++) {
        if (!is_finite(vote_m(j, k))) {
          continue;
        }
        vec y_star_vec = {y_star_m_1(j, k), 
                          y_star_m_2(j, k), 
                          y_star_m_3(j, k)};
        vec out_v = sample_y_star_m(
          y_star_vec, vote_m(j, k), 
          current_param_val_v(alpha_v_1_start_ind + k),
          current_param_val_v(alpha_v_2_start_ind + k),
          current_param_val_v(leg_start_ind + j), 
          current_param_val_v(delta_v_1_start_ind + k), 
          current_param_val_v(delta_v_2_start_ind + k));
        y_star_m_1(j, k) = out_v(0);  
        y_star_m_2(j, k) = out_v(1);
        y_star_m_3(j, k) = out_v(2);
      }
    }
    
    if (sample_beta) {
      for (unsigned int j = 0; j < vote_m.n_rows; j++) {
        uvec current_ind = {j};
        uvec interested_inds = find_finite(vote_m.row(j).t());
        current_param_val_v(leg_start_ind + j) =
          sample_three_utility_probit_beta(
            y_star_m_1.submat(current_ind, interested_inds),
            y_star_m_3.submat(current_ind, interested_inds),
            current_param_val_v(alpha_v_1_start_ind + interested_inds).t(),
            current_param_val_v(alpha_v_2_start_ind + interested_inds).t(),
            current_param_val_v(delta_v_1_start_ind + interested_inds).t(),
            current_param_val_v(delta_v_2_start_ind + interested_inds).t(),
            leg_mean, leg_sd);
      }
    }
    
    vec match_var_v(vote_m.n_cols);
    for (unsigned int j = 0; j < vote_m.n_cols; j++) {
      uvec current_ind = {j};
      uvec interested_inds = find_finite(vote_m.col(j));
      vec delta_v = {current_param_val_v(delta_v_1_start_ind + j),
                     current_param_val_v(delta_v_2_start_ind + j)};
      vec out_v =
        sample_three_utility_probit_matched_alpha(
          y_star_m_1.submat(interested_inds, current_ind), 
          y_star_m_3.submat(interested_inds, current_ind),  
          current_param_val_v(leg_start_ind + interested_inds), 
          delta_v, alpha_mean_v, alpha_cov_s,
          delta_mean_v, delta_cov_s); 
      
      current_param_val_v(alpha_v_1_start_ind + j) = out_v(0);
      current_param_val_v(alpha_v_2_start_ind + j) = out_v(1);
      match_var_v(j) = out_v(2);
    }
    
    for (unsigned int j = 0; j < vote_m.n_cols; j++) {
      uvec current_ind = {j};
      uvec interested_inds = find_finite(vote_m.col(j));
      vec alpha_v = {current_param_val_v(alpha_v_1_start_ind + j),
                     current_param_val_v(alpha_v_2_start_ind + j)};
      vec out_v =
        sample_three_utility_probit_matched_delta(
          y_star_m_1.submat(interested_inds, current_ind), 
          y_star_m_3.submat(interested_inds, current_ind),
          alpha_v, current_param_val_v(leg_start_ind + interested_inds), 
          match_var_v(j), delta_mean_v, delta_cov_s); 
      current_param_val_v(delta_v_1_start_ind + j) = out_v(0);
      current_param_val_v(delta_v_2_start_ind + j) = out_v(1);
    }
    
    if (pos_ind > -1 && (current_param_val_v(leg_start_ind + pos_ind) < 0)) {
      current_param_val_v = -current_param_val_v;
    }
    
    if (neg_ind > -1 && pos_ind < 0 && (current_param_val_v(leg_start_ind + neg_ind) > 0)) {
      current_param_val_v = -current_param_val_v;
    }
    
    int post_burn_i = i - start_iter + 1;
    if (i >= start_iter && (fmod(post_burn_i, keep_iter) == 0)) {
      int keep_iter_ind = post_burn_i / keep_iter - 1;
      all_param_draws.row(keep_iter_ind) = current_param_val_v.t();
    }
  }
  
  return(List::create(Named("param_draws") = all_param_draws, 
                      Named("y_star_m_1") = y_star_m_1, 
                      Named("y_star_m_2") = y_star_m_2, 
                      Named("y_star_m_3") = y_star_m_3));
}

vec adjust_all_judge_ideology(
    vec current_param_val_v, 
    uvec judge_start_ind,
    uvec case_year_v, uvec case_judge_year_v,
    int alpha_v_1_start_ind, int alpha_v_2_start_ind,
    int delta_v_1_start_ind, int delta_v_2_start_ind,
    uvec pos_judge_ind, uvec pos_judge_year,
    uvec neg_judge_ind, uvec neg_judge_year) {
  
  
  for (int i = 0; i < pos_judge_ind.n_elem; i++) {
    if (current_param_val_v(pos_judge_ind(i)) < 0) {
      uvec judge_year = find(case_judge_year_v == pos_judge_year(i));
      uvec cases = find(case_year_v == pos_judge_year(i));
      current_param_val_v(judge_year) =
        -current_param_val_v(judge_year);
      current_param_val_v(alpha_v_1_start_ind + cases) = 
        -current_param_val_v(alpha_v_1_start_ind + cases);
      current_param_val_v(alpha_v_2_start_ind + cases) = 
        -current_param_val_v(alpha_v_2_start_ind + cases);
      current_param_val_v(delta_v_1_start_ind + cases) = 
        -current_param_val_v(delta_v_1_start_ind + cases);
      current_param_val_v(delta_v_2_start_ind + cases) = 
        -current_param_val_v(delta_v_2_start_ind + cases);
    }
  }
  for (int i = 0; i < neg_judge_ind.n_elem; i++) {
    if (current_param_val_v(neg_judge_ind(i)) > 0) {
      uvec judge_year = find(case_judge_year_v == neg_judge_year(i));
      uvec cases = find(case_year_v == neg_judge_year(i));
      current_param_val_v(judge_year) =
        -current_param_val_v(judge_year);
      current_param_val_v(alpha_v_1_start_ind + cases) = 
        -current_param_val_v(alpha_v_1_start_ind + cases);
      current_param_val_v(alpha_v_2_start_ind + cases) = 
        -current_param_val_v(alpha_v_2_start_ind + cases);
      current_param_val_v(delta_v_1_start_ind + cases) = 
        -current_param_val_v(delta_v_1_start_ind + cases);
      current_param_val_v(delta_v_2_start_ind + cases) = 
        -current_param_val_v(delta_v_2_start_ind + cases);
    }
  }
  return(current_param_val_v);
}

// [[Rcpp::export]]
List sample_three_utility_probit_gp(
    mat vote_m, mat all_param_draws, mat y_star_m_1, mat y_star_m_2, mat y_star_m_3,
    uvec judge_start_inds, uvec judge_end_inds, uvec case_years, 
    umat case_judge_years_ind_m, uvec judge_year_v,
    int alpha_v_1_start_ind, int alpha_v_2_start_ind, 
    int delta_v_1_start_ind, int delta_v_2_start_ind, int rho_ind,
    vec alpha_mean_v, mat alpha_cov_s, vec delta_mean_v, mat delta_cov_s, 
    double rho_mean,double rho_sigma, double rho_sd, int num_iter, int start_iter, 
    int keep_iter, uvec pos_judge_ind, uvec neg_judge_ind,
    uvec pos_judge_year, uvec neg_judge_year) {
  
  
  vec current_param_val_v = all_param_draws.row(0).t();
  // vec accept_count(zeta_param_start_ind - psi_param_start_ind);
  // accept_count.zeros();
  for (int i = 0; i < num_iter; i++) {
    if (i % 100 == 0) {
      Rcout << i << "\n";
    }
  
    for (int j = 0; j < vote_m.n_rows; j++) {
      for (int k = 0; k < vote_m.n_cols; k++) {
        if (!is_finite(vote_m(j, k))) {
          continue;
        }
        vec y_star_vec = {y_star_m_1(j, k), 
                          y_star_m_2(j, k), 
                          y_star_m_3(j, k)};
        vec out_v = sample_y_star_m(
          y_star_vec, vote_m(j, k), 
          current_param_val_v(alpha_v_1_start_ind + k),
          current_param_val_v(alpha_v_2_start_ind + k),
          current_param_val_v(judge_start_inds(j) + case_judge_years_ind_m(j, k)), 
          current_param_val_v(delta_v_1_start_ind + k), 
          current_param_val_v(delta_v_2_start_ind + k));
        y_star_m_1(j, k) = out_v(0);  
        y_star_m_2(j, k) = out_v(1);
        y_star_m_3(j, k) = out_v(2);
      }
    }
    
    for (unsigned int j = 0; j < vote_m.n_rows; j++) {
      uvec current_ind = {j};
      uvec interested_inds = find_finite(vote_m.row(j).t());
      rowvec y_star_m_1_v = y_star_m_1.row(j);
      rowvec y_star_m_3_v = y_star_m_3.row(j);
      uvec judge_years_v = case_judge_years_ind_m.row(j).t();
      current_param_val_v(span(
          judge_start_inds(j), judge_end_inds(j))) =
        sample_three_utility_probit_beta_gp(
          y_star_m_1.submat(current_ind, interested_inds),
          y_star_m_3.submat(current_ind, interested_inds),
          current_param_val_v(alpha_v_1_start_ind + interested_inds).t(),
          current_param_val_v(alpha_v_2_start_ind + interested_inds).t(),
          current_param_val_v(delta_v_1_start_ind + interested_inds).t(),
          current_param_val_v(delta_v_2_start_ind + interested_inds).t(),
          judge_years_v(interested_inds), current_param_val_v(rho_ind));
    }
    
    vec match_var_v(vote_m.n_cols);
    for (unsigned int j = 0; j < vote_m.n_cols; j++) {
      uvec current_ind = {j};
      uvec interested_inds = find_finite(vote_m.col(j));
      vec delta_v = {current_param_val_v(delta_v_1_start_ind + j),
                     current_param_val_v(delta_v_2_start_ind + j)};
      uvec judge_years_v = case_judge_years_ind_m.col(j);
      vec out_v =
        sample_three_utility_probit_matched_alpha(
          y_star_m_1.submat(interested_inds, current_ind), 
          y_star_m_3.submat(interested_inds, current_ind),  
          current_param_val_v(
            judge_start_inds(interested_inds) + 
            judge_years_v(interested_inds)), 
          delta_v, alpha_mean_v, alpha_cov_s,
          delta_mean_v, delta_cov_s); 
      
      current_param_val_v(alpha_v_1_start_ind + j) = out_v(0);
      current_param_val_v(alpha_v_2_start_ind + j) = out_v(1);
      match_var_v(j) = out_v(2);
    }
    
    for (unsigned int j = 0; j < vote_m.n_cols; j++) {
      uvec current_ind = {j};
      uvec interested_inds = find_finite(vote_m.col(j));
      vec alpha_v = {current_param_val_v(alpha_v_1_start_ind + j),
                     current_param_val_v(alpha_v_2_start_ind + j)};
      uvec judge_years_v = case_judge_years_ind_m.col(j);
      vec out_v =
        sample_three_utility_probit_matched_delta(
          y_star_m_1.submat(interested_inds, current_ind), 
          y_star_m_3.submat(interested_inds, current_ind),
          alpha_v, current_param_val_v(
              judge_start_inds(interested_inds) + 
                judge_years_v(interested_inds)), 
          match_var_v(j), delta_mean_v, delta_cov_s); 
      current_param_val_v(delta_v_1_start_ind + j) = out_v(0);
      current_param_val_v(delta_v_2_start_ind + j) = out_v(1);
    }
    
    if (pos_judge_ind.n_elem > 0 || neg_judge_ind.n_elem > 0) {
      current_param_val_v(span(0, rho_ind - 1)) =
        adjust_all_judge_ideology(
          current_param_val_v(span(0, rho_ind - 1)), 
          judge_start_inds, case_years, judge_year_v,
          alpha_v_1_start_ind, alpha_v_2_start_ind,
          delta_v_1_start_ind, delta_v_2_start_ind,
          pos_judge_ind, pos_judge_year,
          neg_judge_ind, neg_judge_year);
    }
    
    current_param_val_v(rho_ind) = sample_rho_pos_logit_gibbs(
      current_param_val_v(rho_ind), 
      current_param_val_v(span(0, alpha_v_1_start_ind - 1)), 
      judge_start_inds, judge_end_inds, rho_mean, rho_sigma, rho_sd);
    
    int post_burn_i = i - start_iter + 1;
    if (i >= start_iter && (fmod(post_burn_i, keep_iter) == 0)) {
      int keep_iter_ind = post_burn_i / keep_iter - 1;
      all_param_draws.row(keep_iter_ind) = current_param_val_v.t();
    }
  }
  
  return(List::create(Named("param_draws") = all_param_draws, 
                      Named("y_star_m_1") = y_star_m_1, 
                      Named("y_star_m_2") = y_star_m_2, 
                      Named("y_star_m_3") = y_star_m_3));
}

double phid(double x) {
  return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}

// BVND calculates the probability that X > DH and Y > DK.
// Note: Prob( X < DH, Y < DK ) = BVND( -DH, -DK, R )
// Code and description is adopted from tvpack.f in the 
// mvtnorm package with help from ChatGPT
// [[Rcpp::export]]
double bvnd(double DH, double DK, double R) {
  
  vec x;
  vec w;
  // double as = 0.0;
  // double a = 0.0;
  double b = 0.0;
  // double c = 0.0;
  // double d = 0.0;
  double rs = 0.0;
  double xs = 0.0;
  double bvn = 0.0;
  // double sn = 0.0;
  // double asr = 0.0;
  double h = DH;
  double k = DK;
  double hk = h * k;
  
  if (std::abs(R) < 0.3) {
    x = {-0.9324695142031522,-0.6612093864662647,
         -0.2386191860831970};
    w = {0.1713244923791705, 
         0.3607615730481384, 0.4679139345726904};
  } else if (std::abs(R) < 0.75) {
    x = {-0.9815606342467191, -0.9041172563704750,
         -0.7699026741943050, -0.5873179542866171,
         -0.3678314989981802, -0.1252334085114692};
    w = {0.4717533638651177e-01, 0.1069393259953183, 
         0.1600783285433464, 0.2031674267230659,
         0.2334925365383547, 0.2491470458134029};
  } else {
    x = {-0.9931285991850949, -0.9639719272779138,
         -0.9122344282513259, -0.8391169718222188,
         -0.7463319064601508, -0.6360536807265150,
         -0.5108670019508271, -0.3737060887154196,
         -0.2277858511416451, -0.7652652113349733e-01};
    w = {0.1761400713915212e-01, 0.4060142980038694e-01, 
         0.6267204833410906e-01, 0.8327674157670475e-01,
         0.1019301198172404, 0.1181945319615184,
         0.1316886384491766, 0.1420961093183821,
         0.1491729864726037, 0.1527533871307259};
  }

  if (std::abs(R) < 0.925) {
    if (std::abs(R) > 0.0) {
      double hs = (h * h + k * k) / 2.0;
      double asr = std::asin(R);
      for (int i = 0; i < x.n_elem; ++i) {
        for (int is = -1; is <= 1; is += 2) {
          double sn = std::sin(asr * (is * x[i] + 1) / 2.0);
          bvn += w[i] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
        }
      }
      bvn = bvn * asr / (2.0 * TWOPI);
    }
    bvn += R::pnorm(-h, 0, 1,-datum::inf, false) * R::pnorm(-k, 0, 1,-datum::inf, false);
  } else {
    if (R < 0.0) {
      k = -k;
      hk = -hk;
    }
    if (std::abs(R) < 1.0) {
      double as = (1.0 - R) * (1.0 + R);
      double a = std::sqrt(as);
      double bs = std::pow(h - k, 2);
      double c = (4.0 - hk) / 8.0;
      double d = (12.0 - hk) / 16.0;
      double asr = -(bs / as + hk) / 2.0;
      if (asr > -100.0) {
        bvn = a * std::exp(asr) * 
          (1.0 - c * (bs - as) * (1.0 - d * bs / 5.0) / 3.0 + 
          c * d * as * as / 5);
      }
      if (-hk < 100) {
        b = sqrt(bs);
        bvn = bvn - exp(-hk/2) * sqrt(TWOPI) * R::pnorm(-b/a, 0, 1,-datum::inf, false) * b
                * (1 - c * bs * (1 - d * bs/5) / 3);
      }
      a = a / 2;
      for (int i = 0; i < x.n_elem; i++) {
        for (int is = -1; is <= 1; is += 2) {
          xs = pow(a * (is * x[i] + 1), 2);
          rs = sqrt(1 - xs);
          asr = -(bs/xs + hk) / 2;
          if (asr > -100) {
            bvn = bvn + a * w[i] * exp(asr)
                   * (exp(-hk * (1 - rs) / (2 * (1 + rs)))/rs
                        - (1 + c * xs * (1 + d * xs)));
            }
          }
        }
        bvn = -bvn/TWOPI;
    }
    if (R > 0) {
      bvn = bvn + R::pnorm(-std::max(h, k), 0, 1,-datum::inf, false);
    } else {
      bvn = -bvn;
      if (k > h) {
        bvn = bvn + R::pnorm(k, 0, 1,-datum::inf, false) - R::pnorm(h, 0, 1,-datum::inf, false);
      }
    }
  }
  return(bvn);
}

// [[Rcpp::export]]
mat calc_probit_bggum_three_utility_post_prob_m(
    mat leg_ideology, mat alpha_m, mat delta_m,
    mat case_vote_m, int num_votes) {
  
  mat post_prob(case_vote_m.n_rows, case_vote_m.n_cols, fill::zeros);
  for (int iter = 0; iter < leg_ideology.n_rows; iter++) {
    for (int j = 0; j < case_vote_m.n_cols; j++) {
      for (int i = 0; i < case_vote_m.n_rows; i++) {
        double mean_1 = 
          alpha_m(iter, 2 * j) * (
              leg_ideology(iter, i) - delta_m(iter, 2 * j));
        double mean_2 = 
          alpha_m(iter, 2 * j + 1) * (
              leg_ideology(iter, i) - delta_m(iter, 2 * j + 1));
        post_prob(i, j) += bvnd(-mean_1 / sqrt(2), -mean_2 / sqrt(2), 0.5);
      }
    }
  }
  return(post_prob);
}

// [[Rcpp::export]]
vec calc_waic_probit_bggum_three_utility(
  mat leg_ideology, mat alpha_m, mat delta_m,
  mat case_vote_m, int num_votes) {
  
  vec mean_prob(num_votes, fill::zeros);
  vec mean_log_prob(num_votes, fill::zeros);
  vec log_prob_var(num_votes, fill::zeros);
  for (int iter = 0; iter < leg_ideology.n_rows; iter++) {
    // if (iter + 1 % 100 == 0) {
    //   Rcout << iter << "\n";
    // }
    Rcout << iter << endl;
    int vote_num = 0;
    // Rcout << vote_num << endl;
    for (int j = 0; j < case_vote_m.n_cols; j++) {
      for (int i = 0; i < case_vote_m.n_rows; i++) {
        if (!is_finite(case_vote_m(i, j))) {
          continue;
        }
        double mean_1 = 
          alpha_m(iter, 2 * j) * (
            leg_ideology(iter, i) - delta_m(iter, 2 * j));
        double mean_2 = 
          alpha_m(iter, 2 * j + 1) * (
              leg_ideology(iter, i) - delta_m(iter, 2 * j + 1));
        double yea_prob = bvnd(-mean_1 / sqrt(2), -mean_2 / sqrt(2), 0.5);
        yea_prob = min(yea_prob, 1 - 1e-9);
        yea_prob = max(yea_prob, 1e-9);
        double log_prob = case_vote_m(i, j) * log(yea_prob) +
          (1 - case_vote_m(i, j)) * log(1 - yea_prob);
        mean_prob(vote_num) += exp(log_prob);
        double next_mean_log_prob = (iter * mean_log_prob(vote_num) + log_prob) / (iter + 1);
        log_prob_var(vote_num) +=
          (log_prob - mean_log_prob(vote_num)) * (log_prob - next_mean_log_prob);
        mean_log_prob(vote_num) = next_mean_log_prob;
        vote_num++;
      }
    }
    // Rcout << vote_num << endl;
  }
  return(
    log(mean_prob / leg_ideology.n_rows) -
      (log_prob_var) / (leg_ideology.n_rows - 1));
}

// [[Rcpp::export]]
vec calc_waic_probit_bggum_three_utility_block(
    mat leg_ideology, mat alpha_m, mat delta_m,
    mat case_vote_m, uvec case_year, mat block_m) {
  
  vec mean_prob(block_m.n_rows);
  mean_prob.fill(-datum::inf);
  vec mean_log_prob(block_m.n_rows, fill::zeros);
  vec log_prob_var(block_m.n_rows, fill::zeros);
  for (int iter = 0; iter < leg_ideology.n_rows; iter++) {
    Rcout << iter << endl;
    for (int ind = 0; ind < block_m.n_rows; ind++) {
      int i = block_m(ind, 0);
      int year = block_m(ind, 1);
      int judge_ind = i + (year - 1) * case_vote_m.n_rows;
      double log_prob = 0;
      uvec interested_cases = find(case_year == year);
      Rcout << interested_cases << endl;
      for (int j : interested_cases) {
        if (!is_finite(case_vote_m(i, j))) {
          continue;
        }
        double mean_1 = 
          alpha_m(iter, 2 * j) * (
              leg_ideology(iter, judge_ind) - delta_m(iter, 2 * j));
        double mean_2 = 
          alpha_m(iter, 2 * j + 1) * (
              leg_ideology(iter, judge_ind) - delta_m(iter, 2 * j + 1));
        double yea_prob = bvnd(-mean_1 / sqrt(2), -mean_2 / sqrt(2), 0.5);
        yea_prob = min(yea_prob, 1 - 1e-9);
        yea_prob = max(yea_prob, 1e-9);
        log_prob += case_vote_m(i, j) * log(yea_prob) +
          (1 - case_vote_m(i, j)) * log(1 - yea_prob);
      }
      mean_prob(ind) = max(mean_prob(ind), log_prob) + 
        log(1 + exp(min(mean_prob(ind), log_prob) - max(mean_prob(ind), log_prob)));
      double next_mean_log_prob = (iter * mean_log_prob(ind) + log_prob) / (iter + 1);
      log_prob_var(ind) +=
        (log_prob - mean_log_prob(ind)) * (log_prob - next_mean_log_prob);
      mean_log_prob(ind) = next_mean_log_prob;
    }
  }
  return(
    mean_prob - log(leg_ideology.n_rows) -
      (log_prob_var) / (leg_ideology.n_rows - 1));
}

// [[Rcpp::export]]
vec calc_waic_probit_bggum_three_utility_block_rcpp(
    mat leg_ideology, mat alpha_m, mat delta_m,
    mat case_vote_m, uvec case_year, mat block_m) {
  
  vec mean_prob(block_m.n_rows);
  mean_prob.fill(-datum::inf);
  vec mean_log_prob(block_m.n_rows, fill::zeros);
  vec log_prob_var(block_m.n_rows, fill::zeros);
  for (int iter = 0; iter < leg_ideology.n_rows; iter++) {
    Rcout << iter << endl;
    for (int ind = 0; ind < block_m.n_rows; ind++) {
      int i = block_m(ind, 0);
      int year = block_m(ind, 1);
      double log_prob = 0;
      uvec interested_cases = find(case_year == year);
      for (int j : interested_cases) {
        if (!is_finite(case_vote_m(i, j))) {
          continue;
        }
        double mean_1 = 
          alpha_m(iter, 2 * j) * (
              leg_ideology(iter, ind) - delta_m(iter, 2 * j));
        double mean_2 = 
          alpha_m(iter, 2 * j + 1) * (
              leg_ideology(iter, ind) - delta_m(iter, 2 * j + 1));
        double yea_prob = bvnd(-mean_1 / sqrt(2), -mean_2 / sqrt(2), 0.5);
        yea_prob = min(yea_prob, 1 - 1e-9);
        yea_prob = max(yea_prob, 1e-9);
        log_prob += case_vote_m(i, j) * log(yea_prob) +
          (1 - case_vote_m(i, j)) * log(1 - yea_prob);
      }
      mean_prob(ind) = max(mean_prob(ind), log_prob) + 
        log(1 + exp(min(mean_prob(ind), log_prob) - max(mean_prob(ind), log_prob)));
      double next_mean_log_prob = (iter * mean_log_prob(ind) + log_prob) / (iter + 1);
      log_prob_var(ind) +=
        (log_prob - mean_log_prob(ind)) * (log_prob - next_mean_log_prob);
      mean_log_prob(ind) = next_mean_log_prob;
    }
  }
  return(
    mean_prob - log(leg_ideology.n_rows) -
      (log_prob_var) / (leg_ideology.n_rows - 1));
}
