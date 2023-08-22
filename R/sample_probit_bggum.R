library(boot)
library(tidyverse)

init_y_star_m <- function(vote_m) {
  y_star_m_1 <- vote_m
  y_star_m_2 <- vote_m
  y_star_m_3 <- vote_m

  y_star_m_1[which(vote_m == 1)] = 0
  y_star_m_2[which(vote_m == 1)] = 1
  y_star_m_3[which(vote_m == 1)] = 0

  no_votes <- which(vote_m == 0)
  sample_upper <- rbernoulli(length(no_votes))
  y_star_m_1[no_votes[which(sample_upper == 1)]] = 0
  y_star_m_2[no_votes[which(sample_upper == 1)]] = 0
  y_star_m_3[no_votes[which(sample_upper == 1)]] = 1

  y_star_m_1[no_votes[which(sample_upper == 0)]] = 1
  y_star_m_2[no_votes[which(sample_upper == 0)]] = 0
  y_star_m_3[no_votes[which(sample_upper == 0)]] = 0

  return(list(y_star_m_1, y_star_m_2, y_star_m_3))
}

init_data_rcpp <- function(vote_m, leg_pos_init, alpha_pos_init, delta_pos_init, 
  y_star_m_1_init, y_star_m_2_init, y_star_m_3_init, total_iter) {
    
    if (!is.null(leg_pos_init)) {
      leg_pos_m <- 
        matrix(leg_pos_init, nrow = total_iter, ncol = nrow(vote_m), byrow = T)
    } else {
      leg_pos_m <- matrix(0, nrow = total_iter, ncol = nrow(vote_m))
    }
    
    if (!is.null(alpha_pos_init)) {
      alpha_pos_m <- 
        matrix(t(alpha_pos_init), nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
    } else {
      alpha_pos_m <- 
        matrix(rep(c(-1, 1), ncol(vote_m)), 
               nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
    }
    
    if (!is.null(delta_pos_init)) {
      delta_pos_m <- 
        matrix(t(delta_pos_init), nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
    } else {
      delta_pos_m <- 
        matrix(0, nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
    }
    
    if (!is.null(y_star_m_1_init)) {
      y_star_m_1 <- y_star_m_1_init
      y_star_m_2 <- y_star_m_2_init
      y_star_m_3 <- y_star_m_3_init
    } else {
      y_star_info <- init_y_star_m(vote_m)
      y_star_m_1 <- y_star_info[[1]]
      y_star_m_2 <- y_star_info[[2]]
      y_star_m_3 <- y_star_info[[3]]
    }
    
    all_params_draw <- cbind(leg_pos_m, alpha_pos_m, delta_pos_m)
    beta_start_ind = 0;
    alpha_start_ind = nrow(vote_m);
    alpha_2_start_ind = alpha_start_ind + ncol(vote_m);
    delta_start_ind = alpha_2_start_ind + ncol(vote_m);
    delta_2_start_ind = delta_start_ind + ncol(vote_m);
  
    return(list(all_params_draw, y_star_m_1, y_star_m_2, 
                y_star_m_3, beta_start_ind,
                alpha_start_ind, alpha_2_start_ind,
                delta_start_ind, delta_2_start_ind))
}

sample_three_utility_probit_rcpp <- function(
  vote_m, leg_mean, leg_s, alpha_mean, alpha_cov_s,
  delta_mean, delta_cov_s, num_iter = 2000, start_iter = 0, keep_iter = 1,
  leg_pos_init = NULL, alpha_pos_init = NULL, delta_pos_init = NULL,
  y_star_m_1_init = NULL, y_star_m_2_init = NULL, y_star_m_3_init = NULL,
  pos_ind = 0, neg_ind = 0, start_val = NULL, sample_beta = T) {
  
  total_iter = (num_iter - start_iter) %/% keep_iter
  init_info <- init_data_rcpp(
    vote_m, leg_pos_init, alpha_pos_init, delta_pos_init, 
    y_star_m_1_init, y_star_m_2_init, y_star_m_3_init, total_iter)
  
  if (!is.null(start_val)) {
    init_info[[1]][1,] <- start_val
  }
  draw_info <- sample_three_utility_probit(
    vote_m, init_info[[1]], init_info[[2]], init_info[[3]], init_info[[4]],
    init_info[[5]], init_info[[6]], init_info[[7]], 
    init_info[[8]], init_info[[9]], leg_mean, leg_s, 
    alpha_mean, alpha_cov_s, delta_mean, delta_cov_s,
    num_iter, start_iter, keep_iter, pos_ind - 1, neg_ind - 1, sample_beta)
  
  all_param_draw = draw_info[[1]]
  leg_names <- sapply(rownames(vote_m), function(name) {paste(name, "beta", sep = "_")})
  if (is.null(colnames(vote_m))) {
    colnames(vote_m) <- sapply(1:ncol(vote_m), function(i) {
      paste("vote", i, sep = "_")
    })
  }
  alpha_vote_names_1 <- sapply(colnames(vote_m), function(name) {
    paste(name, "alpha", "1", sep = "_")
  })
  alpha_vote_names_2 <- sapply(colnames(vote_m), function(name) {
    paste(name, "alpha", "2", sep = "_")
  })
  delta_vote_names_1 <- sapply(colnames(vote_m), function(name) {
    paste(name, "delta", "1", sep = "_")
  })
  delta_vote_names_2 <- sapply(colnames(vote_m), function(name) {
    paste(name, "delta", "2", sep = "_")
  })
  colnames(all_param_draw) <- 
    c(leg_names, alpha_vote_names_1, alpha_vote_names_2, delta_vote_names_1, delta_vote_names_2)
  
  return(c(list("param_draws" = all_param_draw), draw_info[-1]))
}

init_data_gp_rcpp <- function(
  vote_m, years_v, leg_pos_init, alpha_pos_init, delta_pos_init, rho_init, 
  y_star_m_1_init, y_star_m_2_init, y_star_m_3_init, total_iter,
  pos_ind_list, neg_ind_list, pos_ind_years_list, neg_ind_years_list) {  
  
  judge_start_inds <- rep(0, nrow(vote_m))
  judge_end_inds <- rep(0, nrow(vote_m))

  case_judge_years_ind_m <- vote_m
  years_served <- rep(0, nrow(vote_m))
  judge_year_v <- c()
  judge_start_inds <- c() 
  judge_end_inds <- c()
  start_ind = 0
  end_ind = 0
  for (i in 1:nrow(vote_m)) {
    interested_inds <- which(!is.na(vote_m[i,]))
    case_judge_years_ind_m[i,interested_inds] <-
      years_v[interested_inds] - min(years_v[interested_inds])
    judge_year_v <- 
      c(judge_year_v, list(min(years_v[interested_inds]):max(years_v[interested_inds])))
    end_ind = start_ind + 
      max(years_v[interested_inds]) - min(years_v[interested_inds])
    judge_start_inds <- c(judge_start_inds, start_ind)
    judge_end_inds <- c(judge_end_inds, end_ind)
    start_ind = end_ind + 1
  }  
  
  if (!is.null(leg_pos_init)) {
    leg_pos_m <- 
      matrix(leg_pos_init, nrow = total_iter, ncol = length(leg_pos_init), byrow = T)
  } else {
    leg_pos_m <- matrix(0, nrow = total_iter, ncol = max(judge_end_inds) + 1)
  }
  
  if (!is.null(alpha_pos_init)) {
    alpha_pos_m <- 
      matrix(t(alpha_pos_init), nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
  } else {
    alpha_pos_m <- 
      matrix(rep(c(-1, 1), ncol(vote_m)), 
             nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
  }
  
  if (!is.null(delta_pos_init)) {
    delta_pos_m <- 
      matrix(t(delta_pos_init), nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
  } else {
    delta_pos_m <- 
      matrix(0, nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
  }
  
  if (!is.null(rho_init)) {
    rho_v <- rep(rho_init, nrow = total_iter)
  } else {
    rho_v <- rep(0.9, nrow = total_iter)
  }
  
  if (!is.null(y_star_m_1_init)) {
    y_star_m_1 <- y_star_m_1_init
    y_star_m_2 <- y_star_m_2_init
    y_star_m_3 <- y_star_m_3_init
  } else {
    y_star_info <- init_y_star_m(vote_m)
    y_star_m_1 <- y_star_info[[1]]
    y_star_m_2 <- y_star_info[[2]]
    y_star_m_3 <- y_star_info[[3]]
  }
  
  all_params_draw <- cbind(leg_pos_m, alpha_pos_m, delta_pos_m, rho_v)
  beta_start_ind = 0;
  alpha_start_ind = max(judge_end_inds) + 1;
  alpha_2_start_ind = alpha_start_ind + ncol(vote_m);
  delta_start_ind = alpha_2_start_ind + ncol(vote_m);
  delta_2_start_ind = delta_start_ind + ncol(vote_m);
  rho_start_ind = delta_2_start_ind + ncol(vote_m);
  
  pos_ind_judge_list <- vector(mode = "integer")
  pos_ind_judge_year_list <- vector(mode = "integer")
  if (length(pos_ind_list) > 0) {
    for (i in 1:length(pos_ind_list)) {
      tmp_judge_list <- 
        judge_start_inds[pos_ind_list[i]]:
        judge_end_inds[pos_ind_list[i]]
      tmp_judge_year_list <- judge_year_v[[pos_ind_list[i]]]
      if (length(pos_ind_years_list) > 0) {
        tmp_judge_list <- tmp_judge_list[pos_ind_years_list[[i]]]
        tmp_judge_year_list <-
          tmp_judge_year_list[pos_ind_years_list[[i]]]
      }
      pos_ind_judge_list <- c(pos_ind_judge_list, tmp_judge_list)
      pos_ind_judge_year_list <- 
        c(pos_ind_judge_year_list, tmp_judge_year_list)
    } 
  }
  neg_ind_judge_list <- vector(mode = "integer")
  neg_ind_judge_year_list <- vector(mode = "integer")
  if (length(neg_ind_list) > 0) {
    for (i in 1:length(neg_ind_list)) {
      tmp_judge_list <- 
        judge_start_inds[neg_ind_list[i]]:
        judge_end_inds[neg_ind_list[i]]
      tmp_judge_year_list <- judge_year_v[[neg_ind_list[i]]]
      if (length(neg_ind_years_list) > 0) {
        tmp_judge_list <- tmp_judge_list[neg_ind_years_list[[i]]]
        tmp_judge_year_list <-
          tmp_judge_year_list[neg_ind_years_list[[i]]]
      } 
      neg_ind_judge_list <- c(neg_ind_judge_list, tmp_judge_list)
      neg_ind_judge_year_list <- 
        c(neg_ind_judge_year_list, tmp_judge_year_list)
    } 
  }
  
  return(list(all_params_draw, y_star_m_1, y_star_m_2, 
              y_star_m_3, beta_start_ind,
              alpha_start_ind, alpha_2_start_ind,
              delta_start_ind, delta_2_start_ind, rho_start_ind,
              judge_start_inds, judge_end_inds, 
              case_judge_years_ind_m, judge_year_v,
              pos_ind_judge_list, neg_ind_judge_list,
              pos_ind_judge_year_list, neg_ind_judge_year_list))
}

sample_three_utility_probit_gp_rcpp <- function(
  vote_m, years_v, leg_mean, leg_s, alpha_mean, alpha_cov_s,
  delta_mean, delta_cov_s, rho_mean, rho_sigma, rho_sd = 0.1,
  num_iter = 2000, start_iter = 0, keep_iter = 1,
  leg_pos_init = NULL, alpha_pos_init = NULL, delta_pos_init = NULL,
  rho_init = NULL, y_star_m_1_init = NULL, y_star_m_2_init = NULL, y_star_m_3_init = NULL,
  pos_ind = 0, neg_ind = 0,  pos_ind_list = NULL, neg_ind_list  = NULL, 
  pos_ind_years_list = NULL, neg_ind_years_list = NULL,
  start_val = NULL) {
  
  total_iter = (num_iter - start_iter) %/% keep_iter
  init_info <- init_data_gp_rcpp(
    vote_m, years_v, leg_pos_init, alpha_pos_init, delta_pos_init, rho_init,
    y_star_m_1_init, y_star_m_2_init, y_star_m_3_init, total_iter, 
    pos_ind_list, neg_ind_list, pos_ind_years_list, neg_ind_years_list)
  
  if (!is.null(start_val)) {
    init_info[[1]][1,] <- start_val
  }
  
  #Init info
  draw_info <- sample_three_utility_probit_gp(
    vote_m, init_info[[1]], init_info[[2]], init_info[[3]], init_info[[4]],
    init_info[[11]], init_info[[12]], years_v,
    init_info[[13]], unlist(init_info[[14]]), init_info[[6]], init_info[[7]], 
    init_info[[8]], init_info[[9]], init_info[[10]], 
    alpha_mean, alpha_cov_s, delta_mean, delta_cov_s,
    rho_mean, rho_sigma, rho_sd,
    num_iter, start_iter, keep_iter, init_info[[15]], init_info[[16]],
    init_info[[17]], init_info[[18]])
  
  all_param_draw = draw_info[[1]]
  leg_names <- unlist(sapply(1:nrow(vote_m), function(i) {
    sapply(init_info[[14]][[i]], function(year) {
      paste(rownames(vote_m)[i], "beta", year, sep = "_")
    })
  }))
  
  if (is.null(colnames(vote_m))) {
    colnames(vote_m) <- sapply(1:ncol(vote_m), function(i) {
      paste("vote", i, sep = "_")
    })
  }
  alpha_vote_names_1 <- sapply(colnames(vote_m), function(name) {
    paste(name, "alpha", "1", sep = "_")
  })
  alpha_vote_names_2 <- sapply(colnames(vote_m), function(name) {
    paste(name, "alpha", "2", sep = "_")
  })
  delta_vote_names_1 <- sapply(colnames(vote_m), function(name) {
    paste(name, "delta", "1", sep = "_")
  })
  delta_vote_names_2 <- sapply(colnames(vote_m), function(name) {
    paste(name, "delta", "2", sep = "_")
  })
  colnames(all_param_draw) <- 
    c(leg_names, alpha_vote_names_1, alpha_vote_names_2, 
      delta_vote_names_1, delta_vote_names_2, "rho")
  
  return(c(list("param_draws" = all_param_draw), draw_info[-1]))
}


