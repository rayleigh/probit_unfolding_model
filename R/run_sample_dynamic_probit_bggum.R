library(Rcpp)

source("sample_probit_bggum.R")
sourceCpp("three_utility_probit_helper_functions.cpp")

args = commandArgs(trailingOnly=TRUE)
data_file = args[1]
result_file = args[2]
#init_file = args[3]

load("data_files/mq_supreme_court_vote_info_2021.Rdata")

chain_run <- sample_three_utility_probit_gp_rcpp(
  mqVotes, mqTime, 0, 1, c(0, 0), 25 * diag(2), 
  c(-2, 10), 10 * diag(2), 0.9, 0.04, 0.1, 
  pos_ind_list = pos_inds, neg_ind_list = neg_inds,
  neg_ind_years_list = neg_year_inds,
  num_iter = 500000, start_iter = 300000, keep_iter = 10)
run_time <- proc.time() - start_time
warnings()
save(chain_run, run_time, file = 
       "result_files/mq_supreme_court_results_2021_update_dynamic_bbgum_three_utility_rcpp.Rdata")
