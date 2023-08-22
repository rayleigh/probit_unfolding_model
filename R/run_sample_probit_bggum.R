library(Rcpp)

source("code/R/sample_probit_bggum.R")
sourceCpp("code/c++/three_utility_probit_helper_functions.cpp")

# args = commandArgs(trailingOnly=TRUE)
# data_file = args[1]
# result_file = args[2]
# init_file = args[3]
# house_ind = as.numeric(args[4])

load(data_file)

theta_constraints_list <-
  c("LOTT \\(R MS-5\\)",
       "GINGRICH \\(R GA-6\\)",
       "GINGRICH \\(R GA-6\\)",
       "GINGRICH \\(R GA-6\\)",
       "DELAY \\(R TX-22\\)",
       "DELAY \\(R TX-22\\)",
       "DELAY \\(R TX-22\\)",
       "DELAY \\(R TX-22\\)",
       "BLUNT \\(R MO-7\\)",
       "BLUNT \\(R MO-7\\)",
       "BLUNT \\(R MO-7\\)",
       "CANTOR \\(R VA-7\\)",
       "MCCARTHY \\(R CA-22\\)",
       "SCALISE \\(R LA-1\\)",
       "SCALISE \\(R LA-1\\)",
       "SCALISE \\(R LA-1\\)",
       "SCALISE \\(R LA-1\\)",
       "SCALISE \\(R LA-1\\)")

for (house_num in 100:117) {
  
  load(paste("data_files/processed_house_votes_", house_num, ".Rdata", sep = ""))
  house_ind = house_num - 100 + 1
  
  load(paste("result_files/house", house_num, "results_MQ.Rdata",
             sep ="_"))
  leg_inits <- colMeans(chain_run[,grep("theta", colnames(chain_run))])
  mq_start_three_up <- 
    sample_three_utility_probit_rcpp(
      house_votes_m, 0, 1, c(0, 0), diag(2) * 25, 
      c(-2, 10), diag(2) * 10, 
      num_iter = 20, start_iter = 0, keep_iter = 1,
      leg_pos_init = leg_inits)
  rm(chain_run)
  
  start_time <- proc.time()
  chain_run <- sample_three_utility_probit_rcpp(
    house_votes_m, 0, 1, c(0, 0), diag(2) * 25, 
    c(-2, 10), diag(2) * 10, 
    num_iter = 400000, start_iter = 200000, keep_iter = 10,
    pos_ind = grep(theta_constraints_list[house_ind], 
                   rownames(house_votes_m)),
    y_star_m_1_init = mq_start_three_up[[2]], 
    y_star_m_2_init = mq_start_three_up[[3]], 
    y_star_m_3_init = mq_start_three_up[[4]],
    start_val = mq_start_three_up[[1]][nrow(mq_start_three_up[[1]]),])
  run_time <- proc.time() - start_time

  save(chain_run, run_time, file = 
         paste("result_files/house", house_num, "results_bbgum_three_utility_mixture_flipped_rcpp.Rdata",
               sep ="_"))  
}

