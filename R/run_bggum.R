library(bggum)

theta_constraints_list <-
  c("LOTT (R MS-5)",
       "GINGRICH (R GA-6)",
       "GINGRICH (R GA-6)",
       "GINGRICH (R GA-6)",
       "DELAY (R TX-22)",
       "DELAY (R TX-22)",
       "DELAY (R TX-22)",
       "DELAY (R TX-22)",
       "BLUNT (R MO-7)",
       "BLUNT (R MO-7)",
       "BLUNT (R MO-7)",
       "CANTOR (R VA-7)",
       "MCCARTHY (R CA-22)",
       "SCALISE (R LA-1)",
       "SCALISE (R LA-1)",
       "SCALISE (R LA-1)",
       "SCALISE (R LA-1)",
       "SCALISE (R LA-1)")

for (house_num in 100:117) {
  
  house_ind = house_num - 100 + 1
  load(paste("data_files/processed_house_votes_", house_num, ".Rdata", sep = ""))
  
  #set.seed(123)
  proposal_sds <- tune_proposals(house_votes_m, tune_iterations = 5000)
  #set.seed(456)
  temps <- tune_temperatures(house_votes_m, n_temps = 6, 
                             proposal_sds = proposal_sds)
  round(temps, 2)
  chain_run <- ggumMC3(data = house_votes_m,
                       sample_iterations = 20000,
                       burn_iterations = 5000,
                       proposal_sds = proposal_sds,
                       temps = temps)
  constraint <- which(rownames(house_votes_m) == 
                        theta_constraints_list[house_ind])
  processed_chains <- post_process(chain_run, constraint = constraint, 
                                   expected_sign = "+")
  save(processed_chains, file = 
    paste("result_files/result_files/house", house_num, "results_bbgum.Rdata",
          sep = "_"))
}

