library(MCMCpack)

theta_constraints_list <- 
  list("LOTT (R MS-5)" = "+",
       "GINGRICH (R GA-6)" = "+",
       "GINGRICH (R GA-6)" = "+",
       "GINGRICH (R GA-6)" = "+",
       "DELAY (R TX-22)" = "+",
       "DELAY (R TX-22)" = "+",
       "DELAY (R TX-22)" = "+",
       "DELAY (R TX-22)" = "+",
       "BLUNT (R MO-7)" = "+",
       "BLUNT (R MO-7)" = "+",
       "BLUNT (R MO-7)" = "+",
       "CANTOR (R VA-7)" = "+",
       "MCCARTHY (R CA-22)" = "+",
       "SCALISE (R LA-1)" = "+",
       "SCALISE (R LA-1)" = "+",
       "SCALISE (R LA-1)" = "+",
       "SCALISE (R LA-1)" = "+",
       "SCALISE (R LA-1)" = "+")

for (house_num in 100:117) {
  
  house_ind = house_num - 100 + 1
  load(paste("data_files/processed_house_votes_", house_num, ".Rdata", sep = ""))
  
  chain_run <- MCMCirt1d(
    datamatrix = house_votes_m, 
    theta.constraints = theta_constraints_list[house_ind],
    T0 = 1/5, AB0 = 1/5, burnin = 10000, mcmc = 20000, verbose = 1000, 
    store.item=T)
  save(processed_chains, file = 
         paste("result_files/house", house_num, "results_MQ.Rdata",
               sep ="_"))
}
