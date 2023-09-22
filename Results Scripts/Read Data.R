## Read in Results Data ##

require(dplyr)
require(purrr)

sim_results_dat_list <- list.files(path = "C:/Users/gwbuh/Desktop/Non-parametric Condor Results/Sim Results 9-20-23", pattern = ".rds", full.names = TRUE) 


sim_results_dat_1 <- map_dfr(sim_results_dat_list[1:9000], readRDS)
gc()
sim_results_dat_2 <- map_dfr(sim_results_dat_list[8001:16000], readRDS)
gc()
sim_results_dat_3 <- map_dfr(sim_results_dat_list[16001:24000], readRDS)
gc()

sim_results_dat <- rbind(sim_results_dat_1)

saveRDS(sim_results_dat_1, "Raw Sim Data/sim_results_dat_9-20-23.rds")
