args <- commandArgs(trailingOnly = TRUE)

cat(args, "\n")

# Parse Command Line Arg --------------------------------------------------

batch <- 
  args |>
  paste(collapse = " ") |>
  stringr::str_extract("batch [0-9]+") |>
  stringr::str_sub(7, -1) |>
  as.integer()

cat(batch, "\n")

# Source Packages and Functions -------------------------------------------

library(dplyr)       # CRAN
library(tidyr)       # CRAN
library(purrr)       # CRAN
library(WeightIt)    # CRAN
library(lme4)        # CRAN
library(BART)        # CRAN
library(dbarts)      # CRAN
library(stan4bart)   # github
library(bcf)         # github
library(multibart)   # tar.gz
library(arm)         # CRAN
library(caret)       # CRAN
library(performance) # CRAN
library(grf)         # CRAN

source("Condor Scripts/nonparametric-simulation-functions.R")

# Create Design -----------------------------------------------------------

design_factors <- list(beta_0 = 20, beta_1 = 1.5, beta_2 = 0.2, beta_3 = 0.01, beta_4 = 1.5,
                       beta_5 = 0.2, beta_6 = -0.25, beta_7 = -5, beta_8 = 0, beta_9 = 1,
                       beta_10 = 1.5, beta_11 = 2.5, beta_12 = 2, beta_13 = 0, beta_14 = -1, 
                       beta_15 = -1.5, beta_16 = -2, beta_17 = -5, beta_18 = 0, 
                       beta_19 = -1.75, beta_20 = -2, beta_21 = -3.5, beta_22 = -5,
                       trt_ceiling = 1.50, RCT = c(TRUE, FALSE), prop_PRIV = 0.05, ICC = c(0.05, 0.10, 0.15), J = c(50, 100, 150, 200), 
                       submit = seq(1:1000), treatment = "trt", testing = FALSE)

params <- 
  expand.grid(design_factors) |>
  dplyr::mutate(iterations = 1,
                seed = 20230728L + 1:n()) %>%
  dplyr::mutate(treatment = as.character(treatment))

# Run Simulations for Specified Batch -------------------------------------

res <- 
  do.call(run_sim, args = as.list(params[batch,])) |>
  dplyr::cross_join(params[batch,])

# Save Results ------------------------------------------------------------

saveRDS(res, file = paste0("nonparametric-project/results/nonparametric_results_batch", batch, ".rds"))



