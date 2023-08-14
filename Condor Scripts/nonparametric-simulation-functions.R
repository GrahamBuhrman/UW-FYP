### Simulation Functions ###


# Packages ----------------------------------------------------------------

require(dplyr)
require(tidyr)
require(WeightIt)
require(lme4)
require(BART)
require(dbarts)
require(stan4bart)
require(bcf)
require(multibart)
require(arm)
require(caret)
require(grf)
require(performance)

# Load Covariates ---------------------------------------------------------

# load TIMSS covariate space
load_covariates <- function() {
  
  # get data
  dat <- readRDS("Simulation Data/TIMSS_covariates.rds")
  
  # process variables
  dat$IDSCH <- as.factor(dat$IDSCH)
  dat$SEX <- dat$STUSEX_1
  dat$PRIVATE <- dat$PRIV_1
  dat$ABSENT_onceweek <- dat$ABSENT_1
  dat$ABSENT_oncetwoweeks <- dat$ABSENT_2
  dat$ABSENT_oncemonth <- dat$ABSENT_3
  dat$ABSENT_oncetwomonths <- dat$ABSENT_4
  dat$ABSENT_never <- dat$ABSENT_5
  dat$PARED_univ <- dat$PARENTED_1
  dat$PARED_postsecondary <- dat$PARENTED_2
  dat$PARED_uppersecondary <- dat$PARENTED_3
  dat$PARED_lowersecondary <- dat$PARENTED_4
  dat$PARED_someprimary <- dat$PARENTED_5
  dat$FREPL_tenth <- dat$FREPL_0
  dat$FREPL_quarter <- dat$FREPL_1
  dat$FREPL_half <- dat$FREPL_2
  dat$FREPL_threequarter <- dat$FREPL_3
  dat$FREPL_overthreequarter <- dat$FREPL_4
  
  # organize variables
  covs_cat_final <- c("STUSEX", "PRIV", "ABSENT", "FREPL", "PARENTED")
  covs_cont <- c("HOMERES", "BELONG", "LIKEMATH", "CONFMATH", "ACADSUC")
  covs_cat_gen <- c("SEX", "PRIVATE", "ABSENT_onceweek", "ABSENT_oncetwoweeks", 
                    "ABSENT_oncemonth", "ABSENT_oncetwomonths", "ABSENT_never",
                    "PARED_univ", "PARED_postsecondary", "PARED_uppersecondary", 
                    "PARED_lowersecondary", "PARED_someprimary", "FREPL_tenth",
                    "FREPL_quarter", "FREPL_half", "FREPL_threequarter", "FREPL_overthreequarter")
  
  # select variables 
  dat <- dat[, c("IDSTU", "IDSCH", covs_cat_final, covs_cat_gen, covs_cont)]
  
  # return data
  return(dat)
  
}

# Data Generation ---------------------------------------------------------

# define data generation mechanism
generate_data <- function(beta_vec = c(20, 1.5, 0.2, 0.01, 1.5, 0.2, -0.25, -5, 0, 1, 1.5, 2.5, 2, 0, -1, -1.5, -2, -5, 0, -1.75, -2, -3.5, -5), 
                          trt_ceiling = 1.50, RCT = FALSE, prop_PRIV = 0.05, sigma_2 = 0.9, tau_2 = 0.1, J = 100, testing = FALSE) {
  
  # load covariates
  dat <- load_covariates()
  
  # organize variables
  covs_cat_final <- c("STUSEX", "PRIV", "ABSENT", "FREPL", "PARENTED")
  covs_cont <- c("HOMERES", "BELONG", "LIKEMATH", "CONFMATH", "ACADSUC")
  covs_cat_gen <- c("SEX", "PRIVATE", "ABSENT_onceweek", "ABSENT_oncetwoweeks", 
                    "ABSENT_oncemonth", "ABSENT_oncetwomonths", "ABSENT_never",
                    "PARED_univ", "PARED_postsecondary", "PARED_uppersecondary", 
                    "PARED_lowersecondary", "PARED_someprimary", "FREPL_tenth",
                    "FREPL_quarter", "FREPL_half", "FREPL_threequarter", "FREPL_overthreequarter")
  
  # run some checks
  stopifnot(J <= length(unique(dat$IDSCH)))
  stopifnot(floor(J * prop_PRIV) <= 14)
  
  # sample clusters
  tab <- subset(data.frame(table(dat$IDSCH, dat$PRIVATE)), Freq != 0) 
  tab <- tab[, c("Var1", "Var2")]
  rownames(tab) <- NULL
  priv_schools <- subset(tab, Var2 == 1)$Var1
  pub_schools <- subset(tab, Var2 == 0)$Var1
  J_priv <- floor(J * prop_PRIV)
  J_pub <-  J - J_priv
  priv_sample <- sample(priv_schools, size = J_priv, replace = TRUE)
  pub_sample <- sample(pub_schools, size = J_pub, replace = TRUE)
  clusters <- c(priv_sample, pub_sample)
  dat_sampled <- subset(dat, IDSCH %in% clusters)
  
  # selection model
  if(isTRUE(RCT)) {
    
    dat_sampled$trt <- rbinom(nrow(dat_sampled), 1, 0.5)
    
  } else {
    
    trt_assign_priv <- ifelse(dat_sampled$PRIVATE == 1, 1, 0)
    trt_assign_absent <- ifelse(dat_sampled$ABSENT_never == 1 | dat_sampled$ABSENT_oncetwomonths == 1, 1, 0)
    trt_assign_vec <- -1.80 + (1.10 * trt_assign_priv) + (1.50 * trt_assign_absent)
    g <- sapply(trt_assign_vec, plogis)
    dat_sampled$trt <- rbinom(n = nrow(dat_sampled), size = 1, prob = g)
    
  }
  
  # treatment heterogeneity model
  likemath <- dat_sampled$LIKEMATH
  
  dat_sampled$ITE <- (ifelse(likemath >= 0 & likemath < 9, (trt_ceiling / 9) * likemath,
                             ifelse(likemath >= 9 & likemath < 12, trt_ceiling,
                                    ifelse(likemath >= 12 & likemath < 14, (7 * trt_ceiling) - ((trt_ceiling / 2) * likemath), 0)))) + 0.5
  
  # error models
  C <- model.matrix(~factor(dat_sampled$IDSCH)-1)
  u <- rnorm(ncol(C), mean = 0, sd = sqrt(tau_2))
  dat_sampled$e <- rnorm(nrow(dat_sampled), mean = 0, sd = sqrt(sigma_2))
  dat_sampled$U <- as.numeric(C %*% u)
  
  # get covariate model matrix
  X <- dat_sampled[, c(covs_cont, covs_cat_gen)]
  N <- nrow(X)
  Xmat <- cbind(rep(1, N), data.matrix(X))
  
  # outcome model
  dat_sampled$y0 <- as.numeric((Xmat %*% beta_vec) + dat_sampled$U + dat_sampled$e)
  dat_sampled$y1 <- as.numeric(dat_sampled$y0 + dat_sampled$ITE)
  dat_sampled$Y <- ifelse(dat_sampled$trt == 1, dat_sampled$y1, dat_sampled$y0)
  
  # return final data
  testing_vars <- c("y0", "y1", "U", "e")
  
  if(isTRUE(testing)) {
    
    dat_final <- dat_sampled
    
  } else {
    
    dat_final <- dat_sampled[, !names(dat_sampled) %in% c(testing_vars, covs_cat_gen)]
    
  }
  
  return(dat_final)
  
}

# Helper Functions --------------------------------------------------------

# define propensity score estimation function
analysis_propensity <- function(dat) {
  
  fit <- WeightIt::weightit(trt ~ STUSEX + ABSENT + PARENTED + FREPL + PRIV + HOMERES + BELONG + LIKEMATH + CONFMATH + ACADSUC,
                            data = dat,
                            method = "bart",
                            estimand = "ATE")
  
  propensity_scores <- fit$ps
  
  return(propensity_scores)
  
}

# define BCF random intercepts and slopes function
random_intercepts_slopes <- function (groups, treatment, original_colnames = TRUE) {
  
  # get dummy codes
  dummies = multibart::get_dummies(groups, original_colnames = original_colnames)
  
  # format data
  random_des = cbind(dummies, dummies * treatment)
  Q = matrix(0, nrow = ncol(random_des), ncol = 2)
  Q[1:ncol(random_des)/2, 1] = 1
  Q[(1 + ncol(random_des)/2):nrow(Q), 2] = 1
  
  # get random intercepts
  intercept_ix = get_group_indices(dummies)
  
  # get random treatments
  treatment_ix = get_group_indices(dummies * treatment)
  
  # package and send out results
  return(list(randeff_design = random_des, randeff_variance_component_design = as.matrix(Q), 
              group_dummies = dummies, intercept_ix = intercept_ix, 
              treatment_ix = treatment_ix))
}

# define BCF random intercepts function
random_intercepts <- function (groups, original_colnames = TRUE) {
  
  # get dummy codes
  dummies <- multibart::get_dummies(groups, original_colnames = original_colnames)
  
  # format data
  random_des <- cbind(dummies)
  Q <- matrix(1, nrow = ncol(random_des), ncol = 1)
  
  # get random intercepts
  intercept_ix = get_group_indices(dummies)
  
  # package and send out results
  return(list(randeff_design = random_des, randeff_variance_component_design = as.matrix(Q), 
              group_dummies = dummies, intercept_ix = intercept_ix))
}

# define summarization function
summarize_analysis <- function(samples_ite) {
  
  # get average ITE for each individual across all tree draws
  ite_pred <- apply(samples_ite, 2, mean)
  
  # return results
  return(ite_pred)
  
  
}

# get ICC of data
get_ICC <- function(dat, treatment) {
  
  # create treatment column name
  trt <- as.symbol(treatment)
  
  # get MLM fit
  fit <- lme4::lmer(formula = Y ~ trt + STUSEX + ABSENT + PARENTED + FREPL + PRIV + HOMERES + BELONG + LIKEMATH + CONFMATH + ACADSUC + (1 | IDSCH),
                    data = dat)
  
  # get ICC
  icc <- as.numeric(performance::icc(fit)[[1]])
  
  return(icc)
  
}

# Methods Functions -------------------------------------------------------

# define causal forest method
analysis_CF <- function(dat, treatment) {
  
  # change variable type for cluster identifier
  dat$IDSCH <- as.character(dat$IDSCH)
  
  # create data objects
  X <- dat %>%
    dplyr::select(-c(Y, IDSTU, as.symbol(treatment), weights, ITE))
  
  X <- model.matrix(~., data = X)
  
  trt <- dat[, treatment]
  Y <- as.numeric(dat$Y)
  
  fit <- grf::causal_forest(X = X, Y = Y, W = trt)
  
  tauhat <- predict(fit)
  
  ite <-  as.numeric(tauhat$predictions)
  
  return(ite)
  
}

# define stan4bart method
analysis_StanBART <- function(dat, treatment) {
  
  # get symbol version of treatment vector
  trt <- as.symbol(treatment)
  
  # change variable type of grouping variable
  dat$IDSCH <- as.character(dat$IDSCH)
  
  # fit stan4bart model to data
  fit <- stan4bart::stan4bart(formula = 
                                Y ~ bart(. - IDSTU - IDSCH - ITE - weights) + 
                                (1 | IDSCH),
                              data = dat,
                              treatment = trt,
                              chains = 1,
                              cores = 1,
                              iter = 2000,
                              warmup = 1000,
                              skip = 1,
                              stan_args = list(verbose = FALSE),
                              bart_args = list(keepTrees = TRUE, 
                                               n.burn = 1000,
                                               n.thin = 1,
                                               n.trees = 200, 
                                               n.threads = 1, 
                                               n.chains = 1))
  
  # extract observed and counterfactual fits
  samples_mu_train <- dbarts::extract(fit)
  samples_mu_test <- dbarts::extract(fit, sample = "test")
  
  # get samples matrix of ite
  samples_ite <- (samples_mu_train - samples_mu_test) * (2 * fit$frame[[fit$treatment]] - 1)
  samples_ite <- t(samples_ite)
  
  return(samples_ite)
  
}

# define BCF method
analysis_BCF <- function(dat, treatment) {
  
  # change variable type for cluster identifier
  dat$IDSCH <- as.character(dat$IDSCH)
  
  # create data objects
  X_control <- dat %>% 
    dplyr::select(-c(Y, IDSTU, as.symbol(treatment), weights, ITE)) %>%
    multibart::mb_modelmatrix()
  
  X_moderate <- dat %>%
    dplyr::select(LIKEMATH) %>%
    multibart::mb_modelmatrix()
  
  trt <- dat[, treatment]
  Y <- as.numeric(dat$Y)
  pihat <- as.numeric(dat$weights)
  
  # fit model
  fit_BCF <- multibart::bcf_binary(y = Y, z = trt, 
                                   x_control = X_control, 
                                   x_moderate = X_moderate,
                                   pihat = pihat,
                                   base_moderate = .95, 
                                   ntree_moderate = 100,
                                   power_moderate = 2, 
                                   sd_moderate = .5, 
                                   nburn = 1000, nsim = 1000, nthin = 1,
                                   include_pi = "control")
  
  # extract mu and tau from model
  fitted_mu <- fit_BCF$control_fit
  fitted_tau <- fit_BCF$moderate_fit
  
  # get samples matrix of ITE
  samples_ite <- multibart::get_forest_fit(fitted_tau, X_moderate)
  
  # return results
  return(samples_ite)
  
}

# define MBCF method
analysis_MBCF <- function(dat, treatment) {
  
  # change variable type for cluster identifier
  dat$IDSCH <- as.character(dat$IDSCH)
  
  # create data objects
  X_control <- dat %>% 
    dplyr::select(-c(Y, IDSTU, IDSCH, as.symbol(treatment), weights, ITE)) %>%
    multibart::mb_modelmatrix()
  
  X_moderate <- dat %>%
    dplyr::select(LIKEMATH) %>%
    multibart::mb_modelmatrix()
  
  trt <- dat[, treatment]
  Y <- as.numeric(dat$Y)
  pihat <- as.numeric(dat$weights)
  
  # random effects
  random_effect_setup <- random_intercepts(groups = dat$IDSCH, 
                                           original_colnames = TRUE)
  
  # fit model
  fit_mBCF <- multibart::bcf_binary(y = Y, z = trt, 
                                    x_control = X_control, 
                                    x_moderate = X_moderate,
                                    pihat = pihat,
                                    randeff_design = random_effect_setup$randeff_design, 
                                    randeff_variance_component_design = random_effect_setup$randeff_variance_component_design, 
                                    base_moderate = .95, 
                                    ntree_moderate = 100,
                                    power_moderate = 2, 
                                    sd_moderate = .5, 
                                    nburn = 1000, nsim = 1000, nthin = 1,
                                    include_pi = "control")
  
  # extract mu and tau from model
  fitted_mu <- fit_mBCF$control_fit
  fitted_tau <- fit_mBCF$moderate_fit
  
  # get samples matrix of ITE
  samples_ite <- multibart::get_forest_fit(fitted_tau, X_moderate)
  
  # return results
  return(samples_ite)
  
}

# Estimation Function -----------------------------------------------------

# define estimation procedure
estimate <- function(dat, treatment) {
  
  ite_mat_list <- list()
  
  # get propensity scores for input data
  dat$weights <- analysis_propensity(dat = dat)
  gc()
  
  # estimate ITEs for Bayesian methods
  ite_mat_list <- list("ITE_stan4bart" = analysis_StanBART(dat = dat, treatment = treatment),
                       "ITE_BCF" = analysis_BCF(dat = dat, treatment = treatment),
                       "ITE_MBCF" = analysis_MBCF(dat = dat, treatment = treatment))
  
  gc()
  
  # estimate ITEs for frequentist methods
  ITE_CF <- analysis_CF(dat = dat, treatment = treatment)
  
  gc()
  
  # create ITE dataframe
  ite_df <- 
    lapply(ite_mat_list, summarize_analysis) %>%
    as.data.frame(., col.names = names(.))
  
  ite_df$ITE_CF <- ITE_CF
  
  # bind ITE and data input
  estimate_dat <- dplyr::bind_cols(dat, ite_df)
  
  # get ICC estimate
  estimate_dat$ICC_est <- rep(get_ICC(dat = estimate_dat, treatment = treatment), nrow(estimate_dat))
  
  # return results
  return(estimate_dat)
  
}

# Driver ------------------------------------------------------------------

run_sim <- function(beta_0, beta_1, beta_2, beta_3, beta_4, beta_5,
                    beta_6, beta_7, beta_8, beta_9, beta_10, beta_11,
                    beta_12, beta_13, beta_14, beta_15, beta_16, beta_17,
                    beta_18, beta_19, beta_20, beta_21, beta_22,
                    trt_ceiling, RCT, prop_PRIV, ICC, J, testing = FALSE, treatment,
                    iterations, submit = submit, seed = seed) {
  
  # check if treatment input is a character object
  stopifnot(is.character(treatment))
  
  # check if seed is set and set seed if not set already
  if (!is.null(seed)) set.seed(seed)
  
  # stop the function if ICC is not between 0 and 1
  stopifnot(ICC >= 0 && ICC < 1)
  
  # create beta vector for data generation input
  beta_vec <- c(beta_0, beta_1, beta_2, beta_3, beta_4, beta_5,
                beta_6, beta_7, beta_8, beta_9, beta_10, beta_11,
                beta_12, beta_13, beta_14, beta_15, beta_16, beta_17,
                beta_18, beta_19, beta_20, beta_21, beta_22)
  
  # run iterations
  purrr::map(1:iterations, ~ {
    
    # data generation
    dat <- generate_data(beta_vec = beta_vec, trt_ceiling = trt_ceiling, RCT = RCT, prop_PRIV = prop_PRIV, sigma_2 = 1 - ICC, tau_2 = ICC, J = J, testing = testing)
    
    # estimation procedure
    estimate(dat = dat, treatment = treatment)
    
  }) %>% dplyr::bind_rows(.id = "iteration")
  
}


