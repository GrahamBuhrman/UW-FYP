### Simulation Analysis Functions ###

# Packages ----------------------------------------------------------------

require(dplyr)
require(tidyr)
require(lme4)
require(simhelpers)
require(ggplot2)

# Summary Statistics ------------------------------------------------------

get_summary_stats <- function(sim_results) {
  
  sim_results %>%
    dplyr::group_by(RCT, ICC, J, submit) %>%
    dplyr::summarize(CATE_TRUE = mean(ITE),
                     CATE_CF = mean(ITE_CF),
                     CATE_stan4bart = mean(ITE_stan4bart),
                     CATE_BCF = mean(ITE_BCF),
                     CATE_MBCF = mean(ITE_MBCF)) %>%
    dplyr::ungroup()
  
}

get_PEHE_stats <- function(sim_results) {

  sim_results %>%
    dplyr::group_by(RCT, J, submit) %>%
    dplyr::summarize(PEHE_CF = sqrt(mean((ITE - ITE_CF)^2)),
                     PEHE_stan4bart = sqrt(mean((ITE - ITE_stan4bart)^2)),
                     PEHE_BCF = sqrt(mean((ITE - ITE_BCF)^2)),
                     PEHE_MBCF = sqrt(mean((ITE - ITE_MBCF)^2))) %>%
    dplyr::ungroup() %>% 
    group_by(RCT, J) %>% 
    summarize(PEHE_CF_EST = mean(PEHE_CF), 
              PEHE_CF_SD = sd(PEHE_CF), 
              PEHE_stan4bart_EST = mean(PEHE_stan4bart), 
              PEHE_stan4bart_SD = sd(PEHE_stan4bart), 
              PEHE_BCF_EST = mean(PEHE_BCF), 
              PEHE_BCF_SD = sd(PEHE_BCF), 
              PEHE_MBCF_EST = mean(PEHE_MBCF), 
              PEHE_MBCF_SD = sd(PEHE_MBCF))
  
}

# Performance Calculations ------------------------------------------------

get_abs_perf <- function(sim_results) {
  
  summary_stats <- get_summary_stats(sim_results)
  
  res_dat <- 
    summary_stats %>%
    tidyr::pivot_longer(cols = c(CATE_CF, CATE_stan4bart, CATE_BCF, CATE_MBCF), 
                        names_prefix = "CATE_", 
                        names_to = "Method", 
                        values_to = "Estimate") %>% 
    dplyr::group_by(Method, RCT, ICC, J) %>% 
    do(simhelpers::calc_absolute(., estimates = Estimate, true_param = CATE_TRUE)) %>%
    dplyr::ungroup()
  
}

get_abs_perf_agg <- function(sim_results) {
  
  summary_stats <- get_summary_stats(sim_results)
  
  res_dat <- 
    summary_stats %>%
    tidyr::pivot_longer(cols = c(CATE_CF, CATE_stan4bart, CATE_BCF, CATE_MBCF), 
                        names_prefix = "CATE_", 
                        names_to = "Method", 
                        values_to = "Estimate") %>% 
    dplyr::group_by(Method, RCT, J) %>% 
    do(simhelpers::calc_absolute(., estimates = Estimate, true_param = CATE_TRUE)) %>%
    dplyr::ungroup()
  
}

# Aggregate Results -------------------------------------------------------

get_agg_results <- function(sim_results) {
  
  sim_results %>%
    tidyr::pivot_longer(cols = c(ITE_stan4bart:ITE_CF),
                        names_prefix = "ITE_",
                        names_to = "method",
                        values_to = "sampleITE") %>%
    dplyr::group_by(RCT, ICC, J, method, IDSTU) %>%
    dplyr::summarize(estITE = mean(sampleITE),
                     lbITE = quantile(sampleITE, 0.025),
                     ubITE = quantile(sampleITE, 0.975)) %>%
    dplyr::ungroup() %>%
    dplyr::select(IDSTU, RCT, ICC, J, method, estITE, lbITE, ubITE)
  
}

get_agg_results_plot <- function(sim_results) {
  
  sim_results %>%
    tidyr::pivot_longer(cols = c(ITE_stan4bart:ITE_CF),
                        names_prefix = "ITE_",
                        names_to = "method",
                        values_to = "sampleITE") %>%
    dplyr::group_by(RCT, J, method, IDSTU) %>%
    dplyr::summarize(estITE = mean(sampleITE),
                     lbITE = quantile(sampleITE, 0.025),
                     ubITE = quantile(sampleITE, 0.975)) %>%
    dplyr::ungroup() %>%
    dplyr::select(IDSTU, RCT, J, method, estITE, lbITE, ubITE)
  
}


# Plot Data ---------------------------------------------------------------

get_plot_dat <- function(sim_results) {
  
  cov_dat <- readRDS("Simulation Data/TIMSS_covariates.rds")
  agg_results <- get_agg_results_plot(sim_results)
  
  gc()
  
  plot_dat <-
    dplyr::left_join(agg_results, cov_dat, by = "IDSTU") %>%
    dplyr::select(IDSCH, IDSTU, RCT, J, method, LIKEMATH, estITE, lbITE, ubITE)
  
  gc()
  
  return(plot_dat)

}


# Visualize Data ----------------------------------------------------------

make_perf_plot <- function(perf_data, grid_spec_1, grid_spec_2, x_axis_spec, perf_spec) {
  
  ggplot(data = perf_data, 
         aes(x = x_axis_spec, 
             y = perf_spec, 
             color = Method, 
             group = Method)) + 
    geom_point() + 
    geom_line() + 
    facet_grid(cols = vars(grid_spec_1), rows = vars(grid_spec_2))
  
}

make_agg_CATE_plot <- function(plot_data, method_spec, color_spec){
  
  # define piecewise functions
  fun1 <- function(x) {
    
    (ifelse(x >= 0 & x < 9, (1.50/9) * x, NA)) + 0.5
    
  }
  
  fun2 <- function(x) {
    
    (ifelse(x >= 9 & x < 12, 1.50, NA)) + 0.5
    
  }
  
  fun3 <- function(x) {
    
    ifelse(x >= 12 & x < 14, (7 * 1.50) - ((1.50 / 2) * x), NA) + 0.5
    
  }
  
  plot_data_filter <- 
    plot_data %>%
    dplyr::filter(., J == 200, RCT == FALSE, method == method_spec)
  
  ggplot(data = plot_data_filter,
         aes(x = LIKEMATH)) +
    ylim(0, 2.25) +
    xlim(5, 14) +
    stat_function(fun = fun1, n = 5001, linewidth = 3, linetype = "solid", color = "black") +
    stat_function(fun = fun2, n = 5001, linewidth = 3, linetype = "solid", color = "black") +
    stat_function(fun = fun3, n = 5001, linewidth = 3, linetype = "solid", color = "black") +
    geom_ribbon(aes(ymin = predict(mgcv::gam(lbITE ~ s(LIKEMATH, bs = "cs"))), 
                    ymax = predict(mgcv::gam(ubITE ~ s(LIKEMATH, bs = "cs")))),
                alpha = 0.25,
                fill = color_spec,
                color = color_spec,
                linewidth = 2,
                linetype = "dotdash") + 
    geom_smooth(aes(y = estITE), color = color_spec, linewidth = 2, method = "gam", se = FALSE) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size = 30, family = "Verdana", color = "black"),
          axis.title = element_blank())
  
}


# get different plots
png(filename="Images/CF_plot_black_3.png", 
    type="cairo",
    units="in", 
    width=16, 
    height=12, 
    pointsize=12, 
    res=300)
make_agg_CATE_plot(plot_data = sim_plot_dat, method_spec = "CF", color_spec = "#0479a8")
dev.off()

png(filename="Images/BCF_plot_black_3.png", 
    type="cairo",
    units="in", 
    width=16, 
    height=12, 
    pointsize=12, 
    res=300)
make_agg_CATE_plot(plot_data = sim_plot_dat, method_spec = "BCF", color_spec = "#f7941e")
dev.off()

png(filename="Images/MBCF_plot_black_3.png", 
    type="cairo",
    units="in", 
    width=16, 
    height=12, 
    pointsize=12, 
    res=300)
make_agg_CATE_plot(plot_data = sim_plot_dat, method_spec = "MBCF", color_spec = "#97b85f")
dev.off()

png(filename="Images/stan4bart_plot_black_3.png", 
    type="cairo",
    units="in", 
    width=16, 
    height=12, 
    pointsize=12, 
    res=300)
make_agg_CATE_plot(plot_data = sim_plot_dat, method_spec = "stan4bart", color_spec = "#c5050c")
dev.off()
