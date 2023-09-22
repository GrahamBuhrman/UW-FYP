### Simulation Analysis Functions ###

# Packages ----------------------------------------------------------------

require(dplyr)
require(tidyr)
require(lme4)
require(simhelpers)
require(ggplot2)
require(ggpubr)
require(cowplot)

# Summary Statistics ------------------------------------------------------

get_PEHE_stats <- function(sim_results) {

  sim_results %>%
    dplyr::group_by(ICC, J, submit) %>%
    dplyr::summarize(PEHE_CF = sqrt(mean((ITE - ITE_CF)^2)),
                     PEHE_stan4bart = sqrt(mean((ITE - ITE_stan4bart)^2)),
                     PEHE_BCF = sqrt(mean((ITE - ITE_BCF)^2)),
                     PEHE_MBCF = sqrt(mean((ITE - ITE_MBCF)^2))) %>%
    dplyr::ungroup() %>% 
    group_by(ICC, J) %>% 
    summarize(PEHE_CF_EST = mean(PEHE_CF), 
              PEHE_CF_SD = sd(PEHE_CF), 
              PEHE_stan4bart_EST = mean(PEHE_stan4bart), 
              PEHE_stan4bart_SD = sd(PEHE_stan4bart), 
              PEHE_BCF_EST = mean(PEHE_BCF), 
              PEHE_BCF_SD = sd(PEHE_BCF), 
              PEHE_MBCF_EST = mean(PEHE_MBCF), 
              PEHE_MBCF_SD = sd(PEHE_MBCF))
  
}

get_PEHE_agg <- function(sim_results) {
  
  sim_results %>%
    dplyr::group_by(ICC, J, submit) %>%
    dplyr::summarize(PEHE_CF = sqrt(mean((ITE - ITE_CF)^2)),
                     PEHE_stan4bart = sqrt(mean((ITE - ITE_stan4bart)^2)),
                     PEHE_BCF = sqrt(mean((ITE - ITE_BCF)^2)),
                     PEHE_MBCF = sqrt(mean((ITE - ITE_MBCF)^2))) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(cols = PEHE_CF:PEHE_MBCF,
                        names_prefix = "PEHE_",
                        names_to = "method",
                        values_to = "PEHE")
  
}

# Aggregate Results -------------------------------------------------------

get_agg_results <- function(sim_results) {
  
  sim_results %>%
    tidyr::pivot_longer(cols = c(ITE_stan4bart:ITE_CF),
                        names_prefix = "ITE_",
                        names_to = "method",
                        values_to = "sampleITE") %>%
    dplyr::group_by(ICC, J, method, IDSTU) %>%
    dplyr::summarize(estITE = mean(sampleITE),
                     lbITE = quantile(sampleITE, 0.025),
                     ubITE = quantile(sampleITE, 0.975)) %>%
    dplyr::ungroup() %>%
    dplyr::select(IDSTU, ICC, J, method, estITE, lbITE, ubITE)
  
}

get_agg_results_plot <- function(sim_results) {
  
  sim_results %>%
    tidyr::pivot_longer(cols = c(ITE_stan4bart:ITE_CF),
                        names_prefix = "ITE_",
                        names_to = "method",
                        values_to = "sampleITE") %>%
    dplyr::group_by(J, method, IDSTU) %>%
    dplyr::summarize(estITE = mean(sampleITE),
                     lbITE = quantile(sampleITE, 0.025),
                     ubITE = quantile(sampleITE, 0.975)) %>%
    dplyr::ungroup() %>%
    dplyr::select(IDSTU, J, method, estITE, lbITE, ubITE)
  
}


# Plot Data ---------------------------------------------------------------

get_plot_dat <- function(sim_results) {
  
  cov_dat <- readRDS("Simulation Data/TIMSS_covariates.rds")
  agg_results <- get_agg_results_plot(sim_results)
  
  gc()
  
  plot_dat <-
    dplyr::left_join(agg_results, cov_dat, by = "IDSTU") %>%
    dplyr::select(IDSCH, IDSTU, J, method, LIKEMATH, estITE, lbITE, ubITE)
  
  gc()
  
  return(plot_dat)

}


# Visualize Data ----------------------------------------------------------

make_PEHE_boxplot <- function(sim_results) {
  
  PEHE_agg_dat <- get_PEHE_agg(sim_results)
  
  PEHE_agg_dat$method <- factor(PEHE_agg_dat$method, levels = c("CF", "stan4bart", "BCF", "MBCF"))
  PEHE_agg_dat$Clusters <- factor(PEHE_agg_dat$J)
    
  ggplot(PEHE_agg_dat,
         aes(x = method,
             y = PEHE,
             fill = method)) +
    geom_boxplot(alpha = 0.5, size = 0.75) +
    geom_hline(yintercept = 0.2, linetype = "dashed", linewidth = 1) +
    scale_fill_manual(values=c("#0479a8", "#c5050c", "#f7941e","#97b85f")) +
    facet_grid(ICC~Clusters, labeller = labeller(ICC = label_both, Clusters = label_both)) +
    theme_bw() +
    labs(x = "Method",
         y = "RPEHE",
         fill = "Method") +
    theme(text = element_text(size = 25, color = "black", face = "plain"),
          legend.position = "none")
  
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
    dplyr::filter(., J == 200, method == method_spec)
  
  ggplot(data = plot_data_filter,
         aes(x = LIKEMATH,
             fill = method_spec)) +
    ylim(0, 2.25) +
    xlim(5, 14) +
    stat_function(fun = fun1, n = 5001, linewidth = 2, linetype = "solid", color = "black") +
    stat_function(fun = fun2, n = 5001, linewidth = 2, linetype = "solid", color = "black") +
    stat_function(fun = fun3, n = 5001, linewidth = 2, linetype = "solid", color = "black") +
    geom_ribbon(aes(ymin = predict(mgcv::gam(lbITE ~ s(LIKEMATH, bs = "cs"))), 
                    ymax = predict(mgcv::gam(ubITE ~ s(LIKEMATH, bs = "cs")))),
                alpha = 0.25,
                fill = color_spec,
                color = color_spec,
                linewidth = 1.5,
                linetype = "dotdash") + 
    geom_smooth(aes(y = estITE), color = color_spec, linewidth = 1.5, method = "gam", se = FALSE) +
    theme_bw() +
    theme(text = element_text(size = 15, color = "black"),
          legend.position = "top",
          axis.title = element_blank()) +
    labs(fill = "Method")
  
}

arrange_CATE_plots <- function(plot_1, plot_2, plot_3, plot_4) {
  
  figure <- ggarrange(plot_1, plot_2, plot_3, plot_4, 
                      ncol = 2, nrow = 2,
                      label.x = 0,
                      label.y = 0,
                      hjust = -1.5,
                      vjust = -5,
                      align = "hv",
                      legend = "top",
                      font.label = list(size = 15, face = "plain"))
  
  annotate_figure(figure, left = ggpubr::text_grob("CATE", rot = 90, size = 30),
                  bottom = ggpubr::text_grob("Like Math", size = 30))
  
}
  

# get arranged cate plots
sim_plot_dat <- get_plot_dat(sim_results_dat)
test_CF_plot <- make_agg_CATE_plot(plot_data = sim_plot_dat, method_spec = "CF", color_spec = "#0479a8")
test_stan4bart_plot <- make_agg_CATE_plot(plot_data = sim_plot_dat, method_spec = "stan4bart", color_spec = "#c5050c")
test_BCF_plot <- make_agg_CATE_plot(plot_data = sim_plot_dat, method_spec = "BCF", color_spec = "#f7941e")
test_MBCF_plot <- make_agg_CATE_plot(plot_data = sim_plot_dat, method_spec = "MBCF", color_spec = "#97b85f")

png(filename="Images/CATE_plots.png",
    type="cairo",
    units="in", 
    width=12, 
    height=12, 
    pointsize=12, 
    res=300)
arrange_CATE_plots(test_CF_plot, test_stan4bart_plot, test_BCF_plot, test_MBCF_plot)
dev.off()

# get boxplot
png(filename="Images/RPEHE_boxplot.png",
    type="cairo",
    units="in", 
    width=16, 
    height=12, 
    pointsize=12, 
    res=300)
make_PEHE_boxplot(sim_results_dat)
dev.off()

##########################################################################################################

# get arranged CATE plots
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
