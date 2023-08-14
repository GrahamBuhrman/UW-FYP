### Poster Visualizations and Reports - TIMSS 2019 Data Simulation ###  


# Packages ----------------------------------------------------------------

require(ggplot2)
require(dplyr)
require(tidyr)
require(Cairo)

# Source Scripts ----------------------------------------------------------

source("Simulation Scripts/Load Covariates.R")
source("Simulation Scripts/Data Generation.R")
source("Simulation Scripts/Estimation.R")

# Generate Data -----------------------------------------------------------

set.seed(2023719)

poster_data <- generate_data(testing = FALSE, J = 100)
poster_data$weights <- analysis_propensity(poster_data)


# Get Estimations ---------------------------------------------------------

MLM_estimates <- analysis_MLM(dat = poster_data, one_run = TRUE)
gc()

# get X-BART estimates
XBART_estimates <- analysis_XBART(dat = poster_data, one_run = TRUE)
gc()

# get stan4bart estimates
StanBART_estimates <- analysis_StanBART(dat = poster_data, one_run = TRUE)
gc()

# get BCF estimates
BCF_estimates <- analysis_BCF(dat = poster_data, one_run = TRUE)
gc()

# get MBCF estimates
MBCF_estimates <- analysis_MBCF(dat = poster_data, one_run = TRUE)
gc()


# Visualize Data ----------------------------------------------------------

# Treatment Model Plot -----------------------------------------------

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


# define plot functions
tau_function_plot <- function(dat) {
  
  ggplot(data = dat,
         aes(x = LIKEMATH)) +
    stat_function(fun = fun1, n = 5001, linewidth = 1.5) +
    stat_function(fun = fun2, n = 5001, linewidth = 1.5) +
    stat_function(fun = fun3, n = 5001, linewidth = 1.5) +
    geom_vline(xintercept = 9, linetype = "dashed", color = "#c5050c", linewidth = 1.5) +
    geom_vline(xintercept = 12, linetype = "dashed", color = "#0479a8", linewidth = 1.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1.5) +
    ylim(0, 2.5) +
    xlim(0, 15) +
    theme_minimal() +
    theme(axis.text = element_text(size = 30, family = "Red Hat Display", color = "black"),
          axis.title = element_text(size = 36, family = "Red Hat Display"),
          axis.title.x = element_text(vjust = -0.5, family = "Red Hat Display"),
          title = element_text(size = 38, face = "bold", family = "Red Hat Display")) +
    labs(y = "ITE",
         x = "Like Math")
  
  
}

# generate treatment model plot
png(filename="Images/Treatment_Curve.png", 
    type="cairo",
    units="in", 
    width=14, 
    height=6, 
    pointsize=12, 
    res=300)
tau_function_plot(poster_data)
dev.off()


# ITE Plot -----------------------------------------------------------

# get plot data
plot_data <- 
  poster_data %>%
  dplyr::mutate(`X-BART` = XBART_estimates$ite_pred,
                `Stan4bart` = StanBART_estimates$ite_pred,
                `BCF` = BCF_estimates$ite_pred,
                `MBCF` = MBCF_estimates$ite_pred,
                `Oracle` = ITE) %>%
  dplyr::select(-ITE) %>%
  dplyr::select(LIKEMATH, `X-BART`, `Stan4bart`, `BCF`, `MBCF`)

# define plot function
ite_plot <- function(dat) {
  
  ggplot(data = dat,
         aes(x = LIKEMATH)) +
    ylim(0, 2.5) +
    xlim(5, 15) +
    geom_smooth(aes(y = `BCF`, linetype = "BCF", color = "BCF"), method = "gam", se = FALSE) +
    geom_smooth(aes(y = `MBCF`, linetype = "MBCF", color = "MBCF"), method = "gam", se = FALSE) +
    geom_smooth(aes(y = `Stan4bart`, linetype = "Stan4bart", color = "Stan4bart"), method = "gam", se = FALSE) +
    geom_smooth(aes(y = `X-BART`, linetype = "X-BART", color = "X-BART"), method = "gam", se = FALSE) +
    stat_function(fun = fun1, n = 5001, linewidth = 1.5, aes(linetype = "Oracle", color = "Oracle")) +
    stat_function(fun = fun2, n = 5001, linewidth = 1.5, aes(linetype = "Oracle", color = "Oracle")) +
    stat_function(fun = fun3, n = 5001, linewidth = 1.5, aes(linetype = "Oracle", color = "Oracle")) +
    theme_minimal() +
    scale_color_manual(name = "Method", 
                       breaks = c("BCF", "MBCF", "Stan4bart", "X-BART", "Oracle"), 
                       values = c("BCF" = "#f7941e", "MBCF" = "#97b85f", "Stan4bart" = "#c5050c", "X-BART" = "#0479a8", "Oracle" = "black")) +
    scale_linetype_manual(name = "Method", 
                       breaks = c("BCF", "MBCF", "Stan4bart", "X-BART", "Oracle"), 
                       values = c("BCF" = "solid", "MBCF" = "solid", "Stan4bart" = "solid", "X-BART" = "solid", "Oracle" = "dashed")) +
    theme(axis.text = element_text(size = 30, family = "Red Hat Display", color = "black"),
          axis.title = element_text(size = 36, family = "Red Hat Display"),
          axis.title.x = element_text(vjust = -0.5, family = "Red Hat Display"),
          title = element_text(size = 38, face = "bold", family = "Red Hat Display"),
          legend.text = element_text(size = 30, family = "Red Hat Display", color = "black"),
          legend.key.size = unit(1, "in")) +
    labs(y = "Individual Treatment Effect",
         x = "Like Math")
  
}

# generate treatment model plot
png(filename="Images/ITE_plot.png", 
    type="cairo",
    units="in", 
    width=20, 
    height=12, 
    pointsize=12, 
    res=300)
ite_plot(plot_data)
dev.off()


# LIKEMATH Histogram -------------------------------------------------

likemath_hist <- function(dat) {
  
  ggplot(data = dat,
         aes(x = LIKEMATH)) +
    geom_histogram(aes(y = after_stat(density)),
                   fill = "#dadfe1",
                   color = "black",
                   binwidth = 1) +
    xlim(5, 15) +
    geom_vline(xintercept = mean(dat$LIKEMATH), linewidth = 1.5, linetype = "dashed", color = "black") +
    geom_density(col = "#c5050c", linewidth = 2) +
    theme_minimal() +
    theme(axis.text = element_text(size = 30, family = "Red Hat Display", color = "black"),
          axis.title = element_text(size = 36, family = "Red Hat Display"),
          axis.title.x = element_text(vjust = -0.5, family = "Red Hat Display"),
          title = element_text(size = 38, face = "bold", family = "Red Hat Display")) +
    labs(y = "Density",
         x = "Like Math")
  
}

# generate treatment model plot
png(filename="Images/LIKEMATH_hist.png", 
    type="cairo",
    units="in", 
    width=10, 
    height=6, 
    pointsize=12, 
    res=300)
likemath_hist(plot_data)
dev.off()


# RMSE Table --------------------------------------------------------------

rmse_function <- function(x) {
  
  sqrt(mean(x^2))
  
}

outcome <- poster_data$Y

table_data <-
  poster_data %>%
  dplyr::mutate(`ITE_MLM` = MLM_estimates$ite_pred,
                `ITE_X-BART` = XBART_estimates$ite_pred,
                `ITE_Stan4bart` = StanBART_estimates$ite_pred,
                `ITE_BCF` = BCF_estimates$ite_pred,
                `ITE_MBCF` = MBCF_estimates$ite_pred,
                `ITE_TRUE` = ITE) %>%
  dplyr::mutate(`diff_MLM` = `ITE_MLM` - `ITE_TRUE`,
                `diff_X-BART` = `ITE_X-BART` - `ITE_TRUE`,
                `diff_Stan4bart` = `ITE_Stan4bart` - `ITE_TRUE`,
                `diff_BCF` = `ITE_BCF` - `ITE_TRUE`,
                `diff_MBCF` = `ITE_MBCF` - `ITE_TRUE`) %>%
  dplyr::select(-ITE) %>%
  tidyr::pivot_longer(cols = c(`diff_MLM`:`diff_MBCF`),
                      names_to = c("Method"),
                      names_prefix = "diff_",
                      values_to = "diff") %>%
  dplyr::group_by(Method) %>%
  summarize(RMSE = rmse_function(diff)) %>%
  mutate(NRMSE = RMSE / sd(outcome))

