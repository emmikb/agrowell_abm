library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
select <- dplyr::select

dr <- "/home/emily/SL_ABM/param_check/"  #directory to write out plots
setwd('/data/emily/param_check/')

############################################################################################################################
#LOAD DATASETS
############################################################################################################################

load_eu_data <- function(risk_model = "rm0") {
  
  filelist <- Sys.glob(paste(risk_model,"*",sep=""))
  datalist = list()

  for (i in 1:length(filelist)) {
    data_full <- read_csv(filelist[i]) 
    data <- data_full[data_full$iter > 9,]  #drop burn in period of 9 seasons
    datalist[[i]] <- data
  }

  rm0_data = do.call(rbind, datalist)

  #prepare columns
  rm0_df <<- rm0_data %>%
    mutate(simulation = ordered(simulation),
           tank_level = ordered(tank, levels = c(1,2,3),
                                labels = c('Reservoir low', 'Reservoir normal', 'Reservoir high')),
           rainfall = ordered(rf, levels = c(1,2,3),
                              labels = c('Dry', 'Normal', 'Wet')),
           agrowell = ordered(aw, levels = c(0,1), labels = c('False', 'True')))
  
  #summary data for each iteration and simulation
  rm0_iter_df <<- rm0_df %>%
    select(-agrowell, -rainfall) %>%
    group_by(iter, simulation, tank_level, r) %>%  #param here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    summarize(paddy = sum(crop == 0) / n(),  ofc = 1 - paddy,
                bethma = sum(bethma) / n(),
                income = mean(profit), agrowell = sum(aw) / n(),
                rainfall = mean(rf)) %>%
    ungroup()
    
  return(list(rm0_iter_df = rm0_iter_df, rm0_df = rm0_df))
}
load_eu_data()

load_regret_data <- function(risk_model = "rm1") {
  
  filelist <- Sys.glob(paste(risk_model,"*",sep=""))
  datalist = list()
  
  for (i in 1:length(filelist)) {
    data_full <- read_csv(filelist[i]) 
    data <- data_full[data_full$iter > 9,]  #drop burn in period of 9 seasons
    datalist[[i]] <- data
  }
  
  rm1_data = do.call(rbind, datalist)
  
  #prepare columns
  rm1_df <<- rm1_data %>%
    mutate(simulation = ordered(simulation),
           tank_level = ordered(tank, levels = c(1,2,3),
                                labels = c('Reservoir low', 'Reservoir normal', 'Reservoir high')),
           rainfall = ordered(rf, levels = c(1,2,3),
                              labels = c('Dry', 'Normal', 'Wet')),
           agrowell = ordered(aw, levels = c(0,1), labels = c('False', 'True')))
  
  #summary data for each iteration and simulation
  rm1_iter_df <<- rm1_df %>%
    select(-agrowell, -rainfall) %>%
    group_by(iter, simulation, tank_level, b, k, rr) %>%  
    summarize(paddy = sum(crop == 0) / n(),  ofc = 1 - paddy,
              bethma = sum(bethma) / n(),
              income = mean(profit), agrowell = sum(aw) / n(),
              rainfall = mean(rf)) %>%
    ungroup()
  
  return(list(rm1_iter_df = rm1_iter_df, rm1_df = rm1_df))
}
load_regret_data()

load_prospect_data <- function(risk_model = "rm2") {
  
  filelist <- Sys.glob(paste(risk_model,"*",sep=""))
  datalist = list()
  
  for (i in 1:length(filelist)) {
    data_full <- read_csv(filelist[i]) 
    data <- data_full[data_full$iter > 9,]  #drop burn in period of 9 seasons
    datalist[[i]] <- data
  }
  
  rm2_data = do.call(rbind, datalist)
  
  #prepare columns
  rm2_df <<- rm2_data %>%
    mutate(simulation = ordered(simulation),
           tank_level = ordered(tank, levels = c(1,2,3),
                                labels = c('Reservoir low', 'Reservoir normal', 'Reservoir high')),
           rainfall = ordered(rf, levels = c(1,2,3),
                              labels = c('Dry', 'Normal', 'Wet')),
           agrowell = ordered(aw, levels = c(0,1), labels = c('False', 'True')))
  
  #summary data for each iteration and simulation
  rm2_iter_df <<- rm2_df %>%
    select(-agrowell, -rainfall) %>%
    group_by(iter, simulation, tank_level, lam, alph) %>%  
    summarize(paddy = sum(crop == 0) / n(),  ofc = 1 - paddy,
              bethma = sum(bethma) / n(),
              income = mean(profit), agrowell = sum(aw) / n(),
              rainfall = mean(rf)) %>%
    ungroup()
  
  return(list(rm2_iter_df = rm2_iter_df, rm2_df = rm2_df))
}
load_prospect_data()

############################################################################################################################
#PLOT FUNCTIONS
############################################################################################################################

############################################################################################################################
#1. EXPECTED UTILITY
############################################################################################################################
#ofc plot below zero

eu_plots = function(dataset = rm0_iter_df, dr = "", save_plots = FALSE, show_plots = FALSE) {
  
  eu_income_plot <- ggplot(dataset, aes(x = r, y = income, group=r)) +  #color, agrowell
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0, alpha = 0.05, size = 0.3) +
    geom_smooth() +
    scale_color_brewer(palette = 'Dark2') +
    labs(y = "Income (Rupees)",
         x = "r") +
    facet_grid(~tank_level, labeller = 'label_value') +
    theme_bw(base_size = 10) +
    theme(panel.margin.x = unit(0.75, 'line'),
          strip.background = element_rect(fill = NA))
  
  eu_bethma_plot <- ggplot(dataset, aes(x = r, y = bethma, group=r)) + 
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0, alpha = 0.05, size = 0.3) +
    geom_smooth() +
    ylim(0,1) +
    scale_color_brewer(palette = 'Dark2') +
    labs(y = "Fraction of field canals choosing bethma",
         x = "r") +
    facet_grid(~tank_level, labeller = 'label_value') +
    theme_bw(base_size = 10) +
    theme(panel.margin.x = unit(0.75, 'line'),
          strip.background = element_rect(fill = NA))

  eu_aw_plot <- ggplot(dataset, aes(x = r, y = agrowell, group=r)) +  
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0, alpha = 0.05, size = 0.3) +
    geom_smooth() +
    ylim(0,1) +
    scale_color_brewer(palette = 'Dark2') +
    labs(y = "Fraction of farmers with agrowell",
         x = "r") +
    facet_grid(~tank_level, labeller = 'label_value') +
    theme_bw(base_size = 10) +
    theme(panel.margin.x = unit(0.75, 'line'),
          strip.background = element_rect(fill = NA))
  
  eu_ofc_plot <- ggplot(dataset, aes(x = r, y = ofc, group=r)) +  
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0, alpha = 0.05, size = 0.3) +
    geom_smooth() +
    ylim(0,1)+
    scale_color_brewer(palette = 'Dark2') +
    labs(y = "Fraction of farmers cultivating OFCs",
         x = "r") +
    facet_grid(~tank_level, labeller = 'label_value') +
    theme_bw(base_size = 10) +
    theme(panel.margin.x = unit(0.75, 'line'),
          strip.background = element_rect(fill = NA))

  eu_all_plots <- grid.arrange(eu_aw_plot, eu_bethma_plot, eu_income_plot, eu_ofc_plot)
  
  
  if (show_plots) {
    x11()
    plot(eu_all_plots)
    
  }

  if (save_plots) {
    # ggsave(paste(dr,"eu_income_plot.pdf",sep=""), plot = eu_income_plot, device = 'pdf', width = 6, height = 4)
    # ggsave(paste(dr,"eu_income_plot.png",sep=""), plot = eu_income_plot, device = 'png', width = 6, height = 4, dpi = 600)
    # ggsave(paste(dr,"eu_bethma_plot.pdf",sep=""), plot = eu_bethma_plot, device = 'pdf', width = 6, height = 4)
    # ggsave(paste(dr,"eu_bethma_plot.png",sep=""), plot = eu_bethma_plot, device = 'png', width = 6, height = 4, dpi = 600)
    # ggsave(paste(dr,"eu_aw_plot.pdf",sep=""), plot = eu_aw_plot, device = 'pdf', width = 6, height = 4)
    # ggsave(paste(dr,"eu_aw_plot.png",sep=""), plot = eu_aw_plot, device = 'png', width = 6, height = 4, dpi = 600)
    # ggsave(paste(dr,"eu_ofc_plot.pdf",sep=""), plot = eu_ofc_plot, device = 'pdf', width = 6, height = 4)
    # ggsave(paste(dr,"eu_ofc_plot.png",sep=""), plot = eu_ofc_plot, device = 'png', width = 6, height = 4, dpi = 600)
    ggsave(paste(dr,"eu_all_plots.pdf",sep=""), plot = eu_all_plots, device = 'pdf', width = 16, height = 8)
    ggsave(paste(dr,"eu_all_plots.png",sep=""), plot = eu_all_plots, device = 'png', width = 16, height = 8, dpi = 600)
  }

  return(list(eu_bethma_plot = eu_bethma_plot, eu_ofc_plot = eu_ofc_plot, eu_aw_plot = eu_aw_plot, eu_income_plot = eu_income_plot))
}
eu_plots_out <- eu_plots(dataset = rm0_iter_df, dr=dr, show_plots = TRUE, save_plots=TRUE)
#eu_plots_out["eu_bethma_plot"]



############################################################################################################################
#2. REGRET-ADJUSTED UTILITY
############################################################################################################################
#fix axes on final large plot

high_tank <- rm1_iter_df[rm1_iter_df$tank_level == "Reservoir high",]
normal_tank <- rm1_iter_df[rm1_iter_df$tank_level == "Reservoir normal",]
low_tank <- rm1_iter_df[rm1_iter_df$tank_level == "Reservoir low",]

regret_plot <- function(dataset = NULL, title = "", y = NULL, ylab = "", dr = "", plot_fname = "", save_plots = FALSE, show_plots = FALSE) {
  regret_plot <-  regret_income_plot_high <- ggplot(dataset, aes(x = rr, y = y, group=rr)) + 
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0, alpha = 0.05, size = 0.3) +
    ylim(0,1) +
    geom_point(color = "grey50") +
    geom_smooth() +
    scale_color_brewer(palette = 'Dark2') +
    labs(y = ylab,
         x = "rr") +
    facet_grid(b~k, labeller = label_both) +
    theme_bw(base_size = 10) +
    theme(panel.margin.x = unit(0.75, 'line'),
          strip.background = element_rect(fill = NA)) +
    labs(title=title)
  
  if (show_plots) {
    x11()
    plot(regret_plot)  
  }
  
  if (save_plots) {
    ggsave(paste(dr,plot_fname,".pdf",sep=""), plot = regret_plot, device = 'pdf', width = 6, height = 4)
    ggsave(paste(dr,plot_fname,".png",sep=""), plot = regret_plot, device = 'png', width = 6, height = 4, dpi = 600)
  }
  
  return(regret_plot)
}
#regret_plot(dataset = high_tank, title = "Reservoir high", y = high_tank$income, ylab = "Income (Rupees)", plot_fname = "test_plot", save_plots = TRUE, show_plots = TRUE)

all_regret_plots <- function(save_plots = FALSE, show_plots = FALSE) {
  
  #AGROWELL
  high_regret_agrowell_plot <- regret_plot(dataset = high_tank, title = "Reservoir high", y = high_tank$agrowell, 
                                           ylab = "Fraction of farmers with agrowell", 
                                           plot_fname = "high_regret_agrowell_plot", 
                                           save_plots = FALSE, show_plots = FALSE)
  normal_regret_agrowell_plot <- regret_plot(dataset = normal_tank, title = "Reservoir normal", y = normal_tank$agrowell, 
                                           ylab = "Fraction of farmers with agrowell", 
                                           plot_fname = "normal_regret_agrowell_plot", 
                                           save_plots = FALSE, show_plots = FALSE)
  low_regret_agrowell_plot <- regret_plot(dataset = low_tank, title = "Reservoir low", y = low_tank$agrowell, 
                                           ylab = "Fraction of farmers with agrowell", 
                                           plot_fname = "low_regret_agrowell_plot", 
                                           save_plots = FALSE, show_plots = FALSE)
  regret_agrowell <- grid.arrange(high_regret_agrowell_plot, normal_regret_agrowell_plot, low_regret_agrowell_plot)
                                  
  
  #BETHMA
  high_regret_bethma_plot <- regret_plot(dataset = high_tank, title = "Reservoir high", y = high_tank$bethma, 
                                           ylab = "Fraction of farmers with bethma", 
                                           plot_fname = "high_regret_bethma_plot", 
                                           save_plots = FALSE, show_plots = FALSE)
  normal_regret_bethma_plot <- regret_plot(dataset = normal_tank, title = "Reservoir normal", y = normal_tank$bethma, 
                                             ylab = "Fraction of farmers with bethma", 
                                             plot_fname = "normal_regret_bethma_plot", 
                                             save_plots = FALSE, show_plots = FALSE)
  low_regret_bethma_plot <- regret_plot(dataset = low_tank, title = "Reservoir low", y = low_tank$bethma, 
                                          ylab = "Fraction of farmers with bethma", 
                                          plot_fname = "low_regret_bethma_plot", 
                                          save_plots = FALSE, show_plots = FALSE)
  regret_bethma <- grid.arrange(high_regret_bethma_plot, normal_regret_bethma_plot, low_regret_bethma_plot)
  
  #income
  high_regret_income_plot <- regret_plot(dataset = high_tank, title = "Reservoir high", y = high_tank$income, 
                                         ylab = "Income (Rupees)", 
                                         plot_fname = "high_regret_income_plot", 
                                         save_plots = FALSE, show_plots = FALSE)
  normal_regret_income_plot <- regret_plot(dataset = normal_tank, title = "Reservoir normal", y = normal_tank$income, 
                                           ylab = "Income (Rupees)", 
                                           plot_fname = "normal_regret_income_plot", 
                                           save_plots = FALSE, show_plots = FALSE)
  low_regret_income_plot <- regret_plot(dataset = low_tank, title = "Reservoir low", y = low_tank$income, 
                                        ylab = "Income (Rupees)", 
                                        plot_fname = "low_regret_income_plot", 
                                        save_plots = FALSE, show_plots = FALSE)
  regret_income <- grid.arrange(high_regret_income_plot, normal_regret_income_plot, low_regret_income_plot)
  
  #ofc
  high_regret_ofc_plot <- regret_plot(dataset = high_tank, title = "Reservoir high", y = high_tank$ofc, 
                                         ylab = "Fraction of farmers cultivating OFC", 
                                         plot_fname = "high_regret_ofc_plot", 
                                         save_plots = FALSE, show_plots = FALSE)
  normal_regret_ofc_plot <- regret_plot(dataset = normal_tank, title = "Reservoir normal", y = normal_tank$ofc, 
                                           ylab = "Fraction of farmers cultivating OFC", 
                                           plot_fname = "normal_regret_ofc_plot", 
                                           save_plots = FALSE, show_plots = FALSE)
  low_regret_ofc_plot <- regret_plot(dataset = low_tank, title = "Reservoir low", y = low_tank$ofc, 
                                        ylab = "Fraction of farmers cultivating OFC", 
                                        plot_fname = "low_regret_ofc_plot", 
                                        save_plots = FALSE, show_plots = FALSE)
  regret_ofc <- grid.arrange(high_regret_ofc_plot, normal_regret_ofc_plot, low_regret_ofc_plot)
  
  
  if (show_plots) {
    x11()
    plot(regret_agrowell)
    x11()
    plot(regret_bethma)  
    x11()
    plot(regret_income) 
    x11()
    plot(regret_ofc) 
  }
  
  if (save_plots) {
    ggsave(paste(dr,"regret_agrowell.pdf",sep=""), plot = regret_agrowell, device = 'pdf', width = 16, height = 8)
    ggsave(paste(dr,"regret_agrowell.png",sep=""), plot = regret_agrowell, device = 'png', width = 16, height = 8, dpi = 600)
    ggsave(paste(dr,"regret_bethma.pdf",sep=""), plot = regret_bethma, device = 'pdf', width = 16, height = 8)
    ggsave(paste(dr,"regret_bethma.png",sep=""), plot = regret_bethma, device = 'png', width = 16, height = 8, dpi = 600)
    ggsave(paste(dr,"regret_income.pdf",sep=""), plot = regret_income, device = 'pdf', width = 16, height = 8)
    ggsave(paste(dr,"regret_income.png",sep=""), plot = regret_income, device = 'png', width = 16, height = 8, dpi = 600)
    ggsave(paste(dr,"regret_ofc.pdf",sep=""), plot = regret_ofc, device = 'pdf', width = 16, height = 8)
    ggsave(paste(dr,"regret_ofc.png",sep=""), plot = regret_ofc, device = 'png', width = 16, height = 8, dpi = 600)
  }
  
  return(list(regret_agrowell = regret_agrowell, regret_bethma = regret_bethma, regret_income = regret_income, regret_ofc = regret_ofc))
}
all_regret_plots(save_plots=TRUE, show_plots=TRUE)
#all_regret_plots["regret_agrowell"]



############################################################################################################################
#3. PROSPECT THEORY
############################################################################################################################
#add lambda in facet grid label

prospect_plots = function(dataset = rm2_iter_df, dr = "", save_plots = FALSE, show_plots = FALSE) {
  
  lambda <- factor(dataset$lam, levels = c(1, 2.25, 3.5), labels = c("Lambda: 1", "Lambda: 2.25", "Lambda: 3.5"))
  dataset$lambda <- lambda
  
  prospect_income_plot <- ggplot(dataset, aes(x = alph, y = income, group=alph)) +  #color, agrowell
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0, alpha = 0.05, size = 0.3) +
    geom_smooth() +
    scale_color_brewer(palette = 'Dark2') +
    labs(y = "Income (Rupees)",
         x = "Alpha") +
    facet_grid(lambda~tank_level, labeller = ('label_value')) +
    theme_bw(base_size = 10) +
    theme(panel.margin.x = unit(0.75, 'line'),
          strip.background = element_rect(fill = NA)) 

  prospect_bethma_plot <- ggplot(dataset, aes(x = alph, y = bethma, group=alph)) + 
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0, alpha = 0.05, size = 0.3) +
    geom_smooth() +
    ylim(0,1) +
    scale_color_brewer(palette = 'Dark2') +
    labs(y = "Fraction of field canals choosing bethma",
         x = "Alpha") +
    facet_grid(lambda~tank_level, labeller = 'label_value') +
    theme_bw(base_size = 10) +
    theme(panel.margin.x = unit(0.75, 'line'),
          strip.background = element_rect(fill = NA))
  
  prospect_aw_plot <- ggplot(dataset, aes(x = alph, y = agrowell, group=alph)) +  
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0, alpha = 0.05, size = 0.3) +
    geom_smooth() +
    ylim(0,1) +
    scale_color_brewer(palette = 'Dark2') +
    labs(y = "Fraction of farmers with agrowell",
         x = "Alpha") +
    facet_grid(lambda~tank_level, labeller = 'label_value') +
    theme_bw(base_size = 10) +
    theme(panel.margin.x = unit(0.75, 'line'),
          strip.background = element_rect(fill = NA))
  
  prospect_ofc_plot <- ggplot(dataset, aes(x = alph, y = ofc, group=alph)) +  
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0, alpha = 0.05, size = 0.3) +
    geom_smooth() +
    ylim(0,1) +
    scale_color_brewer(palette = 'Dark2') +
    labs(y = "Fraction of farmers cultivating OFCs",
         x = "Alpha") +
    facet_grid(lambda~tank_level, labeller = 'label_value') +
    theme_bw(base_size = 10) +
    theme(panel.margin.x = unit(0.75, 'line'),
          strip.background = element_rect(fill = NA))
  
  prospect_all_plots <- grid.arrange(prospect_aw_plot, prospect_bethma_plot, prospect_income_plot, prospect_ofc_plot)

  if (show_plots) {
    x11()
    plot(prospect_all_plots)
    
  }
  
  if (save_plots) {
    # ggsave(paste(dr,"prospect_income_plot.pdf",sep=""), plot = prospect_income_plot, device = 'pdf', width = 6, height = 4)
    # ggsave(paste(dr,"prospect_income_plot.png",sep=""), plot = prospect_income_plot, device = 'png', width = 6, height = 4, dpi = 600)
    # ggsave(paste(dr,"prospect_bethma_plot.pdf",sep=""), plot = prospect_bethma_plot, device = 'pdf', width = 6, height = 4)
    # ggsave(paste(dr,"prospect_bethma_plot.png",sep=""), plot = prospect_bethma_plot, device = 'png', width = 6, height = 4, dpi = 600)
    # ggsave(paste(dr,"prospect_aw_plot.pdf",sep=""), plot = prospect_aw_plot, device = 'pdf', width = 6, height = 4)
    # ggsave(paste(dr,"prospect_aw_plot.png",sep=""), plot = prospect_aw_plot, device = 'png', width = 6, height = 4, dpi = 600)
    # ggsave(paste(dr,"prospect_ofc_plot.pdf",sep=""), plot = prospect_ofc_plot, device = 'pdf', width = 6, height = 4)
    # ggsave(paste(dr,"prospect_ofc_plot.png",sep=""), plot = prospect_ofc_plot, device = 'png', width = 6, height = 4, dpi = 600)
    ggsave(paste(dr,"prospect_all_plots.pdf",sep=""), plot = prospect_all_plots, device = 'pdf', width = 16, height = 8)
    ggsave(paste(dr,"prospect_all_plots.png",sep=""), plot = prospect_all_plots, device = 'png', width = 16, height = 8, dpi = 600)
  }
  
  return(list(prospect_all_plots = prospect_all_plots, prospect_bethma_plot = prospect_bethma_plot, prospect_ofc_plot = prospect_ofc_plot, prospect_aw_plot = prospect_aw_plot, prospect_income_plot = prospect_income_plot))
}
prospect_plots_out <- prospect_plots(dataset = rm2_iter_df, dr=dr, show_plots = TRUE, save_plots=TRUE)
#prospect_plots_out["prospect_bethma_plot"]






















#UNDER CONSTRUCTION, regret plot for loop
# all_regret_plots <- function(dataset_list = NULL, title = NULL, ylab = NULL, plot_fname = NULL, save_plots = FALSE, show_plots = FALSE) {
#   
#   for (i in length(dataset_list)) {
#     
#     ds = as.data.frame(dataset_list[i])
#     t = title[i]
#     
#     for (d in length(ylab)) {
#       
#       fname = paste(dr,t,plot_fname[d],sep="")
#       
#       #must be a better way to do this
#       if (d == 1) { y <- ds$agrowell }
#       if (d == 2) { y <- ds$bethma }
#       if (d == 3) { y <- ds$income }
#       if (d == 4) { y <- ds$ofc}
#       
#       #regret_plot(dataset = ds, title = t, y = y, ylab = ylab[d], plot_fname = fname)
#       if ((save_plots) & (show_plots)) {
#         regret_plot(dataset = ds, title = t, y = y, ylab = ylab[d], plot_fname = fname, save_plots = TRUE, show_plots = TRUE)
#       }
#       if ((show_plots) & (!save_plots)) {
#         regret_plot(dataset = ds, title = t, y = y, ylab = ylab[d], plot_fname = fname, save_plots = FALSE, show_plots = TRUE)
#       }
#       if ((!show_plots) & (save_plots)) {
#         regret_plot(dataset = ds, title = t, y = y, ylab = ylab[d], plot_fname = fname, save_plots = TRUE, show_plots = FALSE)
#       }
#       if ((!show_plots) & (!save_plots)) {
#         regret_plot(dataset = ds, title = t, y = y, ylab = ylab[d], plot_fname = fname, save_plots = FALSE, show_plots = FALSE)
#       } }}}
# 
# all_regret_plots(dataset_list = list(high_tank, normal_tank, low_tank),
#                  title = list("High", "Normal", "Low"),
#                  ylab = list("Fraction of farmers with agrowell",
#                              "Fraction of field canals choosing bethma",
#                              "Income (Rupees)",
#                              "Fraction of farmers cultivating OFCs"),
#                  plot_fname = list("_regret_agrowell_plot",
#                                    "_regret_bethma_plot",
#                                    "_regret_income_plot",
#                                    "_regret_ofc_plot"),  
#                  save_plots = FALSE, 
#                  show_plots = TRUE)
