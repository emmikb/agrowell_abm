library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(feather)

select <- dplyr::select

load_data <- function(file = 'jg_iter_data.csv') {
  # abm_df <<- read_csv(file)[,-1] %>%
  abm_df <<- read_feather(file) %>%
    mutate(simulation = ordered(simulation),
           tank_level = ordered(tank, levels = c(1,2,3),
                                labels = c('Reservoir low', 'Reservoir normal', 'Reservoir high')),
           rainfall = ordered(rf, levels = c(1,2,3),
                              labels = c('Dry', 'Normal', 'Wet')),
           agrowell = ordered(aw, levels = c(0,1), labels = c('False', 'True')))

  abm_sdf <<- abm_df %>%
    select(-agrowell, -rainfall) %>%
    group_by(iter, simulation, tank_level) %>%
    summarize(paddy = sum(crop == 0) / n(),  ofc = 1 - paddy,
              bethma = sum(bethma) / n(),
              income = mean(profit), agrowell = sum(aw) / n(),
              rainfall = mean(rf)) %>%
    ungroup()

  # abm_gdf <<- abm_sdf %>% gather(key = key, value = value, -iter, -simulation)

  abm_df_rank <<- abm_df %>% group_by(iter, aw, simulation) %>%
    summarize_each(funs(mean(.)), paddy_b_rank, paddy_nb_rank, ofc_b_rank,
                   ofc_nb_rank, tank) %>%
    ungroup() %>%
    gather(key = crop, value = rank, -iter, -simulation, -aw, -tank)

  abm_df_x <- abm_df %>% select(iter, simulation, farmer, agrowell, tank, ends_with('_rank')) %>%
    gather(key = crop, value = rank, ends_with('_rank')) %>%
    group_by(iter, simulation, crop, agrowell) %>%
    summarize(tank_level = mean(tank),
              first = sum(rank == 4) / n(), second = sum(rank == 3) / n(),
              third = sum(rank == 2) / n(), fourth = sum(rank == 1) / n()) %>%
    ungroup() %>%
    gather(key = rank, value = frac, first, second, third, fourth)
}

if (FALSE) {
  ggplot(abm_df, aes(x = iter, y = profit)) +
    geom_smooth(color = "dark green", size = 1, span = 0.001, n = 200)

  ggplot(abm_df, aes(x = iter, y = profit, color = simulation)) +
    geom_point(size = 1)


  ggplot(abm_sdf, aes(x = iter, y = rainfall, color = simulation)) +
    geom_line(size = 1, alpha = 0.2) + ylim(1,3)

  ggplot(abm_sdf, aes(x = iter, y = tank_level, color = simulation)) +
    geom_line(size = 1, alpha = 0.2) + ylim(1,3)

  ggplot(abm_sdf, aes(x = iter, y = bethma, color = simulation)) +
    geom_line(alpha = 0.2, size = 1) + ylim(0,1)

  ggplot(abm_sdf, aes(x = iter, y = paddy, group = simulation)) +
    geom_line(color = "dark green", size = 1, alpha = 0.3) +
    geom_line(aes(y = ofc), color = "dark blue", size = 1, alpha = 0.3)

  ggplot(filter(abm_gdf, key %in% c('bethma','tank_level', 'paddy', 'ofc', 'agrowell')),
         aes(x = iter, y = value)) +
    geom_step(size = 1, color = "dark blue") +
    facet_grid(key~simulation, scales = 'free_y')

  ggplot(abm_df_rank, aes(x = aw, y = rank, color = crop)) +
    geom_jitter(width = 0.25, height = 0.05, alpha = 0.3) +
    scale_color_brewer(palette = 'Dark2') +
    facet_wrap(~tank)

  ggplot(abm_sdf, aes(x = agrowell, y = bethma)) + geom_jitter(width = 0, height = 0.05) + facet_wrap(~tank_level)
}

plot_agrowell_vs_bethma = function(risk_model = c('eu', 'regret', 'prospect', 'mixed')) {
  risk_model <- match.arg(risk_model)
  if (risk_model == 'eu') {
    load_data('eu_iter_data.csv')
    title <- "Risk-averse expected utility"
    figure_name <- "eu_bethma_vs_agrowell.pdf"
  }
  if (risk_model == 'regret') {
    load_data('regret_iter_data.csv')
    title <- "Regret-averse expected utility"
    figure_name <- "regret_bethma_vs_agrowell.pdf"
  }
  if (risk_model == 'prospect') {
    load_data('prospect_iter_data.csv')
    title <- "Prospect theory loss-aversion"
    figure_name <- "prospect_bethma_vs_agrowell.pdf"
  }
  if (risk_model == 'mixed') {
    load_data('mixed_iter_data.csv')
    title <- "Mixed risk models"
    figure_name <- "mixed_bethma_vs_agrowell.pdf"
  }

  print(names(abm_sdf))

  p <- ggplot(abm_sdf, aes(x = agrowell, y = bethma)) +
    geom_jitter(width = 0.02, height = 0.02, alpha = 0.3) +
    facet_wrap(~tank_level, labeller = 'label_value', ncol = 1) +
    labs(x = "Fraction of  farmers with agrowell",
         y = "Fraction of field canals choosing bethma",
         title = title) +
    theme_bw()

  pdf(figure_name, width = 6, height = 6)
  plot(p)
  dev.off()
  invisible(p)
}

calc_grand_df = function() {
  load_data('eu_iter_data_p.feather')
  grand_sdf <- mutate(abm_sdf, risk_model = "Risk-averse")
  grand_df <- mutate(abm_df, risk_model = "Risk-averse")
  load_data('regret_iter_data_p.feather')
  grand_sdf <- abm_sdf %>%
    mutate(risk_model = "Regret-averse")  %>%
    bind_rows(grand_sdf, .)
  grand_df <- abm_df %>%
    mutate(risk_model = "Regret-averse")  %>%
    bind_rows(grand_df, .)
  load_data('prospect_iter_data_p.feather')
  grand_sdf <- abm_sdf %>%
    mutate(risk_model = "Prospect theory")  %>%
    bind_rows(grand_sdf, .)
  grand_df <- abm_df %>%
    mutate(risk_model = "Prospect theory")  %>%
    bind_rows(grand_df, .)
  load_data('mixed_iter_data_p.feather')
  grand_sdf <- abm_sdf %>%
    mutate(risk_model = "Mixed")  %>%
    bind_rows(grand_sdf, .)
  grand_df <- abm_df %>%
    mutate(risk_model = "Mixed")  %>%
    bind_rows(grand_df, .)

  grand_sdf <- grand_sdf %>% mutate(risk_model = ordered(risk_model,
                                                         levels = c('Risk-averse', 'Regret-averse', 'Prospect theory', 'Mixed')),
                                    tank_level = ordered(tank_level,
                                                         levels = c("Reservoir high", "Reservoir normal", "Reservoir low")))
  grand_df <- grand_df %>% mutate(risk_model = ordered(risk_model,
                                                       levels = c('Risk-averse', 'Regret-averse', 'Prospect theory', 'Mixed')),
                                  tank_level = ordered(tank_level,
                                                       levels = c("Reservoir high", "Reservoir normal", "Reservoir low")))
  # levels(grand_sdf$tank_level) <- rev(levels(grand_sdf$tank_level))
  invisible(list(grand_df = grand_df, grand_sdf = grand_sdf))
}

grand_plot = function(save_plot = FALSE, grand_sdf = NULL, grand_df = NULL) {
  if (is.null(grand_sdf) || is.null(grand_df)) {
    gdf <- calc_grand_df()

    grand_sdf <- gdf$grand_sdf
    grand_sdf <<- grand_sdf
    print(names(grand_sdf))
    grand_df <- gdf$grand_df
    grand_sdf <- grand_sdf
    print(names(grand_df))
  }

  p1 <- ggplot(grand_sdf, aes(x = agrowell, y = bethma)) +
    geom_jitter(width = 0.02, height = 0.02, alpha = 0.1, size = 0.5) +
    geom_smooth(color = alpha("blue",0.5), fill = alpha("blue", 0.1)) +
    scale_x_continuous(limits = c(0,1), expand = c(0,0), breaks = seq(0,1,0.25),
                       labels = c('0', '0.25', '0.50', '0.75', '1')) +
    scale_y_continuous(breaks = seq(0,1.25,0.25)) +
    labs(x = "Fraction of  farmers with agrowell",
         y = "Fraction of field canals choosing bethma") +
    facet_grid(tank_level~risk_model, labeller = 'label_value') +
    theme_bw(base_size = 10) +
    theme(panel.margin.x = unit(0.75, 'line'),
          strip.background = element_rect(fill = NA))
  if (save_plot) {
    ggsave('article/grand_abm_bethma.pdf', plot = p1, device = 'pdf', width = 6, height = 4)
    ggsave('article/grand_abm_bethma.png', plot = p1, device = 'png', width = 6, height = 4, dpi = 600)
  }

  p2 <- ggplot(grand_df, aes(x = iter, y = profit, color = agrowell)) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.05, size = 0.3) +
    geom_smooth(aes(linetype = agrowell, size = agrowell), fill = NA) +
    scale_color_brewer(palette = 'Dark2',
                       name = "Agrowell", breaks = c('True','False')) +
    scale_linetype_manual(values = c(True = 'dashed', False = 'solid'),
                          name = "Agrowell", breaks = c('True','False')) +
    scale_size_manual(values = c(True = 0.5, False = 1),
                      name = "Agrowell", breaks = c('True','False')) +
    labs(y = "Income (Rupees)",
         x = "Growing Season") +
    facet_grid(tank_level~risk_model, labeller = 'label_value') +
    theme_bw(base_size = 10) +
    theme(panel.margin.x = unit(0.75, 'line'),
          legend.background = element_rect(fill = "white"),
          strip.background = element_rect(fill = NA))

  if (save_plot) {
    ggsave('article/grand_abm_profit.pdf', plot = p2, device = 'pdf', width = 6, height = 4)
    ggsave('article/grand_abm_profit.png', plot = p2, device = 'png', width = 6, height = 4, dpi = 600)
  }
  return(list(bethma = p1, profit = p2))
}
