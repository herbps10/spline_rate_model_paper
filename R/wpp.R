library(tidyverse)
library(fpemplus)

source("R/extract_tfr_data.R")

dat <- tfr_dataset() %>%
  pivot_longer(`1953`:`2018`, names_to = "period", values_to = "tfr") %>%
  mutate(period = as.numeric(period)) %>%
  filter(!is.na(tfr)) %>%
  group_by(country_name) %>%
  arrange(period) %>%
  mutate(tfr_lead = lead(tfr),
         tfr_change = tfr_lead - tfr) %>%
  ungroup() %>%
  arrange(country_name, period)

countries_phase2 <- dat %>%
  distinct(country_name, start_phase3) %>%
  filter(start_phase3 == 14)

dat <- dat %>% filter(country_name %in% countries)

dat %>%
  ggplot(aes(x = period, y = tfr)) +
  geom_point() +
  facet_wrap(~country_name)

fit <- tfrplus(
  dat,
  start_year = 1953,
  end_year = 2100,
  y = "tfr",
  year = "period",
  area = "country_name",
  num_knots = 7,
  spline_degree = 2,
  hierarchical_splines = c("intercept", "country_name"),
  max_treedepth = 12,
  adapt_delta = 0.99,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500
)

country <- "Burkina Faso"
plot_indicator(fit, country)
plot_indicator(fit)

plot_transition(fit)
plot_transition(fit, "country_name", country, n = 50) +
  geom_point(
    data = dat %>% filter(country_name == country) %>% group_by(country_name) %>% mutate(first_tfr = first(tfr)) %>% rename(name = country_name),
    aes(x = (tfr - 1) / (first_tfr - 1), y = tfr_change)
  )

plot_transition(fit, "country_name", "Afghanistan", n = 50) +
  scale_x_reverse() +
  scale_y_reverse()

plot_transition(fit, "country_name", "Nigeria", n = 50) +
  scale_x_reverse() +
  scale_y_reverse()

P_tilde <- tidybayes::spread_draws(fit$samples, P_tilde[c], P_tilde2[c]) %>%
  left_join(fit$country_index)

scaled_transition <- fit$posteriors$transition$country_name %>%
  left_join(P_tilde, by = c(".chain", ".iteration", name = "country_name")) %>%
  mutate(x = P_tilde + P_tilde2 * x)

scaled_transition <- scaled_transition %>%
  left_join(distinct(scaled_transition, .chain, .iteration) %>%
              mutate(index = 1:n()))

scaled_transition_summary <- scaled_transition %>%
  mutate(x = round(x, 1)) %>%
  group_by(name, x) %>%
  median_qi(Y, .width = c(0.5, 0.8, 0.95))

plot_scaled_transition <- function(country) {
  indices <- sample(unique(scaled_transition$index), 25)
  scaled_transition_summary %>%
    filter(name %in% country) %>%
    ggplot(aes(x = x, y = Y)) +
    ggdist::geom_lineribbon(aes(ymin = .lower, ymax = .upper), size = 0.5) +
    geom_line(data = filter(scaled_transition, name %in% country, index %in% indices), aes(group = index), alpha = 0.25) +
    geom_point(data = rename(dat, name = country_name) %>% filter(name %in% country), aes(x = tfr, y = tfr_change)) +
    scale_fill_brewer() +
    facet_wrap(~name) 
}

plot_scaled_transition("Burkina Faso")

tidybayes::spread_draws(fit$samples, transition_function[c,t]) %>%
  left_join(fit$country_index) %>%
  left_join(fit$time_index) %>%
  group_by(country_name, period) %>%
  median_qi(transition_function, .width = c(0.5, 0.8, 0.95)) %>%
  ggplot(aes(x = period, y = transition_function)) +
  ggdist::geom_lineribbon(aes(ymin = .lower, ymax = .upper), size = 0.5) +
  scale_fill_brewer() +
  facet_wrap(~country_name)

tidybayes::spread_draws(fit$samples, a[c, i]) %>%
  left_join(fit$country_index) %>%
  group_by(country_name, i) %>%
  median_qi(a, .width = c(0.5, 0.8, 0.95)) %>%
  ggplot(aes(x = i, y = a)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  scale_fill_brewer() +
  facet_wrap(~country_name)

plot_hierarchical_sd(fit, "splines")

compare_transitions <- function(country, overlap = FALSE) {
  country_data <- filter(dat, country_name == country)
  
  bayestfr_transition <- extract_bayesTFR_transition("D:/spline_rate_model_paper/bayesTFR.output/", country)
  bayestfr_transition_summary <- bayestfr_transition %>%
    group_by(x) %>% 
    mutate(decrement = -decrement) %>%
    median_qi(decrement, .width = c(0.5, 0.8, 0.95)) %>%
    mutate(model = "BayesTFR")
  
  bayestfr_draws <- sample(unique(bayestfr_transition$draw), 10)
  spline_indices <- sample(unique(scaled_transition$index), 10)
  
  if(overlap == TRUE) {
    scaled_transition_summary %>%
      mutate(model = "Splines") %>%
      filter(name %in% country, .width == 0.95) %>%
      rename(decrement = Y) %>%
      bind_rows(bayestfr_transition_summary %>% filter(.width == 0.95)) %>%
      ggplot(aes(x = x, y = decrement)) +
      geom_line(aes(color = model), size = 1) +
      geom_line(aes(color = model, y = .lower)) +
      geom_line(aes(color = model, y = .upper)) +
      geom_point(data = country_data, aes(x = tfr, y = tfr_change))
  }
  else {
    scaled_transition_summary %>%
      mutate(model = "Splines") %>%
      filter(name %in% country) %>%
      rename(decrement = Y) %>%
      bind_rows(bayestfr_transition_summary) %>%
      ggplot(aes(x = x, y = decrement)) +
      ggdist::geom_lineribbon(aes(ymin = .lower, ymax = .upper), size = 0.5) +
      geom_line(data = filter(bayestfr_transition, draw %in% bayestfr_draws) %>% mutate(model = "BayesTFR"), aes(x = x, y = -decrement, group = draw), alpha = 0.25) +
      geom_line(data = filter(scaled_transition, name %in% country, index %in% spline_indices) %>% mutate(model = "Splines"), aes(y = Y, group = index), alpha = 0.25) +
      scale_fill_brewer() +
      geom_point(data = country_data, aes(x = tfr, y = tfr_change)) +
      facet_wrap(~model) +
      labs(x = "TFR", y = "Change in TFR") +
      scale_x_reverse() +
      scale_y_reverse() +
      cowplot::theme_cowplot()
  }
}

eta_pred2 <- spread_draws(fit$samples$draws("eta_pred2"), eta_pred2[c,t]) %>%
  left_join(fit$country_index) %>%
  left_join(fit$time_index) %>%
  left_join(fit$phases)

compare_projections <- function(country, overlap = TRUE) {
  pred <- get.tfr.prediction("D:/spline_rate_model_paper/bayesTFR.output")
  
  proj <- tfr.trajectories.table(pred, country) %>%
    as_tibble(rownames = "period") %>%
    mutate(period = as.numeric(period))
  
  if(overlap == TRUE) {
    eta_pred2 %>%
      filter(t >= start_phase2, country_name == country) %>%
      filter(!is.nan(eta_pred2)) %>%
      group_by(period, country_name) %>%
      median_qi(eta_pred2, .width = c(0.5, 0.8, 0.95)) %>%
      ggplot(aes(x = period, y = eta_pred2)) +
      geom_lineribbon(aes(ymin = .lower, ymax = .upper), size = 0.1) +
      geom_line(data = proj, aes(x = period, y = `0.975`, color = "BayesTFR 95% CI"), lty = 2) +
      geom_line(data = proj, aes(x = period, y = `0.025`, color = "BayesTFR 95% CI"), lty = 2) +
      geom_line(data = proj, aes(x = period, y = median, color = "BayesTFR Median"), lty = 2) +
      geom_point(data = filter(dat, country_name == country), aes(y = tfr), size = 1) +
      scale_fill_brewer() +
      scale_color_manual(values = c("red", "black")) +
      labs(x = "Period", y = "TFR") +
      cowplot::theme_cowplot()
  }
  else {
    bayestfr_traj <- get.tfr.trajectories(pred, country)
    bayestfr_sample <- bayestfr_traj[, sample(1:ncol(bayestfr_traj), 25)]
    
    spline_draws <- eta_pred2 %>%
      ungroup() %>%
      distinct(.chain, .iteration) %>%
      mutate(index = 1:n())
    
    spline_traj <- eta_pred2 %>%
      filter(country_name == country) %>%
      ungroup() %>%
      left_join(spline_draws, by = c(".chain", ".iteration")) %>%
      filter(index %in% sample(unique(spline_draws$index), 25)) %>%
      mutate(model = "Splines") %>%
      filter(!is.nan(eta_pred2))
    
    as_tibble(bayestfr_sample) %>%
      mutate(period = as.numeric(rownames(bayestfr_sample)), model = "BayesTFR") %>%
      filter(period > 2020) %>%
      pivot_longer(starts_with("V"), names_prefix = "V", names_transform = as.numeric, names_to = "draw", values_to = "tfr") %>%
      ggplot(aes(x = period, y = tfr)) +
      geom_line(aes(group = draw)) +
      geom_line(aes(y = eta_pred2, group = index), data = spline_traj) +
      facet_wrap(~model)
  }
}

country <- "Burkina Faso"
compare_projections(country, TRUE)
compare_transitions(country, FALSE)

fit$posteriors$generated_quantities %>%
  left_join(mutate(dat, i = 1:n())) %>%
  ungroup() %>%
  sample_n(1e4) %>%
  ggplot(aes(x = tfr, y = resid)) +
  geom_point() +
  geom_smooth()


# Plot point estimates at 2050 and 2010
compare_point_estimates <- function(period, label_cutoff = 2.25) {
  period_ <- period
  splines <- fit$posteriors$temporal %>%
    filter(variable == "eta", period == period_, country_name %in% countries_phase2$country_name) %>%
    select(period, country_name, `50%`) %>%
    rename(spline_median = `50%`)
  
  pred <- get.tfr.prediction("D:/spline_rate_model_paper/bayesTFR.output")
  
  proj <- map(unique(dat$country_name), function(x) {
    tfr.trajectories.table(pred, x) %>% 
      as_tibble(rownames = "period") %>%
      mutate(period = as.numeric(period),
             country_name = x) %>%
      filter(country_name %in% countries_phase2$country_name) %>%
      select(period, median, country_name)
  }) %>%
    bind_rows()
  
  plot_data <- left_join(splines, proj)
  
  plot_data %>%
    ggplot(aes(x = median, y = spline_median)) +
    geom_point() +
    ggrepel::geom_label_repel(aes(label = country_name), data = filter(plot_data, spline_median > label_cutoff)) +
    geom_abline(slope = 1, lty = 2) +
    labs(x = paste0("BayesTFR median ", period),
         y = paste0("Splines median ", period))
}

compare_point_estimates(2048, label_cutoff = 3.5) +
  ggtitle("Comparison 2048", subtitle = "Phase II countries only")
ggsave("D:/spline_rate_model_paper/plots/point_estimate_comparison_2048.pdf", width = 8, height = 8)

compare_point_estimates(2098, label_cutoff = 2.25) +
  ggtitle("Comparison 2098", subtitle = "Phase II countries only")
ggsave("D:/spline_rate_model_paper/plots/point_estimate_comparison_2098.pdf", width = 8, height = 8)

compare_projections("Rwanda")
compare_transitions("Rwanda")

pdf("D:/spline_rate_model_paper/plots/bayestfr_splines_comparison_phase2.pdf", width = 18, height = 6)
for(country in sort(unique(countries_phase2$country_name))) {
  print(country)
  p1 <- compare_transitions(country)
  p2 <- compare_projections(country)
  
  p <- p1 + p2 &
    plot_annotation(title = country)
  print(p)
}
dev.off()

p1 <- compare_transitions("Angola", overlap = FALSE) + ggtitle("Angola") + theme(legend.position = "none")
p2 <- compare_transitions("Afghanistan", overlap = FALSE) + ggtitle("Afghanistan")
p3 <- compare_transitions("Cameroon", overlap = FALSE) + ggtitle("Cameroon") + theme(legend.position = "none")
p4 <- compare_transitions("Congo", overlap = FALSE) + ggtitle("Congo") + theme(legend.position = "none")
p5 <- compare_transitions("Chad", overlap = FALSE) + ggtitle("Chad") + theme(legend.position = "none")
p6 <- compare_transitions("Niger", overlap = FALSE) + ggtitle("Niger") + theme(legend.position = "none")

(p1 + p2) / (p3 + p4) / (p5 + p6)
ggsave("D:/spline_rate_model_paper/plots/tfr_transition_examples.pdf", width = 12, height = 10)

p1 <- compare_projections("Angola") + ggtitle("Angola") + theme(legend.position = "none")
p2 <- compare_projections("Afghanistan") + ggtitle("Afghanistan")
p3 <- compare_projections("Cameroon") + ggtitle("Cameroon") + theme(legend.position = "none")
p4 <- compare_projections("Congo") + ggtitle("Congo") + theme(legend.position = "none")
p5 <- compare_projections("Chad") + ggtitle("Chad") + theme(legend.position = "none")
p6 <- compare_projections("Niger") + ggtitle("Niger") + theme(legend.position = "none")

(p1 + p2) / (p3 + p4) / (p5 + p6)
ggsave("D:/spline_rate_model_paper/plots/tfr_projection_examples.pdf", width = 10, height = 10)
