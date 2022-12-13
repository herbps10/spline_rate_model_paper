library(tidyverse)
library(tidybayes)
library(BayesTransitionModels)
library(patchwork)

source("R/extract_tfr_data.R")

compare_transitions <- function(country, overlap = FALSE) {
  country_data <- filter(dat, country_name == country)
  
  bayestfr_transition <- extract_bayesTFR_transition("bayesTFR.output/", country)
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


compare_projections <- function(country, overlap = TRUE) {
  pred <- get.tfr.prediction("bayesTFR.output")
  
  proj <- tfr.trajectories.table(pred, country) %>%
    as_tibble(rownames = "period") %>%
    mutate(period = as.numeric(period))
  
  if(overlap == TRUE) {
    eta %>%
      filter(t >= start_phase2, country_name == country) %>%
      filter(!is.nan(eta)) %>%
      group_by(period, country_name) %>%
      median_qi(eta, .width = c(0.5, 0.8, 0.95)) %>%
      ggplot(aes(x = period, y = eta)) +
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
    
    spline_draws <- eta %>%
      ungroup() %>%
      distinct(.chain, .iteration) %>%
      mutate(index = 1:n())
    
    spline_traj <- eta %>%
      filter(country_name == country) %>%
      ungroup() %>%
      left_join(spline_draws, by = c(".chain", ".iteration")) %>%
      filter(index %in% sample(unique(spline_draws$index), 25)) %>%
      mutate(model = "Splines") %>%
      filter(!is.nan(eta))
    
    as_tibble(bayestfr_sample) %>%
      mutate(period = as.numeric(rownames(bayestfr_sample)), model = "BayesTFR") %>%
      filter(period > 2020) %>%
      pivot_longer(starts_with("V"), names_prefix = "V", names_transform = as.numeric, names_to = "draw", values_to = "tfr") %>%
      ggplot(aes(x = period, y = tfr)) +
      geom_line(aes(group = draw)) +
      geom_line(aes(y = eta, group = index), data = spline_traj) +
      facet_wrap(~model)
  }
}


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

countries <- c("Angola", "Afghanistan", "Cameroon", 
               "Congo", "Burkina Faso", "Kenya", 
               "Mozambique", "South Africa", "Botswana",
               "Cuba", "Mongolia", "United Kingdom")

#countries <- c("Cuba", "United Kingdom", "Kenya", "Mongolia")
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
  adapt_delta = 0.999,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 14623
)

fit_summary <- fit$samples$summary(NULL, posterior::default_convergence_measures())

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

eta <- spread_draws(fit$samples$draws("eta"), eta[c,t]) %>%
  left_join(fit$country_index) %>%
  left_join(fit$time_index) %>%
  left_join(fit$phases) %>%
  filter(t > end_phase2)

pdf("plots/bayestfr_splines_comparison_phase2.pdf", width = 18, height = 6)
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
ggsave("plots/tfr_transition_examples.pdf", width = 12, height = 10)

p1 <- compare_projections("Angola") + ggtitle("Angola") + theme(legend.position = "none")
p2 <- compare_projections("Afghanistan") + ggtitle("Afghanistan")
p3 <- compare_projections("Cameroon") + ggtitle("Cameroon") + theme(legend.position = "none")
p4 <- compare_projections("Congo") + ggtitle("Congo") + theme(legend.position = "none")
p5 <- compare_projections("Chad") + ggtitle("Chad") + theme(legend.position = "none")
p6 <- compare_projections("Niger") + ggtitle("Niger") + theme(legend.position = "none")

(p1 + p2) / (p3 + p4) / (p5 + p6)
ggsave("plots/tfr_projection_examples.pdf", width = 10, height = 10)
