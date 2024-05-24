library(tidyverse)
library(tidybayes)
library(BayesTransitionModels)
library(patchwork)

source("R/extract_tfr_data.R")

compare_transitions <- function(scaled_transition, scaled_transition_summary, bayestfr_dir, country, overlap = FALSE) {
  country_data <- filter(dat, country_name == country)
  
  bayestfr_transition <- extract_bayesTFR_transition(bayestfr_dir, country)
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


compare_projections <- function(eta, country, bayestfr_dir, overlap = TRUE) {
  pred <- get.tfr.prediction(bayestfr_dir)
  
  proj <- tfr.trajectories.table(pred, country) %>%
    as_tibble(rownames = "period") %>%
    mutate(period = as.numeric(period)) %>%
    filter(period <= 2038)
  
  if(overlap == TRUE) {
    eta %>%
      filter(t >= start_phase2, country_name == country) %>%
      filter(!is.nan(eta)) %>%
      filter(period <= 2038) %>%
      group_by(period, country_name) %>%
      median_qi(eta, .width = c(0.5, 0.8, 0.95)) %>%
      ggplot(aes(x = period, y = eta)) +
      geom_lineribbon(aes(ymin = .lower, ymax = .upper), size = 0.1) +
      geom_line(data = proj, aes(x = period, y = `0.9`, color = "BayesTFR 80% CI"), lty = 2) +
      geom_line(data = proj, aes(x = period, y = `0.1`, color = "BayesTFR 80% CI"), lty = 2) +
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
      xlim(c(1950, 2035)) +
      facet_wrap(~model)
  }
}

summarize <- function(fit) {
  #fit_summary <- fit$samples$summary(NULL, posterior::default_convergence_measures())
  
  P_tilde <- tidybayes::spread_draws(fit$samples, P_tilde[c], P_tilde2[c]) %>%
    left_join(fit$country_index)
  
  scaled_transition <- fit$posterior$transition$country_name %>%
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
  
  list(
    P_tilde = P_tilde,
    eta = eta,
    scaled_transition = scaled_transition,
    scaled_transition_summary = scaled_transition_summary
  )
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
               "Cuba", "Mongolia")

#countries <- c("Cuba", "United Kingdom", "Kenya", "Mongolia")
#dat <- dat %>% filter(country_name %in% countries)
#dat <- dat %>% filter(country_name %in% countries_phase2$country_name)

dat %>%
  ggplot(aes(x = period, y = tfr)) +
  geom_point() +
  facet_wrap(~country_name)

fit_tfr <- function(num_knots, cutoff) {
  tfrplus(
    dat,
    held_out = dat$period > cutoff,
    start_year = 1953,
    end_year = 2100,
    y = "tfr",
    year = "period",
    area = "country_name",
    num_knots = num_knots,
    spline_degree = 2,
    hierarchical_splines = c("intercept", "country_name"),
    max_treedepth = 12,
    adapt_delta = 0.999,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 500,
    seed = 14623
  )
}

fit <- fit_tfr(4, Inf) 
plot_indicator(fit, europe)
s <- summarize(fit)

countries <- c("Paraguay", "Cambodia", "Yemen", "Kenya")

p1 <- plot_indicator(fit2, countries) + ggtitle("2 knots") + facet_wrap(~country_name, nrow = 1)
p2 <- plot_indicator(fit3, countries) + ggtitle("3 knots") + facet_wrap(~country_name, nrow = 1)
p3 <- plot_indicator(fit7, countries) + ggtitle("7 knots") + facet_wrap(~country_name, nrow = 1)
p1 / p2 / p3

p1 <- plot_transition(fit2, "country_name", "Yemen") + ggtitle("2 knots") + facet_wrap(~name, nrow = 1)
p2 <- plot_transition(fit3, "country_name", "Yemen") + ggtitle("3 knots") + facet_wrap(~name, nrow = 1)
p3 <- plot_transition(fit7, "country_name", "Yemen") + ggtitle("7 knots") + facet_wrap(~name, nrow = 1)
p1 / p2 / p3

fit3$posteriors$temporal %>% filter(period == 2018) %>% select(`50%`, country_name) %>%
  left_join(fit7$posteriors$temporal %>% filter(period == 2018) %>% select(`50%`, country_name), by = c("country_name")) %>%
  mutate(diff = `50%.x` - `50%.y`) %>%
  arrange(diff)

residuals <- function(fit) {
  left_join(fit$posteriors$temporal, mutate(fit$data, held_out = fit$held_out)) %>%
    filter(held_out == 1) %>%
    group_by(country_name) %>%
    filter(period == 2018) %>%
    mutate(residual = tfr - `50%`, covered_95 = `2.5%` < tfr & `97.5%` > tfr)
}

fits <- fits %>% bind_rows(expand_grid(
  num_knots = c(2:7),
  cutoff = c(2018)
) %>%
  mutate(
    resid = map2(num_knots, cutoff, function(num_knots, cutoff) {
      fit_tfr(num_knots, cutoff) %>%
        residuals()
    })
  ))
  
validation <- fits %>%
  #mutate(resid = map(fit, residuals)) %>%
  #select(-fit) %>%
  unnest(resid) %>%
  filter(period == 2018) %>%
  filter(country_name %in% countries_phase2$country_name) %>%
  group_by(num_knots, cutoff) %>%
  mutate(
    error = tfr - `50%`, 
    below10 = tfr < `10%`, 
    below25 = tfr < `25%`,
    below50 = tfr < `50%`,
    above50 = tfr < `50%`,
    above75 = tfr > `75%`,
    above90 = tfr > `90%`, 
    covered_80 = `10%` <= tfr & `90%` >= tfr, 
    covered_95 = `2.5%` <= tfr & `97.5%` >= tfr
  ) %>%
  summarise(
    mean_error = mean(error), 
    median_error = median(error), 
    median_absolute_error = median(abs(error)), 
    coverage95 = mean(covered_95), 
    coverage80 = mean(covered_80), 
    below10 = mean(below10), 
    below25 = mean(below25),
    below50 = mean(below50),
    above50 = mean(above50),
    above75 = mean(above75),
    above90 = mean(above90), 
    ci_width = mean(`90%` - `10%`),
    interval_score = mean(`90%` - `10%` + 2/0.2 * (`10%` - tfr) * (tfr < `10%`) + 2/0.2 * (tfr - `90%`) * (tfr > `90%`))) %>%
  arrange(cutoff, num_knots) %>%
  mutate(model = "Spline")

validation %>% 
  bind_rows(tfr_validation) %>%
  arrange(cutoff) %>%
  ungroup() %>%
  mutate(name = glue::glue("{model} (K = {num_knots})")) %>%
  group_by(cutoff) %>%
  select(name, below10, coverage80, above90, ci_width, median_error, median_absolute_error) %>%
  mutate_at(vars(ci_width, median_error, median_absolute_error), scales::number_format(scale = 100, accuracy = 0.1)) %>%
  #mutate_at(vars(ci_width, median_error, median_absolute_error), star_best) %>%
  #mutate_at(vars(interval_score), scales::number_format(scale = 1, accuracy = 0.01)) %>%
  mutate_at(vars(below10, coverage80, above90), scales::percent_format(accuracy = 0.01)) %>%
  ungroup() %>%
  select(-cutoff) %>%
  knitr::kable(format = "latex") %>%
  str_replace_all("\\n\\\\hline", "")  %>%
  cat()

validation %>% 
  bind_rows(tfr_validation) %>%
  arrange(cutoff) %>%
  ungroup() %>%
  mutate(name = glue::glue("{model} (K = {num_knots})")) %>%
  group_by(cutoff) %>%
  select(name, below10, below25, below50, above50, above75, above90) %>%
  mutate_at(vars(below10, below25, below50, above50, above75, above90), scales::percent_format(accuracy = 0.1)) %>%
  ungroup() %>%
  select(-cutoff) %>%
  knitr::kable(format = "latex") %>%
  str_replace_all("\\n\\\\hline", "")  %>%
  cat()

sim.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
m <- get.tfr.mcmc(sim.dir)
areas <- tibble(
  country_name = m$meta$regions$country_name,
  area_name = m$meta$regions$area_name
)

validation_by_area <- fits %>%
  unnest(resid) %>%
  filter(country_name %in% countries_phase2$country_name) %>%
  left_join(areas) %>%
  filter(period == 2018) %>%
  group_by(num_knots, cutoff, area_name) %>%
  mutate(error = tfr - `50%`, covered_80 = `10%` <= tfr & `90%` >= tfr, below10 = tfr < `10%`, above90 = tfr > `90%`, covered_95 = `2.5%` <= tfr & `97.5%` >= tfr) %>%
  summarise(n = n(), mean_error = mean(error), median_error = median(error), median_absolute_error = median(abs(error)), coverage95 = mean(covered_95), coverage80 = mean(covered_80), above90 = mean(above90), below10 = mean(below10), ci_width = mean(`90%` - `10%`),
            interval_score = mean(`90%` - `10%` + 2/0.2 * (`10%` - tfr) * (tfr < `10%`) + 2/0.2 * (tfr - `90%`) * (tfr > `90%`))) %>%
  arrange(cutoff, num_knots) %>%
  mutate(model = "Spline")

remove_dups <- function(x) {
  x <- as.character(x)
  x[x == lag(x)] <- ""
  x
}

star_best <- function(x) { 
  x[abs(as.numeric(x)) == min(abs(as.numeric(x)))] <- paste0(x[abs(as.numeric(x)) == min(abs(as.numeric(x)))], "*") 
  x
}

validation_by_area %>% 
  bind_rows(tfr_validation_by_area) %>%
  arrange(cutoff) %>%
  ungroup() %>%
  mutate(name = glue::glue("{model} (K = {num_knots})")) %>%
  arrange(cutoff, area_name, num_knots) %>%
  group_by(cutoff, area_name) %>%
  select(cutoff, area_name, name, below10, coverage80, above90, ci_width, median_error, median_absolute_error) %>%
  mutate_at(vars(ci_width, median_error, median_absolute_error), scales::number_format(scale = 100, accuracy = 0.1)) %>%
  #mutate_at(vars(ci_width, median_error, median_absolute_error), star_best) %>%
  #mutate_at(vars(interval_score), scales::number_format(scale = 1, accuracy = 0.01)) %>%
  mutate_at(vars(below10, coverage80, above90), scales::percent_format(accuracy = 0.01)) %>%
  mutate_at(vars(cutoff, area_name), remove_dups) %>%
  ungroup() %>%
  select(-cutoff) %>%
  knitr::kable(format = "latex") %>%
  str_replace_all("\\n\\\\hline", "")  %>%
  cat()

summaries <- fits %>%
  mutate(summary = map(fit, summarize))

pdf("plots/tfr_knot_comparisons.pdf", width = 8, height = 7)
for(country in sort(unique(dat$country_name))) {
  print(country)
  p1 <- compare_transitions(
    summaries$summary[[1]]$scaled_transition,
    summaries$summary[[1]]$scaled_transition_summary, 
    country, 
    overlap = FALSE
  ) +
    ggtitle(glue::glue("{country}: 4 knots"))
  
  p2 <- compare_transitions(
    summaries$summary[[2]]$scaled_transition,
    summaries$summary[[2]]$scaled_transition_summary, 
    country, 
    overlap = FALSE
  ) +
    ggtitle(glue::glue("{country}: 7 knots"))
  
  print(p1 / p2)
}
dev.off()


pdf("plots/tfr_knot_projection_comparisons.pdf", width = 8, height = 7)
for(country in sort(unique(dat$country_name))) {
  print(country)
  p1 <- compare_projections(
    summaries$summary[[1]]$eta,
    country
  ) +
    ggtitle(glue::glue("{country}: 4 knots"))
  
  p2 <- compare_projections(
    summaries$summary[[2]]$eta,
    country
  ) +
    ggtitle(glue::glue("{country}: 7 knots"))
  
  print(p1 / p2)
}
dev.off()

pdf("plots/bayestfr_splines_comparison_phase2.pdf", width = 18, height = 6)
for(country in sort(unique(countries_phase2$country_name))) {
  print(country)
  p1 <- compare_transitions(s$scaled_transition, s$scaled_transition_summary, "bayesTFR.output", country, overlap = FALSE) + ggtitle(country) + theme(legend.position = "none")
  p2 <- compare_projections(s$eta, country, "bayesTFR.output") + ggtitle(country) + theme(legend.position = "none")
  
  p <- p1 + p2 &
    plot_annotation(title = country)
  print(p)
}
dev.off()

p1 <- compare_transitions(s$scaled_transition, s$scaled_transition_summary, "bayesTFR.output", "Angola", overlap = FALSE) + ggtitle("Angola") + theme(legend.position = "none")
p2 <- compare_transitions(s$scaled_transition, s$scaled_transition_summary, "bayesTFR.output", "Afghanistan", overlap = FALSE) + ggtitle("Afghanistan")
p3 <- compare_transitions(s$scaled_transition, s$scaled_transition_summary, "bayesTFR.output", "Cameroon", overlap = FALSE) + ggtitle("Cameroon") + theme(legend.position = "none")
p4 <- compare_transitions(s$scaled_transition, s$scaled_transition_summary, "bayesTFR.output", "Congo", overlap = FALSE) + ggtitle("Congo") + theme(legend.position = "none")
p5 <- compare_transitions(s$scaled_transition, s$scaled_transition_summary, "bayesTFR.output", "Chad", overlap = FALSE) + ggtitle("Chad") + theme(legend.position = "none")
p6 <- compare_transitions(s$scaled_transition, s$scaled_transition_summary, "bayesTFR.output", "Zambia", overlap = FALSE) + ggtitle("Zambia") + theme(legend.position = "none")

(p1 + p2) / (p3 + p4) / (p5 + p6)
ggsave("plots/tfr_transition_examples_k4.pdf", width = 12, height = 10)

p1 <- compare_projections(s$eta, "Angola", "bayesTFR.output") + ggtitle("Angola") + theme(legend.position = "none")
p2 <- compare_projections(s$eta, "Afghanistan", "bayesTFR.output") + ggtitle("Afghanistan")
p3 <- compare_projections(s$eta, "Cameroon", "bayesTFR.output") + ggtitle("Cameroon") + theme(legend.position = "none")
p4 <- compare_projections(s$eta, "Congo", "bayesTFR.output") + ggtitle("Congo") + theme(legend.position = "none")
p5 <- compare_projections(s$eta, "Chad", "bayesTFR.output") + ggtitle("Chad") + theme(legend.position = "none")
p6 <- compare_projections(s$eta, "Zambia", "bayesTFR.output") + ggtitle("Zambia") + theme(legend.position = "none")

(p1 + p2) / (p3 + p4) / (p5 + p6)
ggsave("plots/tfr_projection_examples_k4.pdf", width = 10, height = 10)

compare_projections(s$eta, "Senegal", "bayesTFR.output") + ggtitle("Senegal") + theme(legend.position = "none")
