library(tidyverse)
library(patchwork)

prior_predictive_splines_plot <- function() {
  knots <- 7
  beta_world <- rnorm(knots - 1, 0, 1)
  sigma_world <- truncnorm::rtruncnorm(1, 0, 1, a = 0, b = Inf)
  sigma_region <- truncnorm::rtruncnorm(1, 0, 1, a = 0, b = Inf)
  
  grid <- seq(0, 1, 0.01)
  regions <- 2
  subregions <- 2
  countries_per_subregion <- 5
  
  beta_region <- tibble(
    r = 1:regions,
    beta_region = rerun(regions, rnorm(knots - 1, mean = beta_world, sd = sigma_world)),
    y = map(beta_region, function(beta) {
      y <- spline_transition_function(0.01 + 0.29 * boot::inv.logit(beta), num_knots = knots)(grid)
      tibble(x = grid, y = y)
    })
  )
  
  beta_subregion <- expand_grid(
    r = 1:regions,
    s = 1:subregions
  ) %>%
    left_join(beta_region) %>%
    mutate(
      subregion_index = (r - 1) * subregions + s,
      beta_subregion = map(beta_region, function(beta_region) {
        rnorm(knots - 1, mean = beta_region, sd = sigma_region)
      }),
      y = map(beta_subregion, function(beta) {
        y <- spline_transition_function(0.01 + 0.29 * boot::inv.logit(beta), num_knots = knots)(grid)
        tibble(x = grid, y = y)
      })
    )
  
  prior_predictive <- expand_grid(
    r = 1:regions,
    s = 1:subregions,
    c = 1:countries_per_subregion,
  ) %>%
    left_join(beta_subregion) %>%
    mutate(
      index = (subregion_index - 1) * countries_per_subregion + c,
      beta_country = map(beta_subregion, function(beta_subregion) {
        rnorm(knots - 1, mean = beta_subregion, sd = sigma_region)
      }),
      y = map(beta_country, function(beta) {
        y <- spline_transition_function(0.01 + 0.29 * boot::inv.logit(beta), num_knots = knots)(grid)
        tibble(x = grid, y = y)
      })
    ) %>%
    unnest(y)
  
  p1 <- prior_predictive %>%
    ggplot(aes(x = x, y = y, group = index)) +
    geom_line(data = tibble(index = 0, x = grid, y = spline_transition_function(0.01 + 0.29 * boot::inv.logit(beta_world), num_knots = knots)(grid))) +
    cowplot::theme_cowplot() +
    theme(legend.position = "none", plot.title = element_text(face = "italic")) +
    labs(x = expression(eta[ct]/P[c]^u), y = expression(f[b]), title = "Global")
  
  p2 <- prior_predictive %>%
    ggplot(aes(x = x, y = y, group = index)) +
    geom_line(data = unnest(beta_region, y), aes(group = r, lty = factor(r))) +
    cowplot::theme_cowplot() +
    theme(legend.position = "none", plot.title = element_text(face = "italic")) +
    labs(x = expression(eta[ct]/P[c]^u), y = expression(f[b]), title = "Regional")
  
  p3 <- prior_predictive %>%
    ggplot(aes(x = x, y = y, group = index)) +
    geom_line(data = unnest(beta_subregion, y), aes(group = subregion_index, lty = factor(r), color = factor(subregion_index))) +
    cowplot::theme_cowplot() +
    theme(legend.position = "none", plot.title = element_text(face = "italic")) +
    labs(x = expression(eta[ct]/P[c]^u), y = expression(f[b]), title = "Subregional")
  
  p4 <- prior_predictive %>%
    ggplot(aes(x = x, y = y, group = index, lty = factor(r), color = factor(subregion_index))) +
    geom_line() +
    cowplot::theme_cowplot() +
    theme(legend.position = "none", plot.title = element_text(face = "italic")) +
    labs(x = expression(eta[ct]/P[c]^u), y = expression(f[b]), title = "Area-specific")
  
  p1 + p2 + p3 + p4
}

set.seed(20)
p1 <- prior_predictive_splines_plot()
p1

ggsave("plots/prior_predictive_transition_functions.pdf", height = 8, width = 8)
