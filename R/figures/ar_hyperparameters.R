library(tidyverse)
library(targets)
library(fpemplus)
library(patchwork)

final_spline <- readRDS("server-output/final_spline")
final_logistic <- readRDS("server-output/final_logistic")

my_theme <- cowplot::theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

bind_rows(
  final_spline$posteriors$ar %>%
    mutate(x = est_tau / sqrt(1 - est_rho^2)) %>% 
    mutate(model = "B-spline (d=2, K=5)") %>%
    select(model, x),
  final_logistic$posteriors$ar %>%
    mutate(x = est_tau / sqrt(1 - est_rho^2)) %>%
    mutate(model = "Logistic") %>%
    select(model, x)
) %>%
  ggplot(aes(x = x, color = model)) +
  geom_density() +
  my_theme +
  labs(x = expression(tau/sqrt(1-rho^2))) +
  theme(legend.position = "bottom")

ggsave("plots/ar_logistic_spline_comparison.pdf", width = 10, height = 7)
