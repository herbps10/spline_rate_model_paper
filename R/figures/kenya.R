library(tidyverse)
library(targets)
library(fpemplus)
library(patchwork)

final_spline <- readRDS("server-output/final_spline")

my_theme <- cowplot::theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

p1 <- plot_transition(final_spline, "name_country", areas = "Kenya") +
  guides(fill = "none") +
  labs(x = expression(eta[ct] / lambda[c]^u), y = expression(f[b])) +
  my_theme

p2 <- plot_smoother(final_spline, areas = "Kenya") +
  guides(fill = "none") +
  labs(x = "Year", y = expression(epsilon[ct])) +
  scale_x_continuous(breaks = seq(1970, 2030, 10)) +
  my_theme

p3 <- plot_indicator(final_spline, areas = "Kenya", color_sources = TRUE) +
  scale_shape_discrete(name = "Data Source", drop = FALSE) +
  guides(fill = "none") +
  labs(x = "Year", y = expression(eta[ct])) +
  scale_x_continuous(breaks = seq(1970, 2030, 10)) +
  my_theme

p3$layers[[5]]$mapping <- aes(ymin = lower, ymax = upper, y = contraceptive_use_modern)
p3$layers[[6]]$mapping <- aes(y = contraceptive_use_modern, shape = data_series_type)
p3

p <- (p1 / p2 / p3) + 
  plot_layout(guides = "collect") &
  plot_annotation(
    title = "Kenya",
    tag_levels = "A", 
    theme = theme(plot.title  = element_text(face = "bold"))
  ) &
  theme(plot.tag = element_text(face = "plain"))

print(p)

ggsave("plots/spline_result_example.pdf", width = 8, height = 10)
