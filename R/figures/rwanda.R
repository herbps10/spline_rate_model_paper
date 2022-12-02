library(tidyverse)
library(targets)
library(fpemplus)
library(patchwork)

final_logistic <- readRDS("server-output/final_logistic")
final_spline <- readRDS("server-output/final_spline")

my_theme <- cowplot::theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

area <- "Rwanda"

#
# Logistic
#
p1 <- plot_transition(final_logistic, "name_country", areas = area) +
  guides(fill = "none") +
  ylim(c(0, 0.25)) +
  labs(x = expression(eta[ct] / lambda[c]^u), y = expression(f[b])) +
  my_theme +
  ggtitle("", subtitle = "Logistic")

p2 <- plot_smoother(final_logistic, areas = area) +
  guides(fill = "none") +
  ylim(c(-0.3, 0.35)) +
  scale_x_continuous(breaks = seq(1970, 2030, 10)) +
  labs(x = "Year", y = expression(epsilon[ct])) +
  my_theme +
  ggtitle("", subtitle = "Logistic")

p3 <- plot_indicator(final_logistic, areas = area, color_sources = TRUE) +
  scale_shape_discrete(name = "Data Source", drop = FALSE) +
  guides(fill = "none") +
  labs(x = "Year", y = expression(eta[ct])) +
  scale_x_continuous(breaks = seq(1970, 2030, 10)) +
  ylim(c(0, 0.8)) +
  my_theme +
  ggtitle("", subtitle = "Logistic")

p3$layers[[5]]$mapping <- aes(ymin = lower, ymax = upper, y = contraceptive_use_modern)
p3$layers[[6]]$mapping <- aes(y = contraceptive_use_modern, shape = data_series_type)


#
# Spline
#
p4 <- plot_transition(final_spline, "name_country", areas = area) +
  guides(fill = "none") +
  ylim(c(0, 0.25)) +
  labs(x = expression(eta[ct] / lambda[c]^u), y = expression(f[b])) +
  my_theme +
  ggtitle("", subtitle = "B-spline (d=2, K=5)")

p5 <- plot_smoother(final_spline, areas = area) +
  guides(fill = "none") +
  ylim(c(-0.3, 0.35)) +
  labs(x = "Year", y = expression(epsilon[ct])) +
  my_theme +
  ggtitle("", subtitle = "B-spline (d=2, K=5)")

p6 <- plot_indicator(final_spline, areas = area, color_sources = TRUE) +
  scale_shape_discrete(name = "Data Source", drop = FALSE) +
  guides(fill = "none") +
  labs(x = "Year", y = expression(eta[ct])) +
  scale_x_continuous(breaks = seq(1970, 2030, 10)) +
  ylim(c(0, 0.8)) +
  my_theme +
  ggtitle("", subtitle = "B-spline (d=2, K=5)")

p6$layers[[5]]$mapping <- aes(ymin = lower, ymax = upper, y = contraceptive_use_modern)
p6$layers[[6]]$mapping <- aes(y = contraceptive_use_modern, shape = data_series_type)

p <- (p1 + p4) / (p2 + p5) / (p3 + p6) + 
  plot_layout(guides = "collect") &
  plot_annotation(
    title = area,
    tag_levels = "A", 
    theme = theme(plot.title  = element_text(face = "bold"))
  ) &
  theme(plot.tag = element_text(face = "plain"))

print(p)
ggsave("plots/rwanda_example.pdf", width = 10, height = 12)
