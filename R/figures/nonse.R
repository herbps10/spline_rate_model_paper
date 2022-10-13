library(tidyverse)
library(targets)
library(fpemplus)
library(patchwork)
library(tidybayes)

final_spline <- readRDS("server-output/final_spline")

my_theme <- cowplot::theme_cowplot() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

final_spline$posteriors$nonse %>%
  group_by(data_series_type) %>%
  median_qi(nonse) %>%
  ggplot(aes(x = nonse, y = reorder(data_series_type, nonse))) +
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0) +
  geom_point() +
  my_theme +
  labs(x = expression(sigma[s]~(posterior~non-sampling~error)), y = "")

ggsave("plots/nonse.pdf", width = 8, height = 8)