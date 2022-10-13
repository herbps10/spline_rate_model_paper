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

sub_regions <- final_spline$posteriors$transition_quantiles[["name_sub_region"]] %>%
  left_join(distinct(final_spline$data, name_region, name_sub_region), by = c(name = "name_sub_region")) %>%
  rename(name_sub_region = name)

regions <- final_spline$posteriors$transition_quantiles[["name_region"]] %>%
  left_join(distinct(final_spline$data, name_region, name_sub_region), by = c(name = "name_region")) %>%
  rename(name_region = name)

countries <- final_spline$posteriors$transition_quantiles[["name_country"]] %>%
  left_join(distinct(final_spline$data, name_country, name_region, name_sub_region), by = c(name = "name_country")) %>%
  rename(name_country = name)

ggplot(sub_regions, aes(x = x, y = Y)) +
  geom_line(data = countries, aes(color = "Country", group = name_country), size = 0.5, alpha = 0.5) +
  geom_line(aes(color = "Sub-region"), size = 0.8, alpha = 0.8) +
  geom_line(data = regions, aes(color = "Region"), size = 1, alpha = 0.5) +
  scale_color_manual(values = c("#7f8c8d", "#2980b9", "#f39c12")) +
  facet_wrap(~fct_reorder(factor(name_sub_region), name_region)) +
  my_theme +
  labs(y = expression(f[b]), x = expression(eta[ct]/P^u[c]))

ggsave("plots/transition_function_hierarchy.pdf", width = 14, height = 10)
