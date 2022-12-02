library(tidyverse)
library(patchwork)

final_spline <- readRDS("server-output/final_spline")

data <- final_spline$data %>%
  mutate(`Data Source` = factor(data_series_type)) 

plot_countries <- function(data, countries) data %>%
  ggplot(aes(x = year, y = contraceptive_use_modern)) +
  geom_point(aes(shape = `Data Source`), size = 1.25) +
  scale_shape_discrete(drop = FALSE) +
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  cowplot::theme_cowplot() +
  cowplot::panel_border() +
  theme(axis.text = element_text(size = 12)) +
  facet_wrap(~name_country, ncol = 4) +
  labs(x = "Survey Year", y = "mCPR") +
  ylim(c(0, 0.7)) +
  scale_x_continuous(limits = c(1970, 2020), breaks = c(1980, 2000, 2020), labels = c(1980, 2000, 2020))

p1 <- data %>% filter(
  name_country %in% c("Bangladesh", "Indonesia", "Kenya")
) %>%
  plot_countries() +
  theme(axis.title.x = element_blank())

p2 <- data %>% filter(
  name_country %in% c("Swaziland", "South Africa", "Guinea")
) %>%
  plot_countries()

(p1 / p2) +
  plot_layout(guides = "collect")

#ggsave("plots/example_data_countries.pdf", height = 8, width = 12)
ggsave("plots/example_data_countries.pdf", height = 5, width = 8)

