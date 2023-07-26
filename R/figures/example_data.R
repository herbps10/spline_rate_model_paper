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

# Rate of change plot


plot_countries_rate_of_change <- function(data, countries) data %>%
  group_by(name_country, year) %>%
  summarize(contraceptive_use_modern = mean(contraceptive_use_modern)) %>%
  group_by(name_country) %>%
  mutate(diff = c(diff(contraceptive_use_modern), NA), diff_time = c(diff(year), NA), avg_diff = diff / diff_time) %>%
  ggplot(aes(x = contraceptive_use_modern, y = avg_diff)) +
  geom_point(aes(), size = 1.25) +
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  cowplot::theme_cowplot() +
  cowplot::panel_border() +
  theme(axis.text = element_text(size = 12)) +
  facet_wrap(~name_country, ncol = 4) +
  labs(x = "mCPR", y = "Change in mCPR")

p1 <- data %>% filter(
  name_country %in% c("Bangladesh", "Indonesia", "Kenya", "China", 
                      "Viet Nam", "Ghana", "Thailand", "Peru", 
                      "Philippines", "Dominican Republic", "Senegal", "Burkina Faso")
) %>%
  plot_countries_rate_of_change() +
  theme(axis.title.x = element_blank())

p1 + geom_smooth() + facet_wrap(~name_country, scales = "free")

p2 <- data %>% filter(
  name_country %in% c("Swaziland", "South Africa", "Guinea")
) %>%
  plot_countries_rate_of_change()

(p1 / p2) +
  plot_layout(guides = "collect")

#
# TFR example data
#
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

dat %>%
  filter(country_name %in% c("India", "Kenya", "Mozambique", "Pakistan", "Thailand", "Zimbabwe")) %>%
  ggplot(aes(x = period, y = tfr)) +
  geom_point() +
  facet_wrap(~country_name) +
  cowplot::theme_cowplot() +
  cowplot::panel_border() +
  theme(axis.text = element_text(size = 12)) +
  labs(x = "Period", y = "TFR") +
  theme(plot.margin = margin(0.1, 1, 0.1, 0.1, "cm"))

ggsave("plots/tfr_plots.pdf", height = 5, width = 8)

# Rate of change version
dat %>%
  filter(country_name %in% c("India", "Kenya", "Mozambique", "Pakistan", "Thailand", "Zimbabwe")) %>%
  group_by(country_name) %>%
  mutate(diff = c(diff(tfr), NA)) %>%
  ggplot(aes(x = tfr, y = -diff)) +
  geom_point() +
  facet_wrap(~country_name) +
  cowplot::theme_cowplot() +
  cowplot::panel_border() +
  theme(axis.text = element_text(size = 12)) +
  labs(x = "Period", y = "TFR") +
  theme(plot.margin = margin(0.1, 1, 0.1, 0.1, "cm"))

ggsave("plots/tfr_plots.pdf", height = 5, width = 8)
