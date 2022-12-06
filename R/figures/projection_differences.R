library(tidyverse)
library(targets)
library(BayesTransitionModels)
library(patchwork)

final_spline <- readRDS("server-output/final_spline")
final_logistic <- readRDS("server-output/final_logistic")

my_theme <- cowplot::theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

eta2030 <- bind_rows(
  final_spline$posteriors$temporal %>%
    filter(variable == "eta", year == 2030) %>%
    select(name_country, `2.5%`, `50%`, `97.5%`) %>%
    mutate(model = "B-spline"),
  final_logistic$posteriors$temporal %>%
    filter(variable == "eta", year == 2030) %>%
    select(name_country, `2.5%`, `50%`, `97.5%`) %>%
    mutate(model = "Logistic")
)

diffs <- eta2030 %>%
  select(name_country, `50%`, model) %>%
  pivot_wider(values_from = `50%`, names_from = model) %>%
  mutate(abs_diff = abs(`B-spline` - Logistic)) %>%
  arrange(-abs_diff)

diffs$abs_diff %>% hist()

eta2030 %>%
  left_join(diffs) %>%
  filter(name_country %in% diffs$name_country[1:15]) %>%
  ggplot(aes(x = reorder(name_country, abs_diff), y = `50%`, color = model)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  labs(x = "Country", y = "mCPR in 2030") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_text(angle = 90))

ggsave("plots/eta_diffs2030.pdf", width = 8, height = 6)
