library(tidyverse)
library(targets)
library(BayesTransitionModels)
library(patchwork)

final_spline <- readRDS("server-output/final_spline")
final_logistic <- readRDS("server-output/final_logistic")

my_theme <- cowplot::theme_cowplot() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

#
# Spline
#

p <- plot_indicator(final_spline, c("Bangladesh", "Indonesia", "Kenya",
                               "Guinea", "South Africa", "Swaziland"),
               color_sources = TRUE) +
  guides(fill = "none") +
  scale_y_continuous(limits = c(0, 0.9)) +
  scale_x_continuous(limits = c(1970, 2030)) +
  labs(x = "Year") +
  my_theme  +
  theme(panel.spacing = unit(0.04, "npc"),
        axis.line = element_line())

p$layers[[5]]$mapping <- aes(ymin = lower, ymax = upper, y = contraceptive_use_modern)
p$layers[[6]]$mapping <- aes(y = contraceptive_use_modern, shape = source)
p$facet$params$free$x <- TRUE
p$facet$params$free$y <- TRUE
  
print(p)
ggsave("plots/example_fits.pdf", width = 8, height = 5)


#
# Logistic
#
p <- plot_indicator(final_logistic, c("Bangladesh", "Indonesia", "Kenya",
                               "Guinea", "South Africa", "Swaziland"),
               color_sources = TRUE) +
  guides(fill = "none") +
  scale_y_continuous(limits = c(0, 0.9)) +
  scale_x_continuous(limits = c(1970, 2030)) +
  labs(x = "Year") +
  my_theme  +
  theme(panel.spacing = unit(0.04, "npc"),
        axis.line = element_line())

p$layers[[5]]$mapping <- aes(ymin = lower, ymax = upper, y = contraceptive_use_modern)
p$layers[[6]]$mapping <- aes(y = contraceptive_use_modern, shape = source)
p$facet$params$free$x <- TRUE
p$facet$params$free$y <- TRUE
  
print(p)
ggsave("plots/example_logistic_fits.pdf", width = 8, height = 5)
