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

pdf("plots/spline_results.pdf", width = 15, height = 5)
for(country in unique(final_spline$data$name_country)) {
  print(country)
  p1 <- plot_transition(final_spline, "name_country", areas = country) +
    labs(x = expression(eta[ct] / lambda^u[c]), y = expression(f[b])) +
    guides(fill = "none") +
    my_theme
    
  p2 <- plot_smoother(final_spline, areas = country) +
    labs(x = "Year") +
    guides(fill = "none") +
    my_theme
  
  p3 <- plot_indicator(final_spline, areas = country, color_sources = TRUE) +
    scale_shape_discrete(name = "Data Source", drop = FALSE) +
    labs(x = "Year") +
    guides(fill = "none") +
    my_theme
  
  p3$layers[[5]]$mapping <- aes(ymin = lower, ymax = upper, y = contraceptive_use_modern)
  p3$layers[[6]]$mapping <- aes(y = contraceptive_use_modern, shape = data_series_type)
  
  if(country == "RÈunion") country <- "Réunion"
  
  p <- (p1 + p2 + p3) + 
    plot_annotation(title = country, tag_levels = "A", theme = theme(plot.title = element_text(face = "bold"))) &
    theme(plot.tag = element_text(face = "plain"))
  print(p)
}
dev.off()
