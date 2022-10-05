library(tidyverse)
library(fpemplus)
library(targets)

tar_option_set(
  packages = c("tidyverse", "fpemplus")
)

analysis_data_target <- tar_target(analysis_data, national_data(
  countries = c(),
  start_year = 1970
) %>% 
  filter(name_country != "Other non-specified areas") %>%
  filter(modern_method_bias == "None", has_geographical_region_bias == "N") %>%
  mutate(name_region = case_when(
    name_region == "Europe" ~ "Europe and Northern America",
    name_region == "Northern America" ~ "Europe and Northern America",
    TRUE ~ name_region
  )) 
)

final_spline_target <- tar_target(final_spline, fpemplus(
  analysis_data,
  y = "contraceptive_use_modern",
  se = "se_modern",
  year = "year",
  source = "data_series_type",
  area = "name_country",
  
  start_year = 1970, 
  end_year = 2030, 
  
  model = "spline",
  
  # Spline setup
  spline_degree = 2,
  num_knots = 5,
  
  # Hierarchical setup
  hierarchical_level     = c("intercept", "name_region", "name_sub_region", "name_country"), 
  hierarchical_splines   = c("intercept", "name_region", "name_sub_region", "name_country"),
  hierarchical_asymptote = c("intercept", "name_region", "name_sub_region", "name_country"),
  
  # Prior settings
  tau_prior = "normal(0, 2)",
  rho_prior_sd = "uniform(0, 1)",
  
  # Stan sampler settings
  adapt_delta = 0.999,
  max_treedepth = 14,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 1482395,
  parallel_chains = 4,
  refresh = 50
))

list(
  analysis_data_target,
  final_spline_target
)
