library(tidyverse)
library(fpemplus)
library(targets)
library(tarchetypes)
library(future)

source("R/fpem_cv.R")

plan(multicore)

tar_option_set(
  packages = c("tidyverse", "fpemplus")
)

#
# Global settings
#
tau_prior <- "normal(0, 2)"
rho_prior <- "uniform(0, 1)"

cv_adapt_delta      <- 0.99
cv_max_treedepth    <- 12
cv_iter_warmup      <- 250
cv_iter_sampling    <- 500

#print("TESTING MODE")
#cv_adapt_delta      <- 0.95
#cv_max_treedepth    <- 8
#cv_iter_warmup      <- 250
#cv_iter_sampling    <- 250


final_adapt_delta   <- 0.999
final_max_treedepth <- 14
final_iter_warmup   <- 500
final_iter_sampling <- 500

output_dir <- "/work/hsusmann_umass_edu/spline_rate_model_paper/output/"
#output_dir <- "d:/spline_rate_model_paper/output"

countries <- c()
#countries <- c("Somalia", "India", "Palau", "Bangladesh", "Zimbabwe",
#             "Indonesia", "Kenya", "Turkey", "Mexico", "Guatemala")

analysis_data_target <- tar_target(analysis_data, national_data(
  countries = countries,
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

cv_setup <- tribble(
  ~model,     ~spline_degree, ~num_knots,
  "spline",   2,              5,
  "spline",   2,              7,
  "spline",   3,              5,
  "spline",   3,              7,
  "logistic", 2,              7
)

cv_fits_target <- tar_map(
  unlist = FALSE,
  values = cv_setup,
  tar_target(cv_fit_cutoff2010, cv_fit_cutoff(
    data = analysis_data,
    cutoff = 2010,
    logistic = model == "logistic",
    
    spline_degree = spline_degree,
    num_knots = num_knots,
    
    tau_prior = tau_prior,
    rho_prior = rho_prior,
    
    adapt_delta   = cv_adapt_delta,
    max_treedepth = cv_max_treedepth,
    iter_warmup   = cv_iter_warmup,
    iter_sampling = cv_iter_sampling,
    seed = ifelse(model == "logistic", 100, 0) + spline_degree + num_knots
  )),
  tar_target(cv_fit_cutoff2015, cv_fit_cutoff(
    data = analysis_data,
    cutoff = 2015,
    logistic = model == "logistic",
    
    spline_degree = spline_degree,
    num_knots = num_knots,
    
    tau_prior = tau_prior,
    rho_prior = rho_prior,
    
    adapt_delta   = cv_adapt_delta,
    max_treedepth = cv_max_treedepth,
    iter_warmup   = cv_iter_warmup,
    iter_sampling = cv_iter_sampling,
		output_dir = output_dir,
    seed = ifelse(model == "logistic", 150, 0) + spline_degree + num_knots
  )),
  tar_target(cv_fit_random5, cv_fit_random(
    data = analysis_data,
    seed = ifelse(model == "logistic", 200, 0) + spline_degree + num_knots,
    prop = 0.2,
    reps = 2,
    
    spline_degree = spline_degree,
    num_knots = num_knots,
    
    tau_prior = tau_prior,
    rho_prior = rho_prior,
    
    adapt_delta   = cv_adapt_delta,
    max_treedepth = cv_max_treedepth,
    iter_warmup   = cv_iter_warmup,
    iter_sampling = cv_iter_sampling,
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
  hierarchical_asymptote = c("intercept", "name_country"),
  
  # Prior settings
  tau_prior = "normal(0, 2)",
  rho_prior = "uniform(0, 1)",
  
  # Stan sampler settings
  adapt_delta = final_adapt_delta,
  max_treedepth = final_max_treedepth,
  iter_warmup = final_iter_warmup,
  iter_sampling = final_iter_sampling,
	output_dir = output_dir,
  seed = 1482395,
  parallel_chains = 4,
  refresh = 50
))

final_logistic_target <- tar_target(final_logistic, fpemplus_logistic(
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
  hierarchical_rate      = c("intercept", "name_region", "name_sub_region", "name_country"),
  hierarchical_asymptote = c("intercept", "name_country"),
  
  # Prior settings
  tau_prior = "normal(0, 2)",
  rho_prior = "uniform(0, 1)",
  
  # Stan sampler settings
  adapt_delta = final_adapt_delta,
  max_treedepth = final_max_treedepth,
  iter_warmup = final_iter_warmup,
  iter_sampling = final_iter_sampling,
  output_dir = output_dir,
  seed = 1482395,
  parallel_chains = 4,
  refresh = 50
))

list(
  analysis_data_target,
  cv_fits_target
  #final_spline_target
)
