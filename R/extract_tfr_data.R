library(bayesTFR)
library(tidyverse)

extract_bayesTFR_transition <- function(sim.dir, country) {
  mcmc <- get.tfr.mcmc(sim.dir)
  mcmc.list <- get.mcmc.list(mcmc)
  meta <- mcmc.list[[1]]$meta
  country <- get.country.object(country, meta)
  xs <- seq(0, 10, 0.1)
  dlcurves <- tfr.get.dlcurves(xs, mcmc.list, country$code, country$index, 0, NULL, FALSE)
  
  bayesTFR:::get.observed.tfr(country$index, meta, "tfr_matrix_observed", 
                              "tfr_matrix_all")[1:meta$T_end_c[country$index]]
  
  dec <- as_tibble(t(dlcurves)) %>%
    mutate(x = xs) %>%
    pivot_longer(starts_with("V"), names_prefix = "V", names_to = "draw", values_to = "decrement", names_transform = as.numeric)
  
  dec
}

extract_bayesTFR_projection <- function(sim.dir, country) {
  mcmc <- get.tfr.mcmc(sim.dir)
  mcmc.list <- get.mcmc.list(mcmc)
  meta <- mcmc.list[[1]]$meta
  country <- get.country.object(country, meta)
  xs <- seq(0, 10, 0.1)
  dlcurves <- tfr.get.dlcurves(xs, mcmc.list, country$code, country$index, 0, NULL, FALSE)
  
  bayesTFR:::get.observed.tfr(country$index, meta, "tfr_matrix_observed", 
                              "tfr_matrix_all")[1:meta$T_end_c[country$index]]
  
  dec <- as_tibble(t(dlcurves)) %>%
    mutate(x = xs) %>%
    pivot_longer(starts_with("V"), names_prefix = "V", names_to = "draw", values_to = "decrement", names_transform = as.numeric)
  
  dec
}

tfr_dataset <- function() {
  # following the example code of run mcmc to get the meta associated with default wpp2019 run
  sim.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
  m <- get.tfr.mcmc(sim.dir)

  # relevant information, stored as vectors with 201 countries:
  # (after finding out about tfr_matrix, this is no longer necessary)
  m$meta$start_c # transition starts at or prior to this index
  m$meta$lambda_c # name suggests that transition ends at this index -1, but given lots of 14s,
  # interpretation should be that this index still needs to be included for rate of change in phase 2

  # country info is here
  m$meta$regions$country_code
  m$meta$regions$country_name

  # the UN estimates (input data) are here
  m$meta$tfr_matrix # this has NAs outside phase 2 (hence start and lambda no longer needed)
  dim(m$meta$tfr_matrix)
  # quick check that codes indeed match up
  cbind(colnames(m$meta$tfr_matrix_observed), m$meta$regions$country_code)

  # combining information
  # actually with NAs in tfr_matrix, no need to store lambda and start
  country_info <- tibble(
    country_code = as.character(m$meta$regions$country_code),
    country_name = m$meta$regions$country_name,
    start_phase2 = m$meta$start_c,
    start_phase3 =  m$meta$lambda_c
  )
  country_data <- t(m$meta$tfr_matrix)#_observed)

  tfr_inputs <- country_info %>%
    left_join(as_tibble(country_data) %>% mutate(country_code = rownames(country_data)))

  tfr_inputs
}
