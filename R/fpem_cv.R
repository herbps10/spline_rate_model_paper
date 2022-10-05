cv_fit <- function(data, held_out, ...) {
  fpemplus(
      data,
      y = "contraceptive_use_modern",
      se = "se_modern",
      year = "year",
      source = "data_series_type",
      area = "name_country",
      held_out = held_out, 
      
      start_year = 1970, 
      end_year = 2030, 
      
      # Hierarchical setup
      hierarchical_level     = c("intercept", "name_region", "name_sub_region", "name_country"), 
      hierarchical_splines   = c("intercept", "name_region", "name_sub_region", "name_country"),
      hierarchical_asymptote = c("intercept", "name_region", "name_sub_region", "name_country"),
      
      # Stan sampler settings
      parallel_chains = 4,
      refresh = 50, 
      ...
    )
}

cv_fit_cutoff <- function(data, cutoff, ...) {
  held_out <- data$year >= cutoff
  cv_fit(data, held_out, ...) 
}

cv_fit_random <- function(data, seed, prop = 0.2, reps = 5, ...) {
  set.seed(seed)
  tibble(
    rep = 1:reps
  ) %>%
    mutate(held_out = map(rep, function(x) {
      indices <- sample(1:nrow(data), size = round(nrow(data) * prop))
      held_out <- rep(FALSE, nrow(data))
      held_out[indices] <- TRUE
      held_out
    })) %>%
    mutate(cv_fit = map(held_out, function(held_out) {
      cv_fit(data, held_out, ...)
    }))
}