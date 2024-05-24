# Install modified bayesTFR from Github that allows
# disabling phase-3 estimation
#
# remotes::install_github("herbps10/bayesTFR")

library(bayesTFR)

mc <- run.tfr.mcmc(output.dir = "bayesTFR.output", iter = 8e3, nr.chains = 2, replace = TRUE)
mc <- get.tfr.mcmc(sim.dir = "bayesTFR.output")

par(mfrow = c(1, 2))
bayesTFR::DLcurve.plot(mc, country = "Thailand", nr.curves = 100)
bayesTFR::DLcurve.plot(mc, country = "Zambia", nr.curves = 100)
par(mfrow = c(1, 1))

pred <- tfr.predict(
  mc, 
  output.dir = "bayesTFR.output", 
  replace.output = TRUE,
  use.tfr3 = FALSE,
  mu = 0,
  rho = 0,
  sigmaAR1 = 0,
  nr.traj = 8e3,
  enable_phase3 = FALSE
)

#for(cutoff in c(2003, 2008, 2013)) {
tfr_results <- tibble(
  cutoff = c(2003, 2008, 2013)
) %>%
  mutate(res = map(cutoff, function(cutoff) {
    output.dir <- glue::glue("bayesTFR.output-{cutoff}")
    mc <- run.tfr.mcmc(output.dir = output.dir, present.year = cutoff, iter = 8e3, nr.chains = 2, replace = TRUE)
    
    pred <- tfr.predict(
      mc, 
      output.dir = output.dir,
      replace.output = TRUE,
      use.tfr3 = FALSE,
      mu = 0,
      rho = 0,
      sigmaAR1 = 0,
      nr.traj = 5e2,
      burnin = 0,
      enable_phase3 = FALSE, 
      end.year = 2025
    )
    
    res <- tibble(
      country_name = mc$meta$regions$country_name,
    ) %>%
      mutate(projections = map(country_name, function(country) {
        s <- summary(pred, country = country)$projections
        tibble(
          q2.5 = s["2018", "2.5%"],
          q5 = s["2018", "5%"],
          q10 = s["2018", "10%"],
          q25 = s["2018", "25%"],
          q50 = s["2018", "50%"],
          q75 = s["2018", "75%"],
          q95 = s["2018", "95%"],
          q97.5 = s["2018", "97.5%"]
        )
      })) %>%
      left_join(filter(dat, period == 2018) %>% select(country_name, tfr)) %>%
      filter(!is.na(tfr))
    
    res
  }))

tfr_validation <- tfr_results %>% unnest(res) %>% unnest(projections) %>%
  mutate(
    error = tfr - q50, 
    covered_95 = q2.5 <= tfr & q97.5 >= tfr, 
    covered_80 = q10 <= tfr & q90 >= tfr, 
    below10 = tfr < q10, 
    above90 = tfr > q90) %>%
  group_by(cutoff) %>%
  summarise(mean_error = mean(error), median_error = median(error), median_absolute_error = median(abs(error)), coverage95 = mean(covered_95), coverage80 = mean(covered_80), below10 = mean(below10), above90 = mean(above90), ci_width = mean(q90 - q10),
            interval_score = mean(q90 - q10 + 2/0.2 * (q10 - tfr) * (tfr < q10) + 2/0.2 * (tfr - q90) * (tfr > q90))) %>%
  mutate(model = "BayesTFR", num_knots = NA)


tfr_validation_by_area <- tfr_results %>% unnest(res) %>% unnest(projections) %>%
  left_join(areas) %>%
  mutate(error = tfr - q50, covered_95 = q2.5 <= tfr & q97.5 >= tfr, covered_80 = q10 <= tfr & q90 >= tfr, below10 = tfr < q10, above90 = tfr > q90) %>%
  group_by(area_name, cutoff) %>%
  summarise(mean_error = mean(error), median_error = median(error), median_absolute_error = median(abs(error)), coverage95 = mean(covered_95), coverage80 = mean(covered_80), below10 = mean(below10), above90 = mean(above90), ci_width = mean(q90 - q10),
            interval_score = mean(q90 - q10 + 2/0.2 * (q10 - tfr) * (tfr < q10) + 2/0.2 * (tfr - q90) * (tfr > q90))) %>%
  mutate(model = "BayesTFR", num_knots = NA)
