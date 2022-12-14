# Install modified bayesTFR from Github that allows
# disabling phase-3 estimation
#
# remotes::install_github("herbps10/bayesTFR")

library(bayesTFR)

mc <- run.tfr.mcmc(output.dir = "bayesTFR.output", iter = 8e3, nr.chains = 2, replace = TRUE)
mc <- get.tfr.mcmc(sim.dir = "bayesTFR.output")

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
