# remotes::install_github("herbps10/bayesTFR")
library(bayesTFR)

mc <- run.tfr.mcmc(output.dir = "D:/spline_rate_model_paper/bayesTFR.output", iter = 8e3, nr.chains = 2, replace = TRUE)
mc <- get.tfr.mcmc(sim.dir = "D:/spline_rate_model_paper/bayesTFR.output")

pred <- tfr.predict(
  mc, 
  output.dir = "D:/spline_rate_model_paper/bayesTFR.output_test", 
  replace.output = TRUE,
  use.tfr3 = FALSE,
  mu = 0,
  rho = 0,
  sigmaAR1 = 0,
  nr.traj = 8e3,
  enable_phase3 = FALSE
)
