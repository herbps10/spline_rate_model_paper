# B-spline Transition Models

Code for "Flexible Modeling of Transition Processes with B-Splines" by Herbert Susmann & Leontine Alkema.

The principal dependency of this code is the `BayesTransitionModels` R package.

To run the main mCPR models:
```
library(targets)
tar_make()
```

Or to run the main models in parallel:
```
library(targets)
tar_make_future()
```

Code to generate plots and tables can be found in the `R/figures` and
`R/tables` folders. The plot generation code assumes that the final logistic
and spline fits can be found in `server-output/final_spline` and
`server-output/final_logistic`, respectively, which may require copying from
the `targets` output folder.

The main code for running the TFR case study can be found in `R/wpp.R`.

This work was supported, in whole or in part, by the Bill & Melinda Gates Foundation (INV-00844).
