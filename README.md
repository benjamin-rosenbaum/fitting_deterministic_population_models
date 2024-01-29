# Tutorial: Fitting deterministic population models in R using RStan

A tutorial accompanying the publication:

B. Rosenbaum, E.A. Fronhofer. 2023. Confronting population models with experimental microcosm data: from trajectory matching to state-space models. Ecosphere 14(4): e4503. [https://doi.org/10.1002/ecs2.4503](https://doi.org/10.1002/ecs2.4503)

For questions, contact <a href= "mailto:benjamin.rosenbaum@idiv.de">benjamin.rosenbaum@idiv.de</a>

## Update 01/2024

Raw R-code is now available in a seperate folder.

I updated deprecated Stan code, concerning arrays and the ODE-solver function. 

While I wrote R code using the popular [RStan](https://mc-stan.org/users/interfaces/rstan) package, I strongly recommend using the [CmdStanR](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) package instead. This increases performance and I observed model fitting is <b>three times faster</b>. While the Stan model code is identical, R functions for compiling, fitting and analyzing the model are slightly different. An example is now given below. 

## Single-species system

[Logistic growth - OBS model](https://benjamin-rosenbaum.github.io/fitting_deterministic_population_models/logistic_obs.html)

[Logistic growth - PROC model](https://benjamin-rosenbaum.github.io/fitting_deterministic_population_models/logistic_proc.html)

[Logistic growth - State-space model](https://benjamin-rosenbaum.github.io/fitting_deterministic_population_models/logistic_ssm.html)

## Two-species system

[Consumer resource - OBS model](https://benjamin-rosenbaum.github.io/fitting_deterministic_population_models/consumer_resource_obs.html)

[Consumer resource - PROC model](https://benjamin-rosenbaum.github.io/fitting_deterministic_population_models/consumer_resource_proc.html)

[Consumer resource - State-space model](https://benjamin-rosenbaum.github.io/fitting_deterministic_population_models/consumer_resource_ssm.html)

## CmdStanR

[Logistic growth - OBS model cmdstanr](https://benjamin-rosenbaum.github.io/fitting_deterministic_population_models/logistic_obs_cmdstanr.html)

