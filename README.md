# CO2projections

This code is designed for generating probabilistic CO2 projections in paper _Less Than 2C Warming by 2100 Unlikely_

## Code Usage
This code is mostly meant to be run sequentially. 

For use of this code, one should load all required packages in the first 30 lines, and load all functionals in constructing the model from line 32 to 1785, including population projection, fit and project the complete model, and evaluation of prediction by summarizing statistics in projections.

Users should also load lines between line 1786 to 1998 for data loading, cleaning and loading of countries.

Population projections are conducted via function **_predict.population_** and the predicted trajectories are extracted via **_find.pop.trajectories_**. For example,
```R
preds.countries <- predict.population(year.present=2010, year.end=2100,
                                               make.new=TRUE)
pred.pop <- preds.countries$pred.pop
names.countries.tmp <- as.character(unique(subset(data.medium, Isocode != "USA")$Isocode))
preds.countries.trajs <- find.pop.trajectories(data.medium, pred.pop, n.trajectories=1000,
                                               c("USA", names.countries.tmp), year.start=2010, year.end=2100)

```

This model is estimated via MCMC and **JAGS software** is used in the estimation process. This is wrapped by function **_fit.project.model_** and summary statistics are grabbed by function **_evaluate.projections_**. Before running these functions, please make sure that file **corr_model_ar1trend_const.bug** for describing models for JAGS is in the directory, and for predicting CO2 emissions from year a to b, please predict population with the above codes and save them with name **"poppreds_formatted_a_b.rda"**. For example, if we want to predict CO2 emissions with this model from year 2010 to 2100, we should have **"poppreds_formatted_2010_2100.rda"** in _NatureData/_ directory, and use:
```R
model.ar1.const.results.2010 <- fit.project.model(data.medium, fit.ar1.const=T)
proj.evals.2010.ar1.const <- evaluate.projections(data.medium, model.ar1.const.results.2010,
                                                  year.start=2010, year.end=2100,
                                                  outofsample.validate=F,
                                                  quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975))

```

Codes between line 2104 to 3145 are for creating plots in the paper. It includes:
..* Past and Present GDP plots by countries.
..* Current Carbon Intensities.
..* Yearly Carbon Emissions Predictions (World and Sub-Saharan Areas).
..* Cumulative Carbon Emissions Predictions (World and Sub-Saharan Areas).
..* IPAT Components Projections (Intensity, Population and GDP per capita)
..* Out of Sample Validations
..* Projections of Carbon Emissions by Regions.

Codes between line 3146 to 3248 are for modeling for temperature changes based on Figure SPM.5 of AR5 Summary Report

Codes between line 3249 to 3635 are for predicting the probability that 6 main countries and regions meet their objective in Paris Agreement, including China, India, European Union, United States, Russia and Japan.

Codes between line 3636 to 3743 are for plots in supplementary data, including predictions for 15 countries.
