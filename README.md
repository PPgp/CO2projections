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

This model is estimated via MCMC and **JAGS software** is used in the estimation process. 
