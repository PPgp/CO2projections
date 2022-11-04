
# Things to run before any session:

# miscellaneous
library(gdata)
library(plyr)
library(countrycode)
library(doBy)
library(reshape2)

# Population prediction libraries
library(wpp2015)
library(bayesTFR)
library(bayesLife)
library(bayesPop)
library(popReconstruct)

# For model fitting with JAGS
library(rjags)
library(coda)

# for plotting, displaying
library(xtable)
library(ggplot2)
library(gridExtra)
library(scales)
library(MASS)

library(car)

#=========================================================

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

find_maxes <- function(data.medium, plot.countries=FALSE) {
  # find_maxes finds when each country peaks in carbon intensity.
  # If the peak for a country is within the last 5 years of the data,
  # we list that country under rejects.late and don't report a peak
  # for it.
  # If there are fewer than 20 years of carbon intensity (Tech) data
  # for a country, we list that country under rejects.insuf and don't
  # report a peak for it.
  max.vals <- list()
  max.vals[["Iso"]] <- c()
  max.vals[["GDP"]] <- c()
  max.vals[["Tech"]] <- c()
  max.vals[["CO2"]] <- c()
  max.vals[["Year"]] <- c()
  rejects.early <- c()
  rejects.late <- c()
  rejects.insuf <- c()
  isolist <- unique(data.medium$Isocode)
  if (plot.countries) i <- 1
  for (iso in isolist) {
    data.tmp <- data.medium[data.medium$Isocode == iso,]
    data.tmp <- data.tmp[!is.na(data.tmp$Tech),]
    if (sum(!is.na(data.tmp$Tech)) < 20) {
      rejects.insuf <- c(rejects.insuf, iso)
      next
    }
    model.smoothed <- loess(Tech ~ Year, data=data.tmp, span=0.25)
    predicted <- predict(model.smoothed)
    max.ind <- which(predicted == max(predicted))
    stopifnot(length(max.ind) == 1)
    
    # Record unsmoothed observations at this point, if we're not near
    # the start or end of the period
    year.last <- max(data.medium$Year) # Allow for out-of-sample validation
    year.last.cutoff <- min(2003, year.last - 5)
    if (1965 <= data.tmp$Year[max.ind] & 
        data.tmp$Year[max.ind] <= year.last.cutoff) {
      max.vals[["Iso"]] <- c(max.vals[["Iso"]], iso)
      max.vals[["GDP"]] <- c(max.vals[["GDP"]], data.tmp$GDP[max.ind])
      max.vals[["Tech"]] <- c(max.vals[["Tech"]], data.tmp$Tech[max.ind])
      max.vals[["CO2"]] <- c(max.vals[["CO2"]], data.tmp$CO2[max.ind])
      max.vals[["Year"]] <- c(max.vals[["Year"]], data.tmp$Year[max.ind])
    } else if (data.tmp$Year[max.ind] <= 1964) {
      rejects.early <- c(rejects.early, iso)
    } else {
      rejects.late <- c(rejects.late, iso)
    }
    
    # Create and save plot
    if (plot.countries) {
      n.perpage <- 8
      if (i %% n.perpage == 1) {
        if (i != 1) dev.off()
        png(paste0(plot.location, "max_plots", floor(i/n.perpage), ".png"),
            width=800, height=1200)
        par(mfrow=c(4,2))
      }
      plot(data.tmp$Year, data.tmp$Tech,
           main=paste0(getCountry(iso), " (", iso, ")"))
      points(data.tmp$Year[!is.na(data.tmp$Tech)], predicted, type="l")
      i <- i + 1
    }
  }
  if (plot.countries) dev.off()
  
  print("Rejects because peak was before 1965:")
  print(rejects.early)
  print(length(rejects.early))
  print(paste(getCountry(rejects.early), collapse=", "))
  print(paste0("Rejects because peak was after ", min(2003, year.last - 5)))
  print(rejects.late)
  print(length(rejects.late))
  print(paste(getCountry(rejects.late), collapse=", "))
  print("Rejects because of insufficient data:")
  print(rejects.insuf)
  print(length(rejects.insuf))
  print(paste(getCountry(rejects.insuf), collapse=", "))
  
  list(max.vals, rejects.early, rejects.late, rejects.insuf)
}

#==========================================


predict.population <- function(year.present=2010, year.end=2100,
                               n.iter=10000, n.burnin=1000,
                               make.new=TRUE) {
  # predict.population uses the packages bayesPop, bayesTFR,
  # and bayesLife to create population projections by country,
  # saving results into sims.location.
  # Set make.new=FALSE if you've already run predict.population
  # for the years of interest to load previous results.
  library(bayesPop)
  if (make.new) {
    library(bayesTFR)
    library(bayesLife)
    sim.dir.tfr <- paste0(sims.location, "Sims_TFR_", year.present, "_", year.end)
    m.tfr <- run.tfr.mcmc(present.year=year.present, 
                          iter=n.iter,
                          burnin=n.burnin,
                          output.dir=sim.dir.tfr, verbose=FALSE,
                          replace.output=TRUE)
    if (year.present != 1980) {
      m3.tfr <- run.tfr3.mcmc(present.year=year.present,
                              sim.dir=sim.dir.tfr,
                              iter=n.iter,
                              thin=1,
                              burnin=n.burnin,
                              verbose=TRUE, replace.output=TRUE)
      pred.tfr <- tfr.predict(end.year=year.end,
                              sim.dir=sim.dir.tfr, m.tfr,
                              burnin=n.burnin,
                              use.tfr3=FALSE,
                              verbose=TRUE,
                              replace.output=TRUE)
    } else {
      pred.tfr <- tfr.predict(end.year=year.end,
                              sim.dir=sim.dir.tfr, m.tfr,
                              burnin=n.burnin, burnin3=n.burnin,
                              verbose=TRUE,
                              replace.output=TRUE)
    }
    
    sim.dir.e0 <- paste0(sims.location, "Sims_E0_", year.present, "_", year.end)
    m.e0 <- run.e0.mcmc(present.year=year.present,
                        my.e0.file=paste0(sims.location, "aids_include.txt"),
                        iter=n.iter,
                        thin=1,
                        verbose=FALSE,
                        output.dir=sim.dir.e0,
                        replace.output=TRUE)
    
    if (year.present < 2000) {
      pred.e0 <- e0.predict(end.year=year.end, 
                            max.e0.eq1 = 78, max.e0.eq1.pred=78,
                            sim.dir=sim.dir.e0, m.e0, burnin=n.burnin, verbose=FALSE)
    } else {
      pred.e0 <- e0.predict(end.year=year.end, 
                            sim.dir=sim.dir.e0, m.e0, burnin=n.burnin, verbose=FALSE)
    }
    summary(pred.e0)
    
    sim.dir.pop <- paste0(sims.location, "Sims_Pop_", year.present, "_", year.end)
    pred.pop <- pop.predict(present.year=year.present, end.year=year.end,
                            output.dir=sim.dir.pop, verbose=FALSE, 
                            inputs = list(tfr.sim.dir=sim.dir.tfr, 
                                          e0F.sim.dir=sim.dir.e0, e0M.sim.dir='joint_'),
                            replace.output=TRUE)
    summary(pred.pop)
  } else {
    # Loading an already-made file
    sim.dir.pop <- paste0(sims.location, "Sims_Pop_", year.present, "_", year.end)
    pred.pop <- get.pop.prediction(sim.dir.pop)
  }
  
  # 9 for pred.pop$quantiles[,9,] is the 0.5 quantile
  preds.countries <- data.frame(cbind(as.character(pred.pop$countries[,2]),
                                      pred.pop$quantiles[,'0.5',]),
                                stringsAsFactors=F)
  names(preds.countries) <- c("Country", paste0("Pop", seq(year.present, year.end, by=5)))
  library(countrycode)
  preds.countries$Isocode <- getIso(preds.countries$Country)
  
  list(preds.countries=preds.countries, pred.pop=pred.pop)
}



find.pop.trajectories <- function(data.medium, pred.pop, n.trajs, isos,
                                  n.iter.predict.population=1000,
                                  year.start=2010, year.end=2100,
                                  make.new=TRUE) {
  # find.pop.trajectories converts the format of the results from
  # predict.population to a format we use in our later analyses.
  trajectories.seq <- seq(from=1, to=n.iter.predict.population, length.out=n.trajs)
  if (make.new) {
    countries.codes <- pred.pop$countries$code
    preds.countries <- list()
    index <- 0
    for (traj in trajectories.seq) {
      index <- index + 1
      preds.countries[[index]] <- data.frame(matrix(ncol=(2+(year.end-year.start)/5), 
                                                    nrow=length(isos)))
      # Picking arbitrary country (the 356 for India) for pulling variable names
      names(preds.countries[[index]]) <- c("Isocode", paste0("Pop",
                                                             names(bayesPop:::get.pop.trajectories(pred.pop, 356)$trajectories[, 1] )))
    }
    row.ind <- 0
    for (iso in isos) {
      row.ind <- row.ind + 1
      country.code <- countrycode(iso, "iso3c", "un")
      if (iso == "TWN") {
        country.code <- 158
      }
      country.trajs <- bayesPop:::get.pop.trajectories(pred.pop, country.code)$trajectories
      index <- 0
      for (traj in trajectories.seq) {
        index <- index + 1
        preds.countries[[index]][row.ind, ] <- c(iso, country.trajs[, traj])
      }
    }
    index <- 0
    for (traj in trajectories.seq) {
      index <- index + 1
      # Convert back to numeric. Silly R!
      preds.countries[[index]][,-1] <- sapply(preds.countries[[index]][,-1], as.numeric)
    }
    # browser()
    save(preds.countries,
         file=paste0(sims.location, "Sims_Pop_", year.start, "_", year.end, "/predictions/converted_pop_", 
                     year.start, "_", year.end, ".RDA"))
  } else {
    load(paste0(sims.location, "Sims_Pop_", year.start, "_", year.end, "/predictions/converted_pop_", 
                year.start, "_", year.end, ".RDA"))
  }
  
  preds.countries
}

convert.projections <- function(tech.projections.log, 
                                gdp.projections.log, frontier.projections.log,
                                tech.USA.projections.log) {
  # convert.projections converts intensity (tech) and GDP projections
  # out of the log scale.
  tech.projections <- cbind(tech.projections.log[,1], exp(tech.projections.log[,-1]))
  names(tech.projections) <- gsub("Log", "", names(tech.projections))
  names(tech.projections)[1] <- "Isocode"
  gdp.projections <- cbind(gdp.projections.log[,1], exp(gdp.projections.log[,-1]))
  names(gdp.projections) <- gsub("Log", "", names(gdp.projections))
  names(gdp.projections)[1] <- "Isocode"
  frontier.projections <- cbind(frontier.projections.log[,1], 
                                exp(frontier.projections.log[,-1]))
  names(frontier.projections) <- gsub("Log", "", names(frontier.projections))
  names(frontier.projections)[1] <- "Isocode"
  tech.USA.projections <- cbind(tech.USA.projections.log[,1], 
                                exp(tech.USA.projections.log[,-1]))
  names(tech.USA.projections) <- gsub("Log", "", names(tech.USA.projections))
  names(tech.USA.projections)[1] <- "Isocode"
  
  list(GDPFrontier=frontier.projections, GDPData=gdp.projections, 
       TechData=tech.projections, TechUSA=tech.USA.projections)
}


estimate.co2.projections <- function(tech.projections, gdp.projections, 
                                     frontier.projections, tech.USA.projections,
                                     pop.projections,
                                     year.range) {
  # estimate.co2.projections combines GDP, CO2 intensity (tech), and
  # population projections to get CO2 emissions.
  # The United States must be handled separately due to its unique
  # role in the GDP model.
  stopifnot(year.range[1] %% 5 == 0 & year.range[2] %% 5 == 0)
  
  # Deal with the USA specially
  pop.projections.USA <- pop.projections[pop.projections$Isocode == "USA", ]
  stopifnot(dim(pop.projections.USA)[1] == 1)
  
  isos.tmp <- tech.projections$Isocode
  # Population projections match up by country with tech and GDP projections
  pop.projections.sorted <- pop.projections[-1, ] # Remove the USA row
  stopifnot(all(pop.projections.sorted$Isocode == isos.tmp))
  
  stopifnot(all(gdp.projections$Isocode == isos.tmp))
  
  year.vals <- seq(year.range[1], year.range[2], by=5)
  emissions.projections <- data.frame(matrix(nrow=length(isos.tmp), ncol=(length(year.vals)+1)))
  names(emissions.projections)[1] <- "Isocode"
  emissions.projections$Isocode <- isos.tmp
  names(emissions.projections)[-1] <- paste0("CO2", year.vals)
  # Add in the US specially
  emissions.projections <- rbind(c("USA", rep(NA, length(year.vals))),
                                 emissions.projections)
  n.countries <- dim(emissions.projections)[1]
  for (year in year.vals) {
    # dividing by 10 comes from intensity which has $10,000 in the denominator,
    # population which is divided by a factor of 1,000, and CO2 is in megatonnes
    var.name <- paste0("CO2", year)
    emissions.projections[-1,var.name] <- gdp.projections[, paste0("GDP", year)] *
      tech.projections[, paste0("Tech", year)] *
      pop.projections.sorted[, paste0("Pop", year)] / 10
    emissions.projections[1, var.name] <- frontier.projections[, paste0("GDP", year)] *
      tech.USA.projections[, paste0("Tech", year)] *
      pop.projections.USA[, paste0("Pop", year)] / 10
  }
  
  # Fuck R's conversion away from numeric type
  emissions.projections[,-1] <- sapply(emissions.projections[,-1], as.numeric)
  emissions.projections
}


get.ipat.pastdata <- function(data.medium, year.sequence, names.countries,
                              na.remove=F) {
  # get.ipat.pastdata gets IPAT variables worldwide for the past data.
  # If na.remove is FALSE, there will be missing data for GDP prior to
  # 1990 (before the end of the Soviet Union) and in 2010 (some GDP data
  # missing), and Tech will also be missing for these countries.
  # If na.remove is TRUE, there will be full observations for GDP and Tech,
  # though this will be misleading, as GDP numbers will be too low for
  # some years, which can also bias intensity upwards.
  # The results of get.ipat.pastdata are used for plotting,
  # and not model fitting.
  library(doBy)
  
  data.tmp1 <- subset(data.medium, Year %% 5 == 0)
  data.tmp <- subset(data.tmp1, Isocode %in% names.countries)
  data.tmp$Pop <- data.tmp$PopTotal * 1000
  pop.data <- summaryBy(Pop ~ Year, data=data.tmp, FUN=sum)
  data.tmp$GDPTotal <- data.tmp$GDP * data.tmp$Pop
  gdp.data <- summaryBy(GDPTotal ~ Year, data=data.tmp, FUN=sum, na.rm=na.remove)
  data.tmp$CO2Total <- data.tmp$CO2 * data.tmp$Pop
  co2.data <- summaryBy(CO2Total ~ Year, data=data.tmp, FUN=sum, na.rm=na.remove)
  gdp.percapita <- gdp.data$GDPTotal.sum / pop.data$Pop.sum
  co2.percapita <- co2.data$CO2Total.sum / pop.data$Pop.sum
  tech.data <- co2.percapita / gdp.percapita * 10^4
  
  data.frame(Year=year.sequence, Pop=pop.data$Pop.sum, GDP=gdp.data$GDPTotal.sum,
             GDPpercapita=gdp.percapita,
             CO2=co2.data$CO2Total.sum, CO2percapita=co2.percapita,
             Tech=tech.data)
}



get.ipat.components.total <- function(tech.projections, gdp.projections, 
                                      frontier.projections, tech.USA.projections,
                                      pop.projections, co2.projections,
                                      year.sequence) {
  # get.ipat.components.total gets IPAT components for the world for our
  # projections.
  
  # US is special
  pop.USA.projections <- pop.projections[1,]
  stopifnot(pop.USA.projections$Isocode == "USA")
  pop.projections.mod <- pop.projections[-1,]
  
  names.tmp <- c("Year", "Pop", "GDP", "GDPpercapita", "CO2", "CO2percapita", "Tech")
  data.world <- data.frame(matrix(nrow=length(year.sequence), ncol=length(names.tmp)))
  names(data.world) <- names.tmp
  data.world[, "Year"] <- year.sequence
  index <- 0
  for (year in year.sequence) {
    index <- index + 1
    var.name.co2 <- paste0("CO2", year)
    var.name.gdp <- paste0("GDP", year)
    var.name.pop <- paste0("Pop", year)
    # Sum over population
    data.world[index, "Pop"] <- sum(pop.projections[, var.name.pop]) * 1000
    # Sum over GDP
    stopifnot(all(gdp.projections$Isocode == pop.projections.mod$Isocode))
    gdp.sum.tmp1 <- sum(gdp.projections[, var.name.gdp] * 
                          pop.projections.mod[, var.name.pop]) * 1000
    gdp.sum.tmp2 <- frontier.projections[, var.name.gdp] * 
      pop.USA.projections[, var.name.pop] * 1000
    data.world[index, "GDP"] <- gdp.sum.tmp1 + gdp.sum.tmp2
    # Sum over emisssions
    data.world[index, "CO2"] <- sum(co2.projections[, var.name.co2])
  }
  # Get GDP per capita
  data.world$GDPpercapita <- data.world$GDP / data.world$Pop
  # Get carbon per capita. Note: this is what is CO2 in data.medium.
  data.world$CO2percapita <- data.world$CO2 / data.world$Pop
  # Get tech
  data.world$Tech <- data.world$CO2percapita / data.world$GDPpercapita * 10^4
  
  data.world
}


get.ipat.components.bycountry <- function(tech.projections, gdp.projections, 
                                          frontier.projections, tech.USA.projections,
                                          pop.projections, co2.projections,
                                          year.sequence) {
  # get.ipat.components.bycountry reformats projections of IPAT components
  # at a country level.
  
  # US is special
  pop.USA.projections <- pop.projections[1,]
  stopifnot(pop.USA.projections$Isocode == "USA")
  pop.projections.mod <- pop.projections[-1,]
  stopifnot(all(pop.projections.mod$Isocode == gdp.projections$Isocode))
  
  #  names.tmp <- c("Year", "Pop", "GDP", "GDPpercapita", "CO2", "CO2percapita", "Tech")
  names.tmp <- c("Year", "Pop", "GDP", "GDPpercapita", "FrontierGap", "CO2", "CO2percapita", "Tech")
  data.ipat.bycountry <- list()
  # US first:
  data.ipat.bycountry[["USA"]] <- data.frame(matrix(nrow=length(year.sequence), ncol=length(names.tmp)))
  names(data.ipat.bycountry[["USA"]]) <- names.tmp
  data.ipat.bycountry[["USA"]]$Year <- year.sequence
  rowindex <- 0
  co2.rowindex.USA <- which(co2.projections$Isocode == "USA")
  tech.rowindex.USA <- which(tech.projections$Isocode == "USA")
  for (year in year.sequence) {
    rowindex <- rowindex + 1
    var.name.pop <- paste0("Pop", year)
    var.name.gdp <- paste0("GDP", year)
    var.name.co2 <- paste0("CO2", year)
    var.name.tech <- paste0("Tech", year)
    data.ipat.bycountry[["USA"]][rowindex, "Pop"] <- pop.USA.projections[var.name.pop]
    data.ipat.bycountry[["USA"]][rowindex, "GDPpercapita"] <- frontier.projections[var.name.gdp]
    data.ipat.bycountry[["USA"]][rowindex, "CO2"] <- co2.projections[co2.rowindex.USA, var.name.co2]
    data.ipat.bycountry[["USA"]][rowindex, "Tech"] <- tech.USA.projections[var.name.tech]
  }
  data.ipat.bycountry[["USA"]][, "FrontierGap"] <- rep(0, length(year.sequence))
  data.ipat.bycountry[["USA"]]$CO2percapita <- data.ipat.bycountry[["USA"]]$CO2 /
    (10^3 * data.ipat.bycountry[["USA"]]$Pop)
  data.ipat.bycountry[["USA"]]$GDP <- data.ipat.bycountry[["USA"]]$GDPpercapita *
    (10^3 * data.ipat.bycountry[["USA"]]$Pop)
  
  for (iso in gdp.projections$Isocode) {
    data.ipat.bycountry[[iso]] <- data.frame(matrix(nrow=length(year.sequence), ncol=length(names.tmp)))
    names(data.ipat.bycountry[[iso]]) <- names.tmp
    data.ipat.bycountry[[iso]]$Year <- year.sequence
    pop.rowindex <- which(pop.projections$Isocode == iso)
    gdp.rowindex <- which(gdp.projections$Isocode == iso)
    co2.rowindex <- which(co2.projections$Isocode == iso)
    tech.rowindex <- which(tech.projections$Isocode == iso)
    rowindex <- 0
    for (year in year.sequence) {
      rowindex <- rowindex + 1
      var.name.pop <- paste0("Pop", year)
      var.name.gdp <- paste0("GDP", year)
      var.name.co2 <- paste0("CO2", year)
      var.name.tech <- paste0("Tech", year)
      data.ipat.bycountry[[iso]][rowindex, "Pop"] <- pop.projections[pop.rowindex, var.name.pop]
      # Confusingly, for hte projections GDP is GDP per capita while CO2 is total co2 emissions
      data.ipat.bycountry[[iso]][rowindex, "GDPpercapita"] <- gdp.projections[gdp.rowindex, var.name.gdp]
      data.ipat.bycountry[[iso]][rowindex, "CO2"] <- co2.projections[co2.rowindex, var.name.co2]
      data.ipat.bycountry[[iso]][rowindex, "Tech"] <- tech.projections[tech.rowindex, var.name.tech]
    }
    data.ipat.bycountry[[iso]]$FrontierGap <- log(data.ipat.bycountry[["USA"]]$GDPpercapita /
                                                    data.ipat.bycountry[[iso]]$GDPpercapita)
    data.ipat.bycountry[[iso]]$CO2percapita <- data.ipat.bycountry[[iso]]$CO2 /
      (10^3 * data.ipat.bycountry[[iso]]$Pop)
    data.ipat.bycountry[[iso]]$GDP <- data.ipat.bycountry[[iso]]$GDPpercapita *
      (10^3 * data.ipat.bycountry[[iso]]$Pop)
  } 
  
  data.ipat.bycountry
}


get.ipat.components.ssa <- function(tech.projections, gdp.projections, 
                                    pop.projections, co2.projections,
                                    year.sequence) {
  # get.ipat.components.total gets IPAT components for Sub-Saharan Africa
  # for our projections. Similar to get.ipat.components.total
  
  # Subset to Sub-Saharan Africa (SSA)
  tech.projections <- subset(tech.projections, Isocode %in% ssa.isos)
  gdp.projections <- subset(gdp.projections, Isocode %in% ssa.isos)
  pop.projections <- subset(pop.projections, Isocode %in% ssa.isos)
  co2.projections <- subset(co2.projections, Isocode %in% ssa.isos)
  
  names.tmp <- c("Year", "Pop", "GDP", "GDPpercapita", "CO2", "CO2percapita", "Tech")
  data.ssa <- data.frame(matrix(nrow=length(year.sequence), ncol=length(names.tmp)))
  names(data.ssa) <- names.tmp
  data.ssa[, "Year"] <- year.sequence
  index <- 0
  for (year in year.sequence) {
    index <- index + 1
    var.name.co2 <- paste0("CO2", year)
    var.name.gdp <- paste0("GDP", year)
    var.name.pop <- paste0("Pop", year)
    # Sum over population
    data.ssa[index, "Pop"] <- sum(pop.projections[, var.name.pop]) * 1000
    # Sum over GDP
    stopifnot(all(gdp.projections$Isocode == pop.projections$Isocode))
    data.ssa[index, "GDP"] <- sum(gdp.projections[, var.name.gdp] * 
                                    pop.projections[, var.name.pop]) * 1000
    # Sum over emisssions
    data.ssa[index, "CO2"] <- sum(co2.projections[, var.name.co2])
  }
  # Get GDP per capita
  data.ssa$GDPpercapita <- data.ssa$GDP / data.ssa$Pop
  # Get carbon per capita. Note: this is what is CO2 in data.medium.
  data.ssa$CO2percapita <- data.ssa$CO2 / data.ssa$Pop
  # Get tech
  data.ssa$Tech <- data.ssa$CO2percapita / data.ssa$GDPpercapita * 10^4
  
  data.ssa
}




get.trajs.quants <- function(co2.projections, year.seq, quantiles=c(0.1,0.5,0.9)) {
  # get.trajs.quants gets quantiles for worldwide CO2 emissions
  trajs.quants <- data.frame(Quantile=quantiles)
  for (year in year.seq) {
    var.name <- paste0("CO2", year)
    trajs.quants <- cbind(trajs.quants,
                          quantile(sapply(1:length(co2.projections), function(i)
                            #  median(sapply(1:length(co2.projections), function(i)
                            sum(co2.projections[[i]][, var.name])), probs=quantiles) / 10^9)
  }
  names(trajs.quants)[-1] <- paste0("CO2", year.seq)
  
  trajs.quants
}



get.trajs.quants.bycountry <- function(co2.projections, year.seq, 
                                       quantiles=c(0.1,0.5,0.9)) {
  # get.trajs.quants.bycountry gets quantiles for CO2 emissions by country.
  # Compare with get.trajs.quants
  isos <- co2.projections[[1]]$Isocode
  trajs.quants <- list()
  for (year in year.seq) {
    trajs.quants[[paste0("CO2", year)]] <- data.frame(Isocode=isos)
    for (q in quantiles) {
      trajs.quants[[paste0("CO2", year)]][, paste0("CO2_", q)] <- rep(NA, length(isos))
    }
  }
  for (indx in 1:length(isos)) {
    iso <- isos[indx]
    for (year in year.seq) {
      var.name <- paste0("CO2", year)
      trajs.quants[[var.name]][indx, -1] <-  t(quantile(sapply(1:length(co2.projections), function(i)
        co2.projections[[i]][indx, var.name]), probs=quantiles, na.rm=T) / 10^9)
    }
  }
  
  trajs.quants
}


get.ipat.medians <- function(trajs.worldwide, quantiles=c(0.1, 0.5, 0.9)) {
  # get.ipat.medians gets worldwide quantiles (not just medians) of IPAT variables,
  # based on projections trajs.worldwide
  trajs.quants <- list()
  for (q in quantiles) {
    trajs.quants[[as.character(q)]] <- trajs.worldwide[[1]]
    trajs.quants[[as.character(q)]][, -1] <- NA
  }
  
  row.index <- 0
  for (year in trajs.quants[[1]]$Year) {
    row.index <- row.index + 1
    for (var.name in names(trajs.worldwide[[1]])[-1]) {
      if (trajs.quants[[1]]$Year[1] >= 1990) {
        quants.tmp <- quantile(sapply(1:length(trajs.worldwide), function(i)
          trajs.worldwide[[i]][row.index, var.name]), probs=quantiles)
      } else {
        quants.tmp <- quantile(sapply(1:length(trajs.worldwide), function(i)
          trajs.worldwide[[i]][row.index, var.name]), probs=quantiles, na.rm=T)
      }
      for (q_ind in 1:length(quantiles)) {
        q <- quantiles[q_ind]
        trajs.quants[[as.character(q)]][row.index, var.name] <- quants.tmp[q_ind]
      }
    }
  }
  
  trajs.quants
}


get.ipat.quantiles.bycountry <- function(ipat.components.bycountry, quantiles=c(0.1, 0.5, 0.9),
                                         na.rm=FALSE) {
  # get.ipat.quantiles.bycountry gives IPAT quantiles on a country level.
  # It's a variation on get.ipat.medians, despite the name change from "medians" to "quantiles"
  
  ipat.quants.bycountry <- list()
  for (iso in names(ipat.components.bycountry[[1]])) {
    ipat.quants.bycountry[[iso]] <- list()
    for (q in quantiles) {
      ipat.quants.bycountry[[iso]][[as.character(q)]] <- ipat.components.bycountry[[1]][[iso]]
      ipat.quants.bycountry[[iso]][[as.character(q)]][, -1] <- NA
    }
    
    row.index <- 0
    for (year in ipat.components.bycountry[[1]][[iso]]$Year) {
      row.index <- row.index + 1
      for (var.name in names(ipat.components.bycountry[[1]][[iso]])[-1]) {
        if (ipat.components.bycountry[[1]][[iso]]$Year[1] >= 1990) {
          quants.tmp <- quantile(sapply(1:length(ipat.components.bycountry), function(i)
            ipat.components.bycountry[[i]][[iso]][row.index, var.name]), probs=quantiles,
            na.rm=na.rm)
        } else {
          # Need to remove missing because of missing data for some countries
          quants.tmp <- quantile(sapply(1:length(ipat.components.bycountry), function(i)
            ipat.components.bycountry[[i]][[iso]][row.index, var.name]), probs=quantiles,
            na.rm=T)
        }
        for (q_ind in 1:length(quantiles)) {
          q <- quantiles[q_ind]
          ipat.quants.bycountry[[iso]][[as.character(q)]][row.index, var.name] <- quants.tmp[q_ind]
        }
      }
    }
  }
  
  ipat.quants.bycountry
}


# =========================================



# Modelling Carbon Intensity:
make.bugs.data.ar1 <- function(data.medium, max.vals, remove.before.peak=TRUE,
                               isos.remove=c()) {
  # make.bugs.data.ar1 creates the variables that will be used by JAGS in fitting
  # the GDP and intensity models given in corr_model_ar1trend_const.bug
  
  # Only removing indeces for the Intensity data, not for GDP data
  indeces.remove <- data.medium$Isocode %in% isos.remove
  indeces.remove.gdp <- data.medium$Isocode == "USA"
  data.tmp <- data.medium[!indeces.remove, c("Year","Isocode","Tech")]
  data.tmp$Tech <- log(data.tmp$Tech)
  data.wide.tech <- dcast(data.tmp[!duplicated(data.tmp[,c("Year","Isocode")]),],
                          Year ~ Isocode, value.var="Tech")
  data.tmp.gdp <- data.medium[!indeces.remove.gdp, c("Year","Isocode","GDP")]
  data.tmp.gdp$GDP <- log(data.tmp.gdp$GDP)
  data.wide.gdp <- dcast(data.tmp.gdp[!duplicated(data.tmp.gdp[,c("Year","Isocode")]),],
                         Year ~ Isocode, value.var="GDP")
  
  #browser()
  tau <- data.wide.tech
  n.years <- dim(tau)[1]
  n.countries <- dim(tau)[2]
  n.countries.gdp <- dim(data.wide.gdp)[2]
  
  # Get isos for GDP and intensity. 1st element is "Year", which we remove
  isos.gdp <- names(data.wide.gdp)[-1]
  isos.tau <- names(tau)[-1]
  
  # Make the frontier:
  frontier.gdp <- log(data.medium$GDP[data.medium$Isocode == "USA"])
  frontier.gap <- frontier.gdp - data.wide.gdp # first column is year and is meaningless
  
  frontier.gap.old.tmp <- frontier.gap[1:(n.years-1), 2:n.countries.gdp]
  frontier.gap.new.tmp <- frontier.gap[2:n.years, 2:n.countries.gdp]
  tau.delta.tmp.gdp <- frontier.gap.new.tmp - frontier.gap.old.tmp
  tau.old.matrix <- tau[1:(n.years-1), 2:n.countries]
  tau.new.matrix <- tau[2:n.years, 2:n.countries]
  n.years <- n.years - 1
  n.countries <- n.countries - 1
  n.countries.gdp <- n.countries.gdp - 1
  
  # Adding in the USA separately for intensity:
  us.tech <- log(data.medium$Tech[data.medium$Isocode == "USA"])
  tau.new.USA <- us.tech[-1]
  tau.old.USA <- us.tech[-length(us.tech)]
  
  tau.delta.tmp <- tau.new.matrix - tau.old.matrix
  # This block removes years before the peak
  if (remove.before.peak) {
    print("Removing data before the peak...")
    c <- 0
    # Assume they peaked at the start by default for now.
    # Using max.vals from find_max.R
    countries.peak <- rep(0, n.countries)
    for (iso.name in names(tau.delta.tmp)) {
      c <- c + 1
      if (sum(max.vals$Iso == iso.name) == 1) {
        # we have a max value for that country
        countries.peak[c] <- max.vals$Year[which(max.vals$Iso == iso.name)]
      }
    }
    smallest.year <- min(data.medium$Year)
    for (c in 1:n.countries) {
      for (t in 1:n.years) {
        if (t <= (countries.peak[c] - smallest.year)) {
          tau.delta.tmp[t,c] <- NA
        }
        # if (t == (countries.peak[c] - smallest.year) && !is.na(tau.delta.tmp[t,c]))
        # {
        #   if (tau.delta.tmp[t,c] < 0) {tau.delta.tmp[t,c] <- NA}
        # }
      }
    }
  } else {
    print("Not removing data before the peak.")
  }
  
  
  
  n.GoodObs <- sum(!is.na(tau.delta.tmp.gdp))
  n.GoodObs.c <- rep(0, n.countries.gdp)  
  frontier.gap.old <- frontier.gap.new <- rep(NA, n.GoodObs)
  tau.new <- tau.old <- years.intensity <- year.inds <- country.inds <- country.gdp.inds <- rep(NA, n.GoodObs)
  tau.new.funny.index <- c()
  intensity.Obs <- c()
  z <- 0
  z.intensity <- 0
  for (c in 1:n.countries.gdp){
    for (t in 1:n.years){
      if (!is.na(tau.delta.tmp.gdp[t,c])) {
        z <- z + 1
        n.GoodObs.c[c] <- n.GoodObs.c[c] + 1
        if (isos.gdp[c] %in% isos.tau) {
          c.tau <- which(isos.tau == isos.gdp[c])
          stopifnot(length(c.tau) == 1)
          if ((!is.na(tau.delta.tmp[t,c.tau])) & (!is.na(tau.new.matrix[t, c.tau]))
              & (!is.na(tau.old.matrix[t, c.tau]))) {
            z.intensity <- z.intensity + 1
            tau.new[z] <- tau.new.matrix[t, c.tau]
            tau.new.funny.index <- c(tau.new.funny.index, tau.new.matrix[t, c.tau])
            tau.old[z] <- tau.old.matrix[t, c.tau]
            years.intensity[z] <- t + 1960 - 1
            intensity.Obs <- c(intensity.Obs, z)
          }
          country.inds[z] <- c.tau
        }
        frontier.gap.old[z] <- frontier.gap.old.tmp[t, c]
        frontier.gap.new[z] <- frontier.gap.new.tmp[t, c]
        
        year.inds[z] <- t # not used now
        country.gdp.inds[z] <- c
      }
    }
    # cat(c, ' : ', n.GoodObs.c[c], '\n')
  }
  n.GoodObsIntensity <- z.intensity
  # browser()
  list(
    tau.new=tau.new, tau.old=tau.old,
    tau.new.funny.index=tau.new.funny.index,
    years.intensity=years.intensity,
    frontier.gap.old=frontier.gap.old, 
    frontier.gap.new=frontier.gap.new, frontier.gdp=frontier.gdp, 
    year.inds=year.inds,
    country.gdp.inds=country.gdp.inds,
    country.inds=country.inds,
    n.countries=n.countries, n.countries.gdp=n.countries.gdp, 
    n.years=n.years, n.GoodObs=n.GoodObs,
    names.countries=isos.tau,
    names.countries.gdp=isos.gdp,
    n.GoodObs.c=n.GoodObs.c,
    tau.new.USA=tau.new.USA,
    tau.old.USA=tau.old.USA,
    n.GoodObsIntensity=n.GoodObsIntensity,
    intensity.Obs=intensity.Obs)
}


fit.corr.model.ar1 <- function(data.medium, isos.remove=c("USA"), max.vals,
                               n.iterations=1000, n.adapt=300, n.chains=5, thin=20,
                               model.name="corr_model_ar1trend.bug",
                               var.list=c('delta', 'sig.eps', 'mu.delta', 'sig.delta',
                                          'phi', 'sig.gap', 'gamma', 'gamma.pre1973',
                                          'sig.frontier', 'mu.phi', 'sig.phi',
                                          'mu.ln.gap', 'sig.ln.gap', 'sig.eps.mean','sig.eps.sd',
                                          'rho',
                                          'eta', 'mu.eta', 'sig.eta',
                                          'beta', 'mu.beta', 'sig.beta',
                                          'delta.USA', 'sig.eps.USA')) {
  # fit.corr.model.ar1 uses JAGS to fit the model described in
  # corr_model_ar1trend_const.bug for GDP and intensity.
  # It calls make.bugs.data.ar1() to get the variables
  # that will be used by JAGS, when calling jags.model() and coda.samples()
  
  library(rjags)
  library(coda)
  library(reshape2)
  print("Calling make.bugs.data")
  jags.input <- make.bugs.data.ar1(data.medium, max.vals,
                                   isos.remove=isos.remove, remove.before.peak=T)
  # browser()
  print("Successfully called make.bugs.data")
  tau.new <- jags.input$tau.new
  tau.new.funny.index <- jags.input$tau.new.funny.index
  tau.old <- jags.input$tau.old
  years.intensity <- jags.input$years.intensity
  frontier.gap.old <- jags.input$frontier.gap.old
  frontier.gap <- jags.input$frontier.gap.new
  frontier.gdp <- jags.input$frontier.gdp
  year.inds <- jags.input$year.inds
  
  country.gdp.inds <- jags.input$country.gdp.inds # todo
  country.inds <- jags.input$country.inds
  n.countries <- jags.input$n.countries
  n.countries.gdp <- jags.input$n.countries.gdp
  n.years <- jags.input$n.years
  n.GoodObs <- jags.input$n.GoodObs
  names.countries <- jags.input$names.countries
  names.countries.gdp <- jags.input$names.countries.gdp
  n.GoodObs.c <- jags.input$n.GoodObs.c
  tau.new.USA <- jags.input$tau.new.USA
  tau.old.USA <- jags.input$tau.old.USA
  # Treat GDP and Intensity separately
  n.GoodObsIntensity <- jags.input$n.GoodObsIntensity
  intensity.Obs <- jags.input$intensity.Obs
  
  frontier.change <- frontier.gdp[2:(n.years+1)] - frontier.gdp[1:n.years]
  jags.input$frontier.change <- frontier.change
  
  print("here")
  print(n.countries)
  print("here")
  #jags.out <- jags.model(model.name, n.chains=n.chains, n.adapt=n.adapt)
  jags.out <- jags.model(model.name, 
                         data=jags.input,
                         n.chains=n.chains, n.adapt=n.adapt)
  jags.output <- coda.samples(jags.out,
                              var.list,
                              n.iter=n.iterations,
                              thin=thin)
  list(jags.input, jags.output)
}




project.corr.ar1.const <- function(data.medium, input.data, jags.corr.output, year.range, 
                                   param.indeces=seq(400, 500, by=10)) {
  # project.corr.ar1.const uses the fitted model for GDP and intensity
  # as fit by fit.corr.model.ar1 and generations projections forward
  
  jags.input <- input.data
  
  tau.change <- jags.input$tau.delta
  frontier.gap.old <- jags.input$frontier.gap.old
  frontier.gap <- jags.input$frontier.gap.new
  frontier.gdp <- jags.input$frontier.gdp
  year.inds <- jags.input$year.inds
  country.gdp.inds <- jags.input$country.gdp.inds
  country.inds <- jags.input$country.inds
  n.countries <- jags.input$n.countries
  n.countries.gdp <- jags.input$n.countries.gdp
  n.years <- jags.input$n.years
  n.GoodObs <- jags.input$n.GoodObs
  names.countries <- jags.input$names.countries
  names.countries.gdp <- jags.input$names.countries.gdp
  n.GoodObs.c <- jags.input$n.GoodObs.c
  tau.change.USA <- jags.input$tau.change.USA
  # Treat GDP and Intensity separately
  n.GoodObsIntensity <- jags.input$n.GoodObsIntensity
  intensity.Obs <- jags.input$intensity.Obs
  
  #CHANGE
  #year.first <- year.range[1]
  year.first <- min(year.range[1], 2008)
  year.last <- year.range[2]
  print("Time to create the base datasets")
  # Create base GDP data:
  gdp.data.tmp1 <- data.medium[data.medium$Year == year.first, c("Isocode", "GDP")]
  gdp.data.tmp2 <- gdp.data.tmp1[gdp.data.tmp1$Isocode %in% names.countries.gdp,
                                 c("Isocode", "GDP")]
  gdp.data.tmp <- cbind(gdp.data.tmp2, t(rep(NA, (year.last - year.first))))
  names(gdp.data.tmp) <- c("Isocode", paste0("LogGDP", year.first:year.last))
  gdp.data.tmp[, 2] <- log(gdp.data.tmp[, 2])
  
  # Create base GDP frontier data:
  stopifnot(all(gdp.data.tmp$Isocode == names.countries.gdp))
  stopifnot(! "USA" %in% names.countries.gdp)
  gdp.frontier <- data.frame(Isocode="USA", 
                             GDP=data.medium$GDP[data.medium$Isocode == "USA" &
                                                   data.medium$Year == year.first],
                             t(rep(NA, (year.last - year.first))))
  names(gdp.frontier) <- c("Isocode", paste0("LogGDP", year.first:year.last))
  gdp.frontier[, 2] <- log(gdp.frontier[, 2])
  # Create base Intensity data:
  tech.data.tmp1 <- data.medium[data.medium$Year == year.first, c("Isocode", "Tech")]
  # Use names.countries.gdp b/c we'll project ahead for all countries.
  tech.data.tmp2 <- tech.data.tmp1[tech.data.tmp1$Isocode %in% names.countries.gdp,
                                   c("Isocode", "Tech")]
  tech.data.tmp <- cbind(tech.data.tmp2, t(rep(NA, (year.last - year.first))))
  names(tech.data.tmp) <- c("Isocode", paste0("LogTech", year.first:year.last))
  tech.data.tmp[, 2] <- log(tech.data.tmp[, 2])
  # Create base USA Intensity data
  tech.USA.tmp <- data.frame(Isocode="USA", 
                             Tech=data.medium$Tech[data.medium$Isocode == "USA" &
                                                     data.medium$Year == year.first],
                             t(rep(NA, (year.last - year.first))))
  names(tech.USA.tmp) <- c("Isocode", paste0("LogTech", year.first:year.last))
  tech.USA.tmp[, 2] <- log(tech.USA.tmp[, 2])
  
  # Project ahead for most countries
  print("Projecting ahead for countries")
  names.excluded <- setdiff(names.countries.gdp,
                            c("USA", names.countries))
  n.excluded <- length(names.excluded)
  n.chains <- length(jags.corr.output)
  list.proj <- list()
  list.index <- 0
  for (chain.num in 1:n.chains) {
    # browser()
    for (index in param.indeces) {
      list.index <- list.index + 1
      nrow.corr <- nrow(jags.corr.output[[chain.num]])
      while(1)
      {
        vars.tmp1 <- jags.corr.output[[chain.num]][index,]
        vars.tmp <- data.frame(names(vars.tmp1), vars.tmp1)
        
        eta.tmp <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 6) == "mu.eta"]
        if (eta.tmp < 0)
        {
          break
        }
        else
        {
          index <- index + 1
          if (index > nrow.corr)
          {
            index <- 1
          }
        }
      }
      
      
      vars.tmp <- data.frame(names(vars.tmp1), vars.tmp1)
      eta.tmp <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 6) == "mu.eta"]
      beta.tmp <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 7) == "mu.beta"]
      deltas.tmp <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 6) == "delta["]
      delta.USA <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 9) == "delta.USA"]
      phis <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 3) == "phi"]
      rho <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 3) == "rho"]
      sig.eps.vals.tmp <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 8) == "sig.eps["]
      sig.eps.USA <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 11) == "sig.eps.USA"]
      sig.gap <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 7) == "sig.gap"]
      sig.frontier <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 12) == "sig.frontier"]
      gamma <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 5) == "gamma"]
      stopifnot(length(gamma) == 1) # make sure not using gamma.pre1973
      
      # Used for prepeak countries. Pull out hyperparameters, then generate country parameters.
      deltas <- rep(NA, n.countries.gdp)
      sig.eps.vals <- rep(NA, n.countries.gdp)
      stopifnot(length(deltas) == length(phis))
      # Fill out deltas with postpeak countries
      gdp.intensity.inds <- names.countries.gdp %in% names.countries
      gdp.excluded.inds <- ! names.countries.gdp %in% names.countries
      stopifnot(names.countries.gdp[gdp.intensity.inds] == names.countries)
      deltas[gdp.intensity.inds] <- deltas.tmp
      sig.eps.vals[gdp.intensity.inds] <- sig.eps.vals.tmp
      # Get hyperdistributions for prepeak countries
      var.names <- varnames(jags.corr.output[[chain.num]])
      mu.delta.tmp <- vars.tmp1[which(var.names == "mu.delta")]
      sig.delta.tmp <- vars.tmp1[which(var.names == "sig.delta")]
      sig.eps.mean.tmp <- vars.tmp1[which(var.names == "sig.eps.mean")]
      sig.eps.sd.tmp <- vars.tmp1[which(var.names == "sig.eps.sd")]
      # Now generate parameters for prepeak countries
      delta.c <- rep(NA, n.excluded)
      sigma.c <- rep(NA, n.excluded)
      for (c.ind in 1:n.excluded) {
        delta.c[c.ind] <- max(0, rnorm(1, mu.delta.tmp, sig.delta.tmp))
        sigma.c[c.ind] <- rlnorm(1, sig.eps.mean.tmp, sig.eps.sd.tmp)
      }
      deltas[gdp.excluded.inds] <- delta.c
      sig.eps.vals[gdp.excluded.inds] <- sigma.c
      
      # Create dataframe for GDP. First year is last year in the data
      # Create dataframe for Intensity
      gdp.frontier.copy <- gdp.frontier
      gdp.data.copy <- gdp.data.tmp
      tech.data.copy <- tech.data.tmp
      tech.USA.copy <- tech.USA.tmp
      #CHANGE
      if (year.range[1] == 2010) {
        # We have missing data for some countries for 2009 and 2010, so we have
        # this very awkward kludge to deal with this inconsistency.
        # In this case we started with a base in 2008, so we now have to fill in
        # 2009 and 2010 data, for
        # gdp.frontier.copy, gdp.data.copy, tech.data.copy, tech.USA.copy
        for (year in 2009:2010) {
          print("Doing special 2009, 2010 stuff for USA")
          gdp.frontier.copy[, paste0("LogGDP", year)] <- log(subset(data.medium, Isocode == "USA" &
                                                                      Year == year)$GDP)
          tech.USA.copy[, paste0("LogTech", year)] <- log(subset(data.medium, Isocode == "USA" &
                                                                   Year == year)$Tech)
        }
        for (year in 2009:2010) {
          print("Doing special 2009, 2010 stuff for non-USA countries")
          # Project non-USA countries as if starting from 2008. We'll replace 
          # simulated values with true values we know after this for loop
          year.old <- year - 1
          gdp.name.now <- paste0("LogGDP", year)
          gdp.name.old <- paste0("LogGDP", year.old)
          tech.name.now <- paste0("LogTech", year)
          tech.name.old <- paste0("LogTech", year.old)
          counter <- 0
          while (counter == 0 ||
                 (max(tech.data.copy[, tech.name.now], na.rm=T) > log(50))) {
            counter <- 1
            # Project other countries GDP
            gdp.epsilon <- rnorm(n.countries.gdp, 0, sig.gap)
            gdp.data.copy[, gdp.name.now] <- gdp.frontier.copy[, gdp.name.now] - phis * 
              (gdp.frontier.copy[, gdp.name.old] - gdp.data.copy[, gdp.name.old]) -
              gdp.epsilon
            # Project Intensity
            tech.epsilon <- rnorm(n.countries.gdp,
                                  rho*(sig.eps.vals/sig.gap)*
                                    gdp.epsilon, 
                                  sqrt(1-rho)*sig.eps.vals)
            tech.data.copy[, tech.name.now] <- eta.tmp * (year - 1990) + 
              beta.tmp * tech.data.copy[, tech.name.old] -
              deltas + tech.epsilon
          }
        }
        for (year in 2009:2010) {
          # Now replace simulated values with true values we know for 2009 and 2010
          gdp.name.now <- paste0("LogGDP", year)
          tech.name.now <- paste0("LogTech", year)
          
          # GDP
          log.gdp.data.tmp <- log(subset(data.medium, Year == year & 
                                           Isocode %in% names.countries.gdp)[, "GDP"])
          which.have.gdp.data <- ! is.na(log.gdp.data.tmp)
          # Replace simulated values with observed values
          gdp.data.copy[which.have.gdp.data, gdp.name.now] <- log.gdp.data.tmp[which.have.gdp.data]
          
          # Intensity
          log.tech.data.tmp <- log(subset(data.medium, Year == year & 
                                            Isocode %in% names.countries.gdp)[, "Tech"])
          which.have.tech.data <- ! is.na(log.tech.data.tmp)
          # Replace simulated values with observed values
          tech.data.copy[which.have.tech.data, tech.name.now] <- log.tech.data.tmp[which.have.tech.data]        
        }
      }
      for (year in (year.range[1]+1):year.range[2]) {
        # year is the year we're projecting towards
        year.old <- year - 1
        gdp.name.now <- paste0("LogGDP", year)
        gdp.name.old <- paste0("LogGDP", year.old)
        tech.name.now <- paste0("LogTech", year)
        tech.name.old <- paste0("LogTech", year.old)
        counter <- 0
        while (counter == 0 ||
               (max(tech.data.copy[, tech.name.now], na.rm=T) > log(50))) {
          #counter <- 1
          counter <- counter + 1
          if (counter > 1) {
            tech.data.copy[, tech.name.now] <- sapply(tech.data.copy[, tech.name.now],
                                                      function(x) min(x, log(50))
            )
            break
          }
          stopifnot(counter < 10)
          # Project frontier GDP
          gdp.epsilon.frontier <- rnorm(1, 0, sig.frontier)
          gdp.frontier.copy[, gdp.name.now] <- gdp.frontier.copy[, gdp.name.old] + 
            gamma + gdp.epsilon.frontier
          # Project other countries GDP
          gdp.epsilon <- rnorm(n.countries.gdp, 0, sig.gap)
          gdp.data.copy[, gdp.name.now] <- gdp.frontier.copy[, gdp.name.now] - phis * 
            (gdp.frontier.copy[, gdp.name.old] - gdp.data.copy[, gdp.name.old]) -
            gdp.epsilon
          # Project Intensity
          tech.epsilon <- rnorm(n.countries.gdp,
                                rho*(sig.eps.vals/sig.gap)*
                                  gdp.epsilon, 
                                sqrt(1-rho)*sig.eps.vals)
          tech.data.copy[, tech.name.now] <- eta.tmp * (year - 1990) + 
            beta.tmp * tech.data.copy[, tech.name.old] -
            deltas + tech.epsilon
        }
        # Project USA Intensity
        tech.epsilon.USA <- rnorm(1,
                                  -rho*(sig.eps.USA/sig.frontier)*gdp.epsilon.frontier, 
                                  sqrt(1-rho)*sig.eps.USA)
        tech.USA.copy[, tech.name.now] <- eta.tmp * (year - 1990) +
          beta.tmp * tech.USA.copy[, tech.name.old] -
          delta.USA + tech.epsilon.USA
      }
      list.proj[[list.index]] <- list(GDPFrontier=gdp.frontier.copy, 
                                      GDPData=gdp.data.copy, TechData=tech.data.copy,
                                      TechUSA=tech.USA.copy)
    }
  }
  print("returning list.proj")
  
  list.proj
}





# =========================================

get.cumulative.emissions <- function(co2.projection) {
  # get.cumulative.emissions is for getting cumulative emissions from 2010 to 2100
  
  # Warning: only works if starting from 2010
  stopifnot(names(co2.projection)[1] == "Isocode")
  emit.yearly <- colSums(co2.projection[,-1])
  cum.emissions <- emit.yearly
  cum.emissions.total <- 0
  for (var.name in names(emit.yearly)) {
    if (! (var.name == "CO22010")) {
      cum.emissions.total <- cum.emissions.total + 2.5*emit.yearly[var.name]
    }
    cum.emissions[var.name] <- cum.emissions.total
    if (! (var.name == "CO22100")) {
      cum.emissions.total <- cum.emissions.total + 2.5*emit.yearly[var.name]
    }
  }
  cum.emissions.total <- unname(cum.emissions.total)
  
  return(list(cum.emissions.total=cum.emissions.total,
              cum.emissions=cum.emissions,
              emit.yearly=emit.yearly))
}


is.carbon.below.limit <- function(co2.projection,
                                  co2.cutoff=11000,
                                  year.list=seq(2015,2100,by=5)) {
  # is.carbon.below.limit checks whether a projection from 2010 to 2100
  # has cumulative CO2 emissions below the cutoff of 11000 gt CO2
  
  stopifnot(names(co2.projection)[1] == "Isocode")
  name.list <- paste0("CO2", year.list)
  carbon.cum.tmp <- 0
  # It's as if we're starting from 2015, as in the paper we cite talking about
  # remaining CO2
  for (name in name.list) {
    if (name == "CO22015" || name == "CO22100") {
      carbon.cum.tmp <- carbon.cum.tmp + 2.5*sum(co2.projection[, name])
    } else {
      carbon.cum.tmp <- carbon.cum.tmp + 5*sum(co2.projection[, name])
    }
  }
  carbon.cum <- carbon.cum.tmp / 10^9 # convert to gigatonnes
  stopifnot(!is.na(carbon.cum))
  print(carbon.cum)
  
  if (carbon.cum < co2.cutoff) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


fit.project.model <- function(data.medium, year.start=2010, year.end=2100,
                              fit.ar1.const=F, n.trajectories=1000) {
  # fit.project.model fits the models and generates forward projections.
  # It uses population predictions already calculated above, calls fit.corr.model.ar1
  # to fit the model for GDP and intensity, and calls project.corr.ar1.const
  # to get projections for GDP and intensity.
  # It then has to combine the population, GDP, and intensity projections.
  
  library(rjags)
  library(coda)
  library(reshape2)
  data.medium.tmp <- subset(data.medium, Year <= year.start)
  
  # Calling find_maxes() from find_max.R
  print("======================================")
  print("Find max intensities using find_max.R")
  print("")
  maxes.countries <- find_maxes(data.medium.tmp)
  max.vals <- maxes.countries[[1]]
  rejects.late <- maxes.countries[[3]]
  rejects.insuf <- maxes.countries[[4]]
  
  # These functions are from estimate_model_corr.R
  print("======================================")
  print("Fit the model! Using fit.corr.model")
  print("")
  n.chains <- 5
  n.iter <- 100000
  #n.iter <- 5000
  n.burnin <- 5000
  #n.burnin <- 1000
  thin <- 20
  #thin <- 1
  stopifnot(n.trajectories %% n.chains == 0)
  param.indeces <- round(seq(1, floor(n.iter / thin), length.out=(n.trajectories / n.chains)))
  #CHANGE
  #year.start.tmp <- min(2008, year.start)
  year.start.tmp <- year.start
  if (! fit.ar1.const) {
    corr.model.out <- fit.corr.model(data.medium.tmp, 
                                     isos.remove=c("USA",rejects.late, rejects.insuf), max.vals,
                                     n.iterations=n.iter, n.adapt=n.burnin, n.chains=5,
                                     thin=thin)
    names.countries <- corr.model.out[[1]]$names.countries.gdp # Now we're including all countries
    names.countries.intensity <- corr.model.out[[1]]$names.countries # For only post-peak countries
    
    print("======================================")
    print("Projecting intensity and GDP using project.corr")
    print("")
    list.proj <- project.corr(data.medium.tmp, input.data=corr.model.out[[1]], 
                              jags.corr.output=corr.model.out[[2]], year.range=c(year.start.tmp,year.end),
                              param.indeces=param.indeces)
  } else {
    # can't get gamma.pre1973, b/c when pattern matching later we'd pull it out along
    # with gamma
    #CHANGE
    #corr.model.out <- fit.corr.model.ar1(data.medium, 
    corr.model.out <- fit.corr.model.ar1(data.medium.tmp, 
                                         isos.remove=c("USA", rejects.late, rejects.insuf), max.vals,
                                         n.iterations=n.iter, n.adapt=n.burnin, n.chains=5,
                                         thin=thin,
                                         model.name="corr_model_ar1trend_const.bug",
                                         var.list=c('delta', 'sig.eps', 'mu.delta', 'sig.delta',
                                                    'phi', 'sig.gap', 'gamma', 
                                                    'sig.frontier', 'mu.phi', 'sig.phi',
                                                    'mu.ln.gap', 'sig.ln.gap', 'sig.eps.mean','sig.eps.sd',
                                                    'rho',
                                                    'mu.eta',
                                                    'mu.beta',
                                                    'delta.USA', 'sig.eps.USA',
                                                    'deviance'))
    names.countries <- corr.model.out[[1]]$names.countries.gdp # Now we're including all countries
    names.countries.intensity <- corr.model.out[[1]]$names.countries # For only post-peak countries
    
    print("======================================")
    print("Projecting intensity and GDP using project.corr")
    print("")
    list.proj <- project.corr.ar1.const(data.medium.tmp, input.data=corr.model.out[[1]], 
                                        jags.corr.output=corr.model.out[[2]], year.range=c(year.start.tmp,year.end),
                                        param.indeces=param.indeces)
  }
  
  # From predict_pop.R
  # Getting population trajectories are slow
  print("======================================")
  print("Fitting the population model and getting population trajectories!")
  print("")
  preds.countries.trajs <- 
    loadRData(paste0(data.location, "poppreds_formatted_", year.start, "_", year.end, ".rda"))
    # preds.countries.tmp <- predict.population(year.present=year.start, year.end=year.end, 
    #                                           make.new=FALSE)
    # #                                            n.iter=500, n.burnin=100)
    # print("Have grabbed population model fit.")
    # preds.countries <- preds.countries.tmp$preds.countries
    # pred.pop <- preds.countries.tmp$pred.pop
    # print(paste0("Converting format of population trajectories, selecting first ", 
    #              n.trajectories, " trajectories."))
    # preds.countries.trajs <- find.pop.trajectories(data.medium.tmp, pred.pop, 1:n.trajectories,
    #                                                c("USA", names.countries), year.start=year.start, year.end=year.end)
  # }
  
  # These functions are from estimate_model_corr.R
  print("======================================")
  print("Combining IPAT projections using list.proj and estimate.co2.projections")
  print("")
  data.proj <- lapply(list.proj, function(x) convert.projections(x$TechData, 
                                                                 x$GDPData, x$GDPFrontier, x$TechUSA))
  co2.projections.tmp <- lapply(1:length(data.proj), function(i) {
    x <- data.proj[[i]]
    preds.countries <- preds.countries.trajs[[i]]
    estimate.co2.projections(x$TechData, x$GDPData, x$GDPFrontier, x$TechUSA,
                             preds.countries, year.range=c(year.start, year.end))
  })
  
  if (year.end == 2100) {
    co2.projections.restrictions <- sapply(co2.projections.tmp, is.carbon.below.limit)
    co2.projections <- co2.projections.tmp[co2.projections.restrictions]
    data.proj.restrict <- data.proj[co2.projections.restrictions]
    preds.countries.trajs.restrict <- preds.countries.trajs[co2.projections.restrictions]
  } else {
    co2.projections.restrictions <- 1:n.trajectories
    co2.projections <- co2.projections.tmp
    data.proj.restrict <- data.proj
    preds.countries.trajs.restrict <- preds.countries.trajs
  }
  
  list(corr.model.out=corr.model.out, 
       data.proj=data.proj.restrict,
       preds.countries.trajs=preds.countries.trajs.restrict,
       co2.projections=co2.projections,
       co2.projections.restrictions=co2.projections.restrictions,
       names.countries.intensity=names.countries.intensity)
}


evaluate.projections <- function(data.medium.full, model.results, 
                                 year.start=2010, year.end=2100,
                                 quantiles=c(0.025, 0.1, 0.5, 0.9, 0.975),
                                 outofsample.validate=FALSE,
                                 restrict.isos=FALSE,
                                 isos.restricted=c()) {
  # evaluate.projections uses the results of fit.project.model(), and
  # creates various summaries of the model and its projections by calling functions
  # like get.ipat.components.total, get.ipat.components.bycountry.
  # If year.start is 1980, 1990, or 2000, we would typically also pass
  # outofsample.validate=TRUE, and checks of the out of sample validation are
  # performed.
  
  if (year.start == 2010) {
    data.proj <- model.results$data.proj
    preds.countries.trajs <- model.results$preds.countries.trajs
    co2.projections <- model.results$co2.projections
  } else {
    co2.proj.tmp <- model.results$co2.projections
    isos.include <- co2.proj.tmp[[1]][
      !is.na(co2.proj.tmp[[1]][, paste0("CO2", year.start)]), "Isocode"]
    data.proj <- lapply(model.results$data.proj, function(x) {
      lapply(x, function(y) {
        subset(y, Isocode %in% isos.include)
      })
    })
    preds.countries.trajs <- lapply(model.results$preds.countries.trajs,
                                    function(x)
                                      subset(x, Isocode %in% isos.include))
    co2.projections <- lapply(model.results$co2.projections, function(x)
      subset(x, Isocode %in% isos.include))
  }
  # Get worldwide data
  trajs.worldwide <- lapply(1:length(data.proj), function(i) {
    x <- data.proj[[i]]
    get.ipat.components.total(x$TechData, x$GDPData, x$GDPFrontier, x$TechUSA,
                              preds.countries.trajs[[i]], co2.projections[[i]], 
                              year.sequence=seq(year.start, year.end,5))
  })
  trajs.quants <- get.trajs.quants(co2.projections,
                                   seq(year.start, year.end, by=5),
                                   quantiles=quantiles)
  
  trajs.quants.bycountry <- get.trajs.quants.bycountry(co2.projections, 
                                                       seq(year.start, year.end, 5),
                                                       quantiles=quantiles)
  
  # Get data for Sub-Saharan Africa (SSA)
  trajs.ssa <- lapply(1:length(data.proj), function(i) {
    x <- data.proj[[i]]
    get.ipat.components.ssa(x$TechData, x$GDPData,
                            preds.countries.trajs[[i]], co2.projections[[i]], 
                            year.sequence=seq(year.start, year.end,5))
  })
  
  ipat.quantiles <- get.ipat.medians(trajs.worldwide, quantiles=quantiles)
  
  ipat.ssa.quantiles <- get.ipat.medians(trajs.ssa, quantiles=quantiles)
  
  if (! outofsample.validate) {
    stopifnot(year.start == 2010)
    # Worldwide cumulative emissions for each scenario,
    # and quantiles of cumulative emissions
    n.projections <- length(co2.projections)
    cum.emissions.tmp <- lapply(1:n.projections,
                                function(i) get.cumulative.emissions(co2.projections[[i]])$cum.emissions)
    cum.emissions <- data.frame(matrix(nrow=n.projections, ncol=length(cum.emissions.tmp[[1]])))
    cum.emissions.quantile <- data.frame(matrix(nrow=length(cum.emissions.tmp[[1]]), ncol=length(quantiles)))
    names(cum.emissions.quantile) <- paste0("Quant", quantiles)
    for (i in 1:n.projections) {
      cum.emissions[i,] <- cum.emissions.tmp[[i]]
    }
    for (year.tmp in 1:dim(cum.emissions)[2]) {
      cum.emissions.quantile[year.tmp,] <- quantile(cum.emissions[, year.tmp], probs=quantiles) / 10^9
    }
    cum.emissions.quantile$Year <- seq(2010,2100,5)
  }
  
  # Carbon emissions for Sub-Saharan Africa (SSA)
  if (year.start == 2010) {
    co2.ssa <- as.data.frame(matrix(nrow=length(co2.projections), 
                                    ncol=length(names(co2.projections[[1]])[-1]),
                                    dimnames=list(c(),names(co2.projections[[1]])[-1])))
    which.ssa.rows <- co2.projections[[1]]$Isocode %in% ssa.isos
    for (i in 1:length(co2.projections)) {
      co2.ssa[i,] <- colSums(co2.projections[[i]][which.ssa.rows, -1])
    }
    # Carbon quantiles for SSA
    co2.ssa.quants <- as.data.frame(
      matrix(nrow=length(quantiles),
             ncol=(1 + dim(co2.ssa)[2]),
             dimnames=list(c(), c("Quantiles", names(co2.ssa)))))
    co2.ssa.quants$Quantiles <- quantiles
    for (col.name in names(co2.ssa)) {
      co2.ssa.quants[, col.name] <- quantile(co2.ssa[, col.name], probs=quantiles) / 10^9
    }
    
    #GDP per capita for Sub-Saharan Africa (SSA)
    gdp.ssa <- pop.ssa <- as.data.frame(matrix(nrow=length(co2.projections), 
                                               ncol=length(seq(2010, 2100, 5)),
                                               dimnames=list(c(), paste0("GDP", seq(2010, 2100, 5)))))
    which.ssa.rows.gdp <- data.proj[[1]][["GDPData"]]$Isocode %in% ssa.isos
    which.ssa.rows.pop <- preds.countries.trajs[[1]]$Isocode %in% ssa.isos
    
    gdp.names.tmp <- paste0("GDP", seq(2010, 2100, 5))
    for (i in 1:length(co2.projections)) {
      gdp.ssa[i,] <- colSums(data.proj[[i]][["GDPData"]][which.ssa.rows.gdp, gdp.names.tmp] *
                               preds.countries.trajs[[i]][which.ssa.rows.pop, -1] * 1000)
    }
    
    # Carbon quantiles for SSA
    co2.ssa.quants <- as.data.frame(
      matrix(nrow=length(quantiles),
             ncol=(1 + dim(co2.ssa)[2]),
             dimnames=list(c(), c("Quantiles", names(co2.ssa)))))
    co2.ssa.quants$Quantiles <- quantiles
    for (col.name in names(co2.ssa)) {
      co2.ssa.quants[, col.name] <- quantile(co2.ssa[, col.name], probs=quantiles) / 10^9
    }
    
    # GDP per capita quantiles for SSA
    gdp.ssa.quants <- as.data.frame(
      matrix(nrow=length(quantiles),
             ncol=(1 + dim(gdp.ssa)[2]),
             dimnames=list(c(), c("Quantiles", names(gdp.ssa)))))
    gdp.ssa.quants$Quantiles <- quantiles
    for (col.name in names(gdp.ssa)) {
      gdp.ssa.quants[, col.name] <- quantile(gdp.ssa[, col.name], probs=quantiles)
    }
  } else {
    # Fill out dummy categories
    co2.ssa <- as.data.frame(matrix(nrow=length(co2.projections), 
                                    ncol=length(names(co2.projections[[1]])[-1]),
                                    dimnames=list(c(),names(co2.projections[[1]])[-1])))
    co2.ssa.quants <- as.data.frame(
      matrix(nrow=length(quantiles),
             ncol=(1 + dim(co2.ssa)[2]),
             dimnames=list(c(), c("Quantiles", names(co2.ssa)))))
    
    gdp.ssa <- as.data.frame(matrix(nrow=length(co2.projections), 
                                    ncol=length(seq(year.start, year.end, 5)),
                                    dimnames=list(c(), paste0("GDP", seq(year.start, year.end, 5)))))
    gdp.ssa.quants <- as.data.frame(
      matrix(nrow=length(quantiles),
             ncol=(1 + dim(gdp.ssa)[2]),
             dimnames=list(c(), c("Quantiles", names(gdp.ssa)))))
  }
  
  # Cumulative SSA emissions
  if (! outofsample.validate) {
    co2.cum.ssa <- co2.ssa # just to size the array
    for (i in 1:dim(co2.ssa)[1]) {
      co2.cum.ssa[i,] <- get.cumulative.emissions(cbind(Isocode="SSA", co2.ssa[i,]))$cum.emissions
    }
    
    cum.emissions.ssa.quantile <- data.frame(matrix(nrow=length(seq(2010,2100,5)), ncol=length(quantiles)))
    names(cum.emissions.ssa.quantile) <- paste0("Quant", quantiles)
    for (year.tmp in 1:dim(cum.emissions)[2]) {
      cum.emissions.ssa.quantile[year.tmp,] <- quantile(co2.cum.ssa[, year.tmp], probs=quantiles) / 10^9
    }
    cum.emissions.ssa.quantile$Year <- seq(2010,2100,5)
  }
  
  # Cumulative emissions, excluding SSA
  if (! outofsample.validate) {
    stopifnot(year.start == 2010)
    # Worldwide cumulative emissions for each scenario,
    # and quantiles of cumulative emissions
    n.projections <- length(co2.projections)
    row.inds.notssa <- which(! co2.projections[[1]]$Isocode %in% ssa.isos)
    cum.emissions.tmp.notssa <- lapply(1:n.projections,
                                       function(i)
                                         get.cumulative.emissions(co2.projections[[i]][row.inds.notssa,])$cum.emissions)
    cum.emissions.notssa <- data.frame(matrix(nrow=n.projections, ncol=length(cum.emissions.tmp.notssa[[1]])))
    cum.emissions.notssa.quantile <- data.frame(matrix(nrow=length(cum.emissions.tmp.notssa[[1]]), ncol=length(quantiles)))
    names(cum.emissions.notssa.quantile) <- paste0("Quant", quantiles)
    for (i in 1:n.projections) {
      cum.emissions.notssa[i,] <- cum.emissions.tmp.notssa[[i]]
    }
    for (year.tmp in 1:dim(cum.emissions)[2]) {
      cum.emissions.notssa.quantile[year.tmp,] <- quantile(cum.emissions.notssa[, year.tmp], probs=quantiles) / 10^9
    }
    cum.emissions.notssa.quantile$Year <- seq(2010,2100,5)
  }
  
  
  # GDP split by OECD/non-OECD countries. Find quantiles
  if (! outofsample.validate) {
    data.proj <- model.results$data.proj
    pop.preds <- model.results$preds.countries.trajs
    which.oecd.not.us <- data.proj[[1]][["GDPData"]][, "Isocode"] %in% oecd.15.iso
    which.world.notoecd <- data.proj[[1]][["GDPData"]][, "Isocode"] %in% oecd.world.minusoecd15
    n.projections <- length(data.proj)
    stopifnot(pop.preds[[1]][1, "Isocode"] == "USA")
    # bc GDP is projected yearly but Pop isn't
    gdp.columns <- paste0("GDP", seq(2010, 2100, by=5))
    gdp.by.oecd <- gdppercapita.by.oecd <- list()
    gdp.tot.names <- paste0("GDPTotal", seq(2010, 2100, by=5))
    gdp.by.oecd[["OECD"]] <- gdp.by.oecd[["NotOECD"]] <-
      gdp.by.oecd[["WorldNotOECD"]] <-
      gdppercapita.by.oecd[["OECD"]] <- gdppercapita.by.oecd[["NotOECD"]] <- as.data.frame(
        matrix(nrow=n.projections, ncol=length(gdp.columns)))
    names(gdp.by.oecd[["OECD"]]) <- names(gdp.by.oecd[["NotOECD"]]) <-
      names(gdp.by.oecd[["WorldNotOECD"]]) <- gdp.tot.names
    names(gdppercapita.by.oecd[["OECD"]]) <- names(gdppercapita.by.oecd[["NotOECD"]]) <- gdp.tot.names
    for (i in 1:n.projections) {
      popdata.nonus <- pop.preds[[i]][-1,]
      # -1 to remove Isocode column
      gdp.oecd <- colSums(data.proj[[i]][["GDPData"]][which.oecd.not.us, gdp.columns] *
                            popdata.nonus[which.oecd.not.us, -1] * 1000) +
        data.proj[[1]][["GDPFrontier"]][, gdp.columns] * pop.preds[[i]][1,-1] * 1000
      gdp.worldnotoecd <- colSums(data.proj[[i]][["GDPData"]][which.world.notoecd, gdp.columns] *
                                    popdata.nonus[which.world.notoecd, -1] * 1000)
      pop.oecd <- colSums(popdata.nonus[which.oecd.not.us, -1] * 1000) +
        pop.preds[[i]][1,-1] * 1000
      gdp.not.oecd <- colSums(data.proj[[i]][["GDPData"]][! which.oecd.not.us, gdp.columns] *
                                popdata.nonus[! which.oecd.not.us, -1] * 1000)
      pop.not.oecd <- colSums(popdata.nonus[! which.oecd.not.us, -1] * 1000)
      gdp.by.oecd[["OECD"]][i, ] <- gdp.oecd
      gdp.by.oecd[["WorldNotOECD"]][i, ] <- gdp.worldnotoecd
      gdp.by.oecd[["NotOECD"]][i, ] <- gdp.not.oecd
      gdppercapita.by.oecd[["OECD"]][i, ] <- gdp.oecd / pop.oecd
      gdppercapita.by.oecd[["NotOECD"]][i, ] <- gdp.not.oecd / pop.not.oecd
    }
    gdp.by.oecd.quant <- gdppercapita.by.oecd.quant <- list()
    gdp.by.oecd.quant[["OECD"]] <- gdp.by.oecd.quant[["NotOECD"]] <- 
      gdp.by.oecd.quant[["WorldNotOECD"]] <- as.data.frame(
        matrix(nrow=length(quantiles), ncol=(1 + length(gdp.columns))))
    gdppercapita.by.oecd.quant[["OECD"]] <- gdppercapita.by.oecd.quant[["NotOECD"]] <- as.data.frame(
      matrix(nrow=length(quantiles), ncol=(1 + length(gdp.columns))))
    names(gdp.by.oecd.quant[["OECD"]]) <- names(gdp.by.oecd.quant[["WorldNotOECD"]]) <-
      names(gdp.by.oecd.quant[["NotOECD"]]) <- c("Quantile", gdp.tot.names)
    names(gdppercapita.by.oecd.quant[["OECD"]]) <- names(gdppercapita.by.oecd.quant[["NotOECD"]]) <- c("Quantile", gdp.tot.names)
    gdp.by.oecd.quant[["OECD"]][, "Quantile"] <- 
      gdp.by.oecd.quant[["WorldNotOECD"]][, "Quantile"] <- 
      gdp.by.oecd.quant[["NotOECD"]][, "Quantile"] <- quantiles
    gdppercapita.by.oecd.quant[["OECD"]][, "Quantile"] <- 
      gdppercapita.by.oecd.quant[["NotOECD"]][, "Quantile"] <- quantiles
    for (var.name in gdp.tot.names) {
      gdp.by.oecd.quant[["OECD"]][, var.name] <- quantile(gdp.by.oecd[["OECD"]][, var.name], 
                                                          probs=quantiles)
      gdp.by.oecd.quant[["WorldNotOECD"]][, var.name] <- quantile(gdp.by.oecd[["WorldNotOECD"]][, var.name], 
                                                                  probs=quantiles)
      gdp.by.oecd.quant[["NotOECD"]][, var.name] <- quantile(gdp.by.oecd[["NotOECD"]][, var.name], 
                                                             probs=quantiles)
      gdppercapita.by.oecd.quant[["OECD"]][, var.name] <- quantile(gdppercapita.by.oecd[["OECD"]][, var.name], 
                                                                   probs=quantiles)
      gdppercapita.by.oecd.quant[["NotOECD"]][, var.name] <- quantile(gdppercapita.by.oecd[["NotOECD"]][, var.name], 
                                                                      probs=quantiles)
    }
  }
  
  # Carbon emissions by 5 IPCC regions
  co2.5regions <- list()
  # Uses global variables rcp.5region.names and reg.5.dataframe
  for (reg in rcp.5region.names) {
    if (reg == "REF" & year.start == 1980) next
    co2.5regions[[reg]] <- as.data.frame(matrix(nrow=length(co2.projections), 
                                                ncol=length(names(co2.projections[[1]])[-1]),
                                                dimnames=list(c(),names(co2.projections[[1]])[-1])))
    reg.isos <- subset(reg.5.dataframe, reg == Reg5)$Isocode
    which.reg.rows <- co2.projections[[1]]$Isocode %in% reg.isos
    for (i in 1:length(co2.projections)) {
      # Remove the isocode column and sum up over each year
      if (year.start < 1990) {
        # Namibia has missing CO2 data from 1980 to 1989. Just exclude
        # them from the count
        co2.5regions[[reg]][i,] <- colSums(co2.projections[[i]][which.reg.rows, -1], na.rm=T)
      } else {
        co2.5regions[[reg]][i,] <- colSums(co2.projections[[i]][which.reg.rows, -1])
      }
    }
  }
  # Carbon quantiles by 5 IPCC regions
  co2.5regions.quants <- list()
  for (reg in rcp.5region.names) {
    if (reg == "REF" & year.start == 1980) next
    print(reg)
    co2.5regions.quants[[reg]] <- as.data.frame(
      matrix(nrow=length(quantiles),
             ncol=(1 + dim(co2.5regions[[reg]])[2]),
             dimnames=list(c(), c("Quantiles", names(co2.5regions[[reg]])))))
    co2.5regions.quants[[reg]]$Quantiles <- quantiles
    for (col.name in names(co2.5regions[[reg]])) {
      if (year.start < 1990) {
        # Namibia has missing CO2 data from 1980 to 1989. Just exclude
        # them from the count
        co2.5regions.quants[[reg]][, col.name] <- 
          quantile(co2.5regions[[reg]][, col.name], probs=quantiles, na.rm=T)
      } else {
        co2.5regions.quants[[reg]][, col.name] <- 
          quantile(co2.5regions[[reg]][, col.name], probs=quantiles)
      }
    }
  }
  
  # IPAT quantiles by country
  ipat.components.bycountry <- lapply(1:length(data.proj), function(i) {
    x <- data.proj[[i]]
    get.ipat.components.bycountry(x$TechData, x$GDPData, x$GDPFrontier, x$TechUSA,
                                  preds.countries.trajs[[i]], co2.projections[[i]],
                                  year.sequence=seq(year.start,year.end,5))
  })
  ipat.quantiles.bycountry <- get.ipat.quantiles.bycountry(ipat.components.bycountry,
                                                           quantiles=quantiles)
  
  # Look at gap from the frontier.
  frontier.gaps.pred <- list()
  for (iso in names(ipat.components.bycountry[[1]])) {
    if (iso != "USA") {
      frontier.gaps.pred[[iso]] <- list()
    }
  }
  for (i in 1:length(ipat.components.bycountry)) {
    us.gdp <- ipat.components.bycountry[[i]]$USA$GDPpercapita
    for (iso in names(ipat.components.bycountry[[1]])) {
      if (iso != "USA") {
        frontier.gaps.pred[[iso]][[i]] <- 
          log(us.gdp / ipat.components.bycountry[[i]][[iso]]$GDPpercapita)
      }
    }
  }
  
  if (outofsample.validate) {
    # Warning: Assuming that the quantiles are symmetric!!!!
    stopifnot(all(abs((quantiles + rev(quantiles)) - 1) < 0.0001))
    stopifnot(year.end < 2015)
    n.intervals <- floor(length(quantiles)/2)
    tests.belowmedian <- prop.belowmedian <- list()
    tests.true <- prop.true.tests <- list()
    tests.true.ipat <- prop.true.tests.ipat <- list()
    n.tests.ipat <- list()
    for (ipat.var in names(ipat.quantiles.bycountry[["USA"]][[1]])[-1]) {
      tests.belowmedian[[ipat.var]] <- prop.belowmedian[[ipat.var]] <- list()
      tests.true.ipat[[ipat.var]] <- prop.true.tests.ipat[[ipat.var]] <- list()
      n.tests.ipat[[ipat.var]] <- list()
    }
    # Have to ignore the first year b/c it's deterministic
    year.index <- 1
    for (year in seq(year.start + 5, year.end, 5)) {
      year.index <- year.index + 1
      year.str <- as.character(year)
      n.tests <- 0
      tests.true[[year.str]] <- rep(0, n.intervals)
      for (ipat.var in names(ipat.quantiles.bycountry[["USA"]][[1]])[-1]) {
        tests.belowmedian[[ipat.var]][[year.str]] <- 0
        tests.true.ipat[[ipat.var]][[year.str]] <- rep(0, n.intervals)
        n.tests.ipat[[ipat.var]][[year.str]] <- 0
      }
      
      trajs.quants.bycountry.tmp <- trajs.quants.bycountry[[paste0("CO2", year)]]
      data.true.tmp <- subset(data.medium.full, Year == year)
      for (i in 1:dim(trajs.quants.bycountry.tmp)[1]) {
        iso <- trajs.quants.bycountry.tmp[i, "Isocode"]
        if ((! restrict.isos) | iso %in% isos.restricted) {
          co2.pred.tmp <- subset(trajs.quants.bycountry.tmp, Isocode == iso)
          data.tmp1 <- subset(data.true.tmp, Isocode == iso)
          co2.true <- data.tmp1$CO2 * (data.tmp1$PopTotal*10^3) / 10^9
          for (q_ind in 1:n.intervals) {
            co2.lower.name <- paste0("CO2_", quantiles[q_ind])
            co2.upper.name <- paste0("CO2_", quantiles[length(quantiles)+1-q_ind])
            co2.lower <- co2.pred.tmp[, co2.lower.name]
            co2.upper <- co2.pred.tmp[, co2.upper.name]
            if ((!is.na(co2.lower)) & (!is.na(co2.upper))) {
              # Need to check is.na b/c of 1980 validation
              stopifnot(co2.lower < co2.upper)
              tests.true[[year.str]][q_ind] <- tests.true[[year.str]][q_ind] + 
                (co2.lower < co2.true & co2.true < co2.upper)
            }
          }
          n.tests <- n.tests + 1
        }
      }
      prop.true.tests[[year.str]] <- tests.true[[year.str]] / n.tests
      
      for (ipat.var in names(ipat.quantiles.bycountry[["USA"]][[1]])[-1]) {
        for (i in 1:dim(trajs.quants.bycountry.tmp)[1]) {
          iso <- as.character(trajs.quants.bycountry.tmp[i, "Isocode"])
          if ((! restrict.isos) | iso %in% isos.restricted) {
            data.tmp1 <- subset(data.true.tmp, Isocode == iso & Year == year)
            
            if (ipat.var == "GDPpercapita") {
              var.true <- data.tmp1["GDP"]
            } else if (ipat.var == "GDP") {
              var.true <- data.tmp1["GDP"] * (10^3 * data.tmp1["PopTotal"])
            } else if (ipat.var == "CO2") {
              var.true <- data.tmp1["CO2"] * (10^3 * data.tmp1["PopTotal"])
            } else if (ipat.var == "CO2percapita") {
              var.true <- data.tmp1["CO2"]
            } else if (ipat.var == "Pop") {
              var.true <- data.tmp1["PopTotal"]
            } else if (ipat.var == "FrontierGap") {
              next
            } else { # Tech
              var.true <- data.tmp1[ipat.var]
            }
            if (!is.na(var.true)) {
              q_med <- n.intervals+1
              var.median <- ipat.quantiles.bycountry[[iso]][[q_med]][year.index, ipat.var]
              if ((!is.na(var.true)) & (!is.na(var.median))) {
                # Have to have this check because of 1980 out of sample validation
                tests.belowmedian[[ipat.var]][[year.str]] <- tests.belowmedian[[ipat.var]][[year.str]] + 
                  (var.true < var.median)
              }
              for (q_ind in 1:n.intervals) {
                var.lower.name <- paste0(ipat.var, "_", quantiles[q_ind])
                var.upper.name <- paste0(ipat.var, "_", quantiles[length(quantiles)+1-q_ind])
                var.lower <- ipat.quantiles.bycountry[[iso]][[q_ind]][year.index, ipat.var]
                var.upper <- ipat.quantiles.bycountry[[iso]][[length(quantiles)+1-q_ind]][year.index, ipat.var]
                if ((!is.na(var.lower)) & (!is.na(var.upper))) {
                  # Need to check is.na b/c of 1980 validation
                  stopifnot(var.lower < var.upper)
                  tests.true.ipat[[ipat.var]][[year.str]][q_ind] <- tests.true.ipat[[ipat.var]][[year.str]][q_ind] + 
                    (var.lower < var.true & var.true < var.upper)
                }
              }
              if ((!is.na(var.true)) & (!is.na(var.median))) {
                n.tests.ipat[[ipat.var]][[year.str]] <- n.tests.ipat[[ipat.var]][[year.str]] + 1
              }
            }
          }
        }
        prop.true.tests.ipat[[ipat.var]][[year.str]] <- tests.true.ipat[[ipat.var]][[year.str]] /
          n.tests.ipat[[ipat.var]][[year.str]]
        prop.belowmedian[[ipat.var]][[year.str]] <- tests.belowmedian[[ipat.var]][[year.str]] /
          n.tests.ipat[[ipat.var]][[year.str]]
      }
    }
  }
  
  if (outofsample.validate) {
    return(list(trajs.worldwide=trajs.worldwide, trajs.quants=trajs.quants, 
                trajs.quants.bycountry=trajs.quants.bycountry,
                trajs.ssa=trajs.ssa,
                ipat.ssa.quantiles=ipat.ssa.quantiles,
                co2.ssa=co2.ssa,
                co2.ssa.quants=co2.ssa.quants,
                gdp.ssa=gdp.ssa,
                gdp.ssa.quants=gdp.ssa.quants,
                co2.5regions=co2.5regions,
                co2.5regions.quants=co2.5regions.quants,
                ipat.quantiles=ipat.quantiles,
                ipat.components.bycountry=ipat.components.bycountry,
                frontier.gaps.pred=frontier.gaps.pred,
                ipat.quantiles.bycountry=ipat.quantiles.bycountry,
                quantiles=quantiles,
                prop.true.tests=prop.true.tests,
                prop.true.tests.ipat=prop.true.tests.ipat,
                prop.belowmedian=prop.belowmedian))
  } else {
    return(list(trajs.worldwide=trajs.worldwide, trajs.quants=trajs.quants, 
                trajs.quants.bycountry=trajs.quants.bycountry,
                trajs.ssa=trajs.ssa,
                ipat.ssa.quantiles=ipat.ssa.quantiles,
                cum.emissions=cum.emissions, cum.emissions.quantile=cum.emissions.quantile,
                co2.ssa=co2.ssa,
                co2.ssa.quants=co2.ssa.quants,
                gdp.ssa=gdp.ssa,
                gdp.ssa.quants=gdp.ssa.quants,
                co2.cum.ssa=co2.cum.ssa,
                cum.emissions.ssa.quantile=cum.emissions.ssa.quantile,
                cum.emissions.notssa=cum.emissions.notssa,
                cum.emissions.notssa.quantile=cum.emissions.notssa.quantile,
                gdp.by.oecd=gdp.by.oecd,
                gdp.by.oecd.quant=gdp.by.oecd.quant,
                gdppercapita.by.oecd=gdppercapita.by.oecd,
                gdppercapita.by.oecd.quant=gdppercapita.by.oecd.quant,
                co2.5regions=co2.5regions,
                co2.5regions.quants=co2.5regions.quants,
                ipat.quantiles=ipat.quantiles,
                ipat.components.bycountry=ipat.components.bycountry,
                frontier.gaps.pred=frontier.gaps.pred,
                ipat.quantiles.bycountry=ipat.quantiles.bycountry,
                quantiles=quantiles))
  }
}

# ==========================================

setwd("/path/to/working/directory")
data.location <- "NatureData/"
sims.location <- "Simulations/"
plot.location <- "NatureData/Plots/"

# Get RCP Data
rcp.list <- list()
for (rcp.num in c("RCP26", "RCP45", "RCP60", "RCP85")) {
  rcp.list[[rcp.num]] <- read.csv(paste0(data.location, "IPCC_", rcp.num, "_data.csv"))
}
rcp.colors <- c("green", "red", "black", "purple")

rcp.carbon.yearly.tmp <- read.csv(paste0(data.location, "rcp_db_carbon_emissions.csv"),
                                  nrows=4)
rcp.carbon.yearly <- rcp.carbon.yearly.tmp[, c("Scenario", "Unit",
                                               paste0("X", c(2000, 2005, seq(2010, 2100, by=10))))]
rcp.carbon.yearly[, -c(1,2)] <- (11/3) * rcp.carbon.yearly[, -c(1,2)]
rcp.carbon.yearly$Unit <- rep("PgCO2/yr", 4)
names(rcp.carbon.yearly) <- gsub("X", "Carbon", names(rcp.carbon.yearly))
rcp.carbon.yearly$Scenario <- c("RCP6.0", "RCP4.5", "RCP2.6", "RCP8.5")

rcp.carbon.cum <- rcp.carbon.yearly[, c("Scenario", "Unit",
                                        paste0("Carbon", seq(2010, 2100, by=10)))]
carbon.cum <- rep(0, 4)
for (year in seq(2010, 2100, by=10)) {
  var.name.tmp <- paste0("Carbon", year)
  if (year == 2010) {
    rcp.carbon.cum[, var.name.tmp] <- carbon.cum
    carbon.cum <- carbon.cum + 5*rcp.carbon.yearly[, var.name.tmp]
  } else if (year == 2100) {
    carbon.cum <- carbon.cum + 5*rcp.carbon.yearly[, var.name.tmp]
    rcp.carbon.cum[, var.name.tmp] <- carbon.cum
  } else {
    carbon.cum <- carbon.cum + 5*rcp.carbon.yearly[, var.name.tmp]
    rcp.carbon.cum[, var.name.tmp] <- carbon.cum
    carbon.cum <- carbon.cum + 5*rcp.carbon.yearly[, var.name.tmp]
  }
}
print(xtable(rcp.carbon.cum), include.rownames=F)
print(xtable(rcp.carbon.yearly), include.rownames=F)

# RCP 8.5, yearly worldwide emissions of CO2 (fossil fuels and industry).
# From https://tntcat.iiasa.ac.at/AR5DB/dsd?Action=htmlpage&page=regions
# Converting to gigatonnes of carbon
rcp.8.5.yearly <- data.frame(Year=c(2005, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100),
                             Carbon=c(29227.000, 32727.200, 42304.533, 50743.000, 61550.867, 74083.533, 86519.400, 95194.733, 100489.033, 103901.233, 105380.367)
                             / ((11/3) * 10^3))
# Starting RCP 8.5 from 2010
rcp.8.5.cumulative <- data.frame(Year=c(2010, seq(2015, 2095, by=10), 2100),
                                 Carbon=NA)
cum.emissions.tmp <- 0
rcp.8.5.cumulative$Carbon[1] <- 0
cum.emissions.tmp <- 5*rcp.8.5.yearly$Carbon[2] # 2010 emissions until 2015
rcp.8.5.cumulative$Carbon[2] <- cum.emissions.tmp
# add on 5 years for 2010 start
for (i in 3:dim(rcp.8.5.cumulative)[1]) {
  if (rcp.8.5.cumulative$Year[i] == 2100) {
    cum.emissions.tmp <- cum.emissions.tmp + 5*rcp.8.5.yearly$Carbon[i]
  } else {
    cum.emissions.tmp <- cum.emissions.tmp + 10*rcp.8.5.yearly$Carbon[i]
  }
  rcp.8.5.cumulative$Carbon[i] <- cum.emissions.tmp
}

# RCP 8.5 GDP data, in billions of 2005USD per year
rcp.8.5.yearly.gdp <- data.frame(Year=c(2005, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100),
                                 GDP=c(42987.656, 48302.922, 63588.657, 83849.204, 110376.903, 139114.084, 165405.548, 190946.622, 215806.785, 239972.143, 262928.537))
# RCP 8.5 Tech data, tonnes of carbon per $10,000
rcp.8.5.yearly.tech <- data.frame(Year=c(2005, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100),
                                  Tech=(rcp.8.5.yearly$Carbon / rcp.8.5.yearly.gdp$GDP * 10^4))
# RCP 8.5 Population data, in millions of people in a year
rcp.8.5.pop <- data.frame(Year=c(2005, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100),
                          PopTotal=c(6469.985, 6884.830, 7814.060, 8712.230, 9548.040, 10245.410, 10864.550, 11392.090, 11847.140, 12202.210, 12386.320))

# OECD world GDP forecasts, taken from https://data.oecd.org/gdp/gdp-long-term-forecast.htm
# In millions of USD per year (2010USD)
oecd.gdp.tmp <- read.csv(paste0(data.location, "oecd_worldgdp_forecast.csv"),
                         sep=",", stringsAsFactors=FALSE)
oecd.gdp <- oecd.gdp.tmp[, c("TIME", "Value")]
names(oecd.gdp) <- c("Year", "GDP")
oecd.gdp$GDP <- as.numeric(oecd.gdp$GDP)

# Load dataset:
load(file=paste0(data.location, "data_medium.Rda"))


# For help in calling the plot function for different countries
countrylist <- as.character(unique(data.medium$Isocode))
isolist <- countrycode(sourcevar=countrylist, origin="iso3c", destination="country.name")
cbind(countrylist, isolist)

getIso <- function(country) {
  library(countrycode)
  # take the first element in case we have a vector of strings
  if (nchar(country[1]) == 3) return(country)
  else return(countrycode(sourcevar=country, origin="country.name", destination="iso3c"))
}

getCountry <- function(iso) {
  library(countrycode)
  # take the first element in case we have a vector of strings
  if (nchar(iso[1]) > 3) return(iso)
  else return(countrycode(sourcevar=iso, origin="iso3c", destination="country.name"))
}

countries.isos <- unique(data.medium$Isocode)


countries.eu <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", 
                  "Czech Republic", "Denmark", "Estonia", "Finland", "France", 
                  "Germany", "Greece", "Hungary", "Ireland", "Italy", "Latvia",
                  "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland",
                  "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", 
                  "Sweden", "United Kingdom")
isos.eu <- sapply(countries.eu, getIso)

rcp.5region.tmp <- read.csv(paste0(data.location, 'rcp_db_5regions.csv'),
                            nrows=20)
rcp.5region <- rcp.5region.tmp[, c("Region", "Scenario", "Unit",
                                   paste0("X", c(2000, 2005, seq(2010, 2100, by=10))))]
names(rcp.5region) <- gsub("X", "Carbon", names(rcp.5region))
rcp.5region[,-c(1:3)] <- (11/3) * rcp.5region[,-c(1:3)]
rcp.5region$Unit <- "PgCO2/yr"
rcp.5region.names <- c("ASIA", "LAM", "MAF", "OECD", "REF")

# World population data:
library(wpp2015)
data(UNlocations)
UNlocations$Isocode <- sapply(as.character(UNlocations$name), getIso)


#=========================================================


oecd.world.iso <- c("AUS", "AUT", "BEL", "BRA", "CAN", "CHE", "CHL", "CHN",
                    "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA",
                    "GBR", "GRC", "HUN", "IDN", "IND", "IRL", "ISL", "ISR",
                    "ITA", "JPN", "KOR", "LUX", "MEX", "NLD", "NOR", "NZL",
                    "POL", "PRT", "RUS", "SVK", "SVN", "SWE", "TUR",
                    "USA", "ZAF")
oecd.world.minusoecd15 <- c("BRA", "CHN", "IDN", "IND", "RUS", "ZAF")
oecd.15.iso <- oecd.world.iso[! oecd.world.iso %in% oecd.world.minusoecd15]

# ========================================================

# Sub-Saharan Africa
ssa.countries <- c("Angola", "Benin", "Botswana",
                   "Burkina Faso", "Burundi", "Cameroon", "Cape Verde", "Central African Republic",
                   "Chad", "Comoros", "Congo, Republic of", "Congo, the Democratic Republic of the", "C么te d'Ivoire",
                   "Djibouti", "Equatorial Guinea", "Eritrea", "Ethiopia", "Gabon", "The Gambia", "Ghana",
                   "Guinea", "Guinea-Bissau", "Kenya", "Lesotho", "Liberia", "Madagascar", "Malawi", "Mali",
                   "Mauritania", "Mauritius", "Mozambique", "Namibia", "Niger", "Nigeria", "Réunion", "Rwanda",
                   "Sao Tome and Principe", "Senegal", "Seychelles", "Sierra Leone", "Somalia", "South Africa", "Sudan",
                   "Swaziland", "Tanzania", "Togo", "Uganda", "Western Sahara", "Zambia", "Zimbabwe")
ssa.isos <- getIso(ssa.countries)


get.rcp.region.isos <- function(countries.str) {
  # get.rcp.region.isos splits a string of countries, gets isocodes
  library(countrycode)
  return(sapply(strsplit(countries.str, ", "), getIso))
}


# IPCC RCP regions
#OECD90 = Includes the OECD 90 countries, therefore encompassing the countries included in the regions Western Europe (Austria, Belgium, Denmark, Finland, France, Germany, Greece, Iceland, Ireland, Italy, Luxembourg, Netherlands, Norway, Portugal, Spain, Sweden, Switzerland, Turkey, United Kingdom), Northern America (Canada, United States of America) and Pacific OECD (Australia, Fiji, French Polynesia, Guam, Japan, New Caledonia, New Zealand, Samoa, Solomon Islands, Vanuatu) .
oecd.90.countries.str <- "Austria, Belgium, Denmark, Finland, France, Germany, Greece, Iceland, Ireland, Italy, Luxembourg, Netherlands, Norway, Portugal, Spain, Sweden, Switzerland, Turkey, United Kingdom, Canada, United States of America, Australia, Fiji, French Polynesia, Guam, Japan, New Caledonia, New Zealand, Samoa, Solomon Islands, Vanuatu"
oecd.90.iso <- get.rcp.region.isos(oecd.90.countries.str)
#REF = Countries from the Reforming Ecomonies region (Albania, Armenia, Azerbaijan, Belarus, Bosnia and Herzegovina, Bulgaria, Croatia, Cyprus, Czech Republic, Estonia, Georgia, Hungary, Kazakhstan, Kyrgyzstan, Latvia, Lithuania, Malta, Poland, Republic of Moldova, Romania, Russian Federation, Slovakia, Slovenia, Tajikistan, TFYR Macedonia, Turkmenistan, Ukraine, Uzbekistan, Yugoslavia).
ref.countries.str <- "Albania, Armenia, Azerbaijan, Belarus, Bosnia and Herzegovina, Bulgaria, Croatia, Cyprus, Czech Republic, Estonia, Georgia, Hungary, Kazakhstan, Kyrgyzstan, Latvia, Lithuania, Malta, Poland, Republic of Moldova, Romania, Russian Federation, Slovakia, Slovenia, Tajikistan, TFYR Macedonia, Turkmenistan, Ukraine, Uzbekistan, Yugoslavia"
ref.iso <- get.rcp.region.isos(ref.countries.str)
#ASIA = The countries included in the regions China + (China, China Hong Kong SAR, China Macao SAR, Mongolia, Taiwan) , India + (Afghanistan, Bangladesh, Bhutan, India, Maldives, Nepal, Pakistan, Sri Lanka) and Rest of Asia (Brunei Darussalam, Cambodia, Democratic People's Republic of Korea, East Timor, Indonesia, Lao People's Democratic Republic, Malaysia, Myanmar, Papua New Guinea, Philippines, Republic of Korea, Singapore, Thailand, Viet Nam) are aggregated into this region.
asia.countries.str <- "China, China Hong Kong SAR, China Macao SAR, Mongolia, Taiwan, India, Afghanistan, Bangladesh, Bhutan, India, Maldives, Nepal, Pakistan, Sri Lanka, Brunei Darussalam, Cambodia, Democratic People's Republic of Korea, East Timor, Indonesia, Lao People's Democratic Republic, Malaysia, Myanmar, Papua New Guinea, Philippines, Republic of Korea, Singapore, Thailand, Viet Nam"
asia.iso <- get.rcp.region.isos(asia.countries.str)
#MAF = This region includes the Middle East (Bahrain, Iran (Islamic Republic of), Iraq, Israel, Jordan, Kuwait, Lebanon, Oman, Qatar, Saudi Arabia, Syrian Arab Republic, United Arab Emirates, Yemen) and African (Algeria, Angola, Benin, Botswana, Burkina Faso, Burundi, Cote d'Ivoire, Cameroon, Cape Verde, Central African Republic, Chad, Comoros, Congo, Democratic Republic of the Congo, Djibouti, Egypt, Equatorial Guinea, Eritrea, Ethiopia, Gabon, Gambia, Ghana, Guinea, Guinea-Bissau, Kenya, Lesotho, Liberia, Libyan Arab Jamahiriya, Madagascar, Malawi, Mali, Mauritania, Mauritius, Morocco, Mozambique, Namibia, Niger, Nigeria, Reunion, Rwanda, Senegal, Sierra Leone, Somalia, South Africa, Sudan, Swaziland, Togo, Tunisia, Uganda, United Republic of Tanzania, Western Sahara, Zambia, Zimbabwe) countries.
maf.countries.str <- "Bahrain, Iran (Islamic Republic of), Iraq, Israel, Jordan, Kuwait, Lebanon, Oman, Qatar, Saudi Arabia, Syrian Arab Republic, United Arab Emirates, Yemen, Algeria, Angola, Benin, Botswana, Burkina Faso, Burundi, Cote d'Ivoire, Cameroon, Cape Verde, Central African Republic, Chad, Comoros, Congo, Democratic Republic of the Congo, Djibouti, Egypt, Equatorial Guinea, Eritrea, Ethiopia, Gabon, Gambia, Ghana, Guinea, Guinea-Bissau, Kenya, Lesotho, Liberia, Libyan Arab Jamahiriya, Madagascar, Malawi, Mali, Mauritania, Mauritius, Morocco, Mozambique, Namibia, Niger, Nigeria, Reunion, Rwanda, Senegal, Sierra Leone, Somalia, South Africa, Sudan, Swaziland, Togo, Tunisia, Uganda, United Republic of Tanzania, Western Sahara, Zambia, Zimbabwe"
maf.iso <- get.rcp.region.isos(maf.countries.str)
#LAM = This region includes the Latin American countries (Argentina, Bahamas, Barbados, Belize, Bolivia, Brazil, Chile, Colombia, Costa Rica, Cuba, Dominican Republic, Ecuador, El Salvador, Guadeloupe, Guatemala, Guyana, Haiti, Honduras, Jamaica, Martinique, Mexico, Netherlands Antilles, Nicaragua, Panama, Paraguay, Peru, Puerto Rico, Suriname, Trinidad and Tobago, Uruguay, Venezuela). 
lam.countries.str <- "Argentina, Bahamas, Barbados, Belize, Bolivia, Brazil, Chile, Colombia, Costa Rica, Cuba, Dominican Republic, Ecuador, El Salvador, Guadeloupe, Guatemala, Guyana, Haiti, Honduras, Jamaica, Martinique, Mexico, Netherlands Antilles, Nicaragua, Panama, Paraguay, Peru, Puerto Rico, Suriname, Trinidad and Tobago, Uruguay, Venezuela"
lam.iso <- get.rcp.region.isos(lam.countries.str)


oecd.world.minusoecd90 <- oecd.world.iso[! oecd.world.iso %in% oecd.90.iso]

reg.5.dataframe <- data.frame(Isocode=c(oecd.90.iso, ref.iso, asia.iso, maf.iso, lam.iso),
                              Reg5=c(rep("OECD", length(oecd.90.iso)),
                                     rep("REF", length(ref.iso)),
                                     rep("ASIA", length(asia.iso)),
                                     rep("MAF", length(maf.iso)),
                                     rep("LAM", length(lam.iso))))

rcp.5region.names <- c("ASIA", "LAM", "MAF", "OECD", "REF")
rcp.5region.names.altorder <- c("OECD", "ASIA", "LAM", "MAF", "REF")


# Note: we're just ignoring missing data here, which is just a few countries.
olddata.5region <- data.frame(Region=rcp.5region.names)
for (year in 1960:2010) {
  var.name <- paste0("Carbon", year)
  olddata.5region[, var.name] <- rep(NA, 5)
  i.row <- 0
  for (reg in rcp.5region.names) {
    i.row <- i.row + 1
    reg.isos <- subset(reg.5.dataframe, Reg5 == reg)$Isocode
    data.tmp <- subset(data.medium, Isocode %in% reg.isos & Year == year)
    olddata.5region[i.row, var.name] <- sum(data.tmp$CO2 * data.tmp$PopTotal * 1000,
                                            na.rm=T)
  }
}

#========================================================

year.start <- 2010
year.end <- 2100
n.trajectories <- 1000

preds.countries.tmp <- predict.population(year.present=year.start, year.end=year.end,
                                               make.new=TRUE)
pred.pop <- preds.countries.tmp$pred.pop
names.countries.tmp <- as.character(unique(subset(data.medium, Isocode != "USA")$Isocode))
preds.countries.trajs <- find.pop.trajectories(data.medium, pred.pop, n.trajectories,
                                               c("USA", names.countries.tmp), year.start=year.start, year.end=year.end)
save(preds.countries.trajs, file=paste0(data.location, "poppreds_formatted_2010_2100.rda"))
load(paste0(data.location, "poppreds_formatted_2010_2100.rda"))

n.trajectories.larger <- 3125
preds.countries.trajs.larger <- find.pop.trajectories(data.medium, pred.pop, n.trajectories.larger,
                                                      c("USA", names.countries.tmp), year.start=year.start, year.end=year.end)
save(preds.countries.trajs.larger, file=paste0(data.location, "poppreds_formatted_2010_2100_larger.rda"))
load(paste0(data.location, "poppreds_formatted_2010_2100_larger.rda"))

preds.countries.tmp.2000 <- predict.population(year.present=2000, year.end=2010,
                                               make.new=TRUE)
pred.pop.2000 <- preds.countries.tmp.2000$pred.pop

preds.countries.trajs.2000 <- find.pop.trajectories(data.medium, pred.pop.2000, n.trajectories,
                                                    c("USA", names.countries.tmp),
                                                    year.start=2000, year.end=2010)

save(preds.countries.trajs.2000, file=paste0(data.location, "poppreds_formatted_2000_2010.rda"))

preds.countries.tmp.1990 <- predict.population(year.present=1990, year.end=2010,
                                               make.new=TRUE)
pred.pop.1990 <- preds.countries.tmp.1990$pred.pop
preds.countries.trajs.1990 <- find.pop.trajectories(data.medium, pred.pop.1990, n.trajectories,
                                                    c("USA", names.countries.tmp),
                                                    year.start=1990, year.end=2010)
save(preds.countries.trajs, file=paste0(data.location, "poppreds_formatted_1990_2010.rda"))

preds.countries.tmp <- predict.population(year.present=1980, year.end=2010,
                                          make.new=TRUE)
pred.pop <- preds.countries.tmp$pred.pop
preds.countries.trajs <- find.pop.trajectories(data.medium, pred.pop, n.trajectories,
                                               c("USA", names.countries.tmp),
                                               #year.start=1980, year.end=2010)
                                               year.start=1980, year.end=2010)

save(preds.countries.trajs, file=paste0(data.location, "poppreds_formatted_1980_2010.rda"))
#rm(preds.countries.trajs)
#rm(pred.pop)

# =======================================


model.ar1.const.results.2010 <- fit.project.model(data.medium, fit.ar1.const=T)

save(model.ar1.const.results.2010, file=paste0(data.location, "model_results_ar1const.rda"))
load(file=paste0(data.location, "model_results_ar1const.rda"))

model.ar1.const.results.2000 <- fit.project.model(data.medium, year.start=2000,
                                                  year.end=2010, fit.ar1.const=T)
#                                                  year.end=2100, fit.ar1.const=T)
save(model.ar1.const.results.2000, file=paste0(data.location, "model_results_ar1const_2000.rda"))
load(file=paste0(data.location, "model_results_ar1const_2000.rda"))
model.ar1.const.results.1990 <- fit.project.model(data.medium, year.start=1990,
                                                  year.end=2010, fit.ar1.const=T)
save(model.ar1.const.results.1990, file=paste0(data.location, "model_results_ar1const_1990.rda"))
#load(file=paste0(data.location, "model_results_ar1const_1990.rda"))
pct <- proc.time()
model.ar1.const.results.1980 <- fit.project.model(data.medium, year.start=1980,
                                                  year.end=2010, fit.ar1.const=T)
#                                                  year.end=2100, fit.ar1.const=T)
save(model.ar1.const.results.1980, file=paste0(data.location, "model_results_ar1const_1980.rda"))
proc.time() - pct
#load(file=paste0(data.location, "model_results_ar1const_1980.rda"))




proj.evals.2010.ar1.const <- evaluate.projections(data.medium, model.ar1.const.results.2010,
                                                  year.start=2010, year.end=2100,
                                                  outofsample.validate=F,
                                                  quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975))
save(proj.evals.2010.ar1.const, file=paste0(data.location, "proj_evals_ar1const_2010.rda"))
load(file=paste0(data.location, "proj_evals_ar1const_2010.rda"))


proj.evals.2000.ar1.const <- evaluate.projections(data.medium, model.ar1.const.results.2000,
                                                  year.start=2000, year.end=2010,
                                                  outofsample.validate=T,
                                                  quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975))
save(proj.evals.2000.ar1.const, file=paste0(data.location, "proj_evals_ar1const_2000.rda"))
load(file=paste0(data.location, "proj_evals_ar1const_2000.rda"))
proj.evals.1990.ar1.const <- evaluate.projections(data.medium, model.ar1.const.results.1990,
                                                  year.start=1990, year.end=2010,
                                                  outofsample.validate=T,
                                                  quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975))
save(proj.evals.1990.ar1.const, file=paste0(data.location, "proj_evals_ar1const_1990.rda"))
#load(file=paste0(data.location, "proj_evals_ar1const_1990.rda"))
proj.evals.1980.ar1.const <- evaluate.projections(data.medium, model.ar1.const.results.1980,
                                                  year.start=1980, year.end=2010,
                                                  outofsample.validate=T,
                                                  quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975))
save(proj.evals.1980.ar1.const, file=paste0(data.location, "proj_evals_ar1const_1980.rda"))
load(file=paste0(data.location, "proj_evals_ar1const_1980.rda"))

# ==============================================================

table.result <- as.data.frame(model.ar1.const.results.2010$corr.model.out[[2]][[1]])
table.result <- table.result[,c('gamma','mu.phi', 'sig.phi', 'sig.frontier', 'mu.ln.gap', 'sig.ln.gap',
                                'mu.beta', 'mu.eta', 'mu.delta', 'sig.delta', 'sig.eps.mean', 'sig.eps.sd', 'rho')]
# Paper plots

#proj.evals.2010.ar1.const$cum.emissions.notssa.quantile[19,]
#proj.evals.2010.ar1.const$cum.emissions.quantile[19,]
#proj.evals.2010.ar1.const$cum.emissions.notssa.quantile[19,] / proj.evals.2010.ar1.const$cum.emissions.quantile[19,]
#   Quant0.025  Quant0.1  Quant0.5  Quant0.9 Quant0.975 Year
#19  0.9596292 0.9483812 0.9469967 0.9614328   0.963868    1

isos.plot <- c("USA", "GBR", "JPN", "BRA", "KOR", "NGA", "IND", "CHN")
countries.plot <- c("USA", "UK", "Japan", "Brazil", "S. Korea", "Nigeria", "India", "China")


subset(data.medium, Year == 2008 & Tech < 0.5)[, c("Isocode", "Country", "Tech", "PopTotal", "GDP")]
dim(subset(data.medium, Year == 2008 & Tech < 0.5)[, c("Isocode", "Country", "Tech", "PopTotal", "GDP")])
sum(subset(data.medium, Year == 2008 & Tech < 0.5)[, c("Isocode", "Country", "Tech", "PopTotal", "GDP")]$GDP < 2000)

subset(data.medium, Year == 2008 & Tech > 2)[, c("Isocode", "Country", "Tech", "PopTotal", "GDP")]
dim(subset(data.medium, Year == 2008 & Tech > 2)[, c("Isocode", "Country", "Tech", "PopTotal", "GDP")])

data.tmp <- subset(data.medium, Isocode %in% isos.plot)
data.tmp.first <- subset(data.tmp, Year == 1960)
data.tmp.first$Country <- gsub("\\.", " ", data.tmp.first$Country)
data.tmp.first$Country[3] <- "UK"
data.tmp.first$Country[6] <- "S. Korea"
data.tmp.first$Country[8] <- "USA"
library(ggplot2)
# using http://colorbrewer2.org/
Palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
font <- "Times"

pdf(paste0(plot.location, "fig-GDP.pdf"))
ggplot(data.tmp, aes(x=Year, y=GDP, group=Isocode,
                     color=Isocode)) +
  geom_line() +
  # Make the USA stand out
  geom_line(data=subset(data.tmp, Isocode == "USA"), aes(x=Year, y=GDP),
            color="black", size=0.7) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.minor=element_line(color="white"),
        panel.grid.major=element_line(color="grey90")) +
  scale_colour_manual(values=Palette, guide=FALSE) +
  scale_y_log10(breaks=c(500, 1000, 2000, 10000, 20000)) +
  ylab("GDP per capita (in 2010 USD, on log scale)") +
  scale_x_continuous(limits=c(1955, 2010), breaks=seq(1960, 2010, by=10)) +
  geom_text(data=data.tmp.first, aes(label=Country), hjust=1, family=font) +
  # Make the USA stand out
  geom_text(data=subset(data.tmp.first, Isocode == "USA"),
            aes(label=Country), hjust=1, fontface="bold", color="black", family=font)
dev.off()

# Now do the same thing with Tech
pdf(paste0(plot.location, "fig-intensity.pdf"))
ggplot(data.tmp, aes(x=Year, y=Tech, group=Isocode,
                     color=Isocode)) +
  geom_line() +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.minor=element_line(color="white"),
        panel.grid.major=element_line(color="grey90")) +
  scale_colour_manual(values=Palette, guide=FALSE) +
  scale_y_continuous(limits=c(0,12)) +
  ylab("Carbon Intensity (tonnes of carbon per 10,000 USD)") +
  scale_x_continuous(limits=c(1955, 2010), breaks=seq(1960, 2010, by=10)) +
  geom_text(data=data.tmp.first, aes(label=Country), hjust=1, family=font)
dev.off()

plot.intensity.trend <- function(data.tmp, data.tmp.first, iso.plot) {
  # Warning: font is a global variable
  data.tmp2 <- subset(data.tmp, Isocode == iso.plot)
  data.tmp.first2 <- subset(data.tmp.first, Isocode == iso.plot)
  intensity.plot <- ggplot(data.tmp2, aes(x=Year, y=Tech, group=Isocode)) +
    geom_line() +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90")) +
    theme(text=element_text(family=font)) +
    scale_y_continuous(limits=c(0,max(data.tmp2$Tech))) +
    ylab("CO2 Intensity") +
    scale_x_continuous(limits=c(1947, 2010), breaks=seq(1960, 2010, by=10)) +
    geom_text(data=data.tmp.first2, aes(label=Country, family=font), hjust=1)
  intensity.plot
}

library(gridExtra)
font <- "Times"
intensity.plots <- list()
for (iso.tmp in c("USA", "CHN", "IND", "NGA")) {
  intensity.plots[[iso.tmp]] <- plot.intensity.trend(data.tmp, data.tmp.first, iso.tmp)
}
pdf(paste0(plot.location, "fig-intensity-grid.pdf"))
do.call(grid.arrange, intensity.plots)
dev.off()


var.decomp.tmp <- find.ipat.vars(proj.evals.2010.ar1.const$trajs.worldwide)
library(xtable)
print(xtable(var.decomp.tmp, digits=5), include.rownames=F)

plot.co2.projections <- function(data.medium,
                                 trajs.quants,
                                 ylabel="Yearly CO2 emissions (gt CO2)",
                                 plot.rcp=TRUE,
                                 ybreaks=seq(0,120,by=20)) {
  # plot.co2.projections plots projections of yearly CO2 emissions worldwide
  
  library(doBy)
  library(ggplot2)
  data.restriction <- data.medium
  
  co2.sums.na.rm <- summaryBy(CO2*PopTotal*1000 ~ Year, data.restriction, FUN=function(x) sum(x, na.rm=T))
  names(co2.sums.na.rm)[2] <- "Carbon"
  co2.sums.na.rm$Carbon <- co2.sums.na.rm$Carbon / 10^9
  
  year.start <- 2010
  year.end <- 2100
  x.vals <- seq(year.start, year.end, by=5)
  
  # Create plot
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=22)) +
    scale_y_continuous(breaks=ybreaks) +
    scale_x_continuous(breaks=seq(1960,2100,by=20))
  
  trajs.quants.tmp <- t(trajs.quants)
  trajs.quants.transp <- data.frame(trajs.quants.tmp[-1, ])
  names(trajs.quants.transp) <- paste0("Quant", trajs.quants.tmp[1,])
  trajs.quants.transp$Year <- seq(2010, 2100, by=5)
  
  # Plot our projections
  color.us <- "#e41a1c"
  plot.obj <- plot.obj +
    ylab(ylabel) +
    geom_ribbon(aes(x=Year, ymin=Quant0.05, ymax=Quant0.95), data=trajs.quants.transp,
                alpha=0.3, fill=color.us) +
    geom_ribbon(aes(x=Year, ymin=Quant0.025, ymax=Quant0.975), data=trajs.quants.transp,
                alpha=0.2, fill=color.us) +
    geom_line(data=trajs.quants.transp, aes(x=Year, y=Quant0.5,
                                            color="Projections"),
              alpha=1, size=1.3, color=color.us)
  
  if (plot.rcp) {
    # Plot RCP numbers
    years.tmp <- seq(2010, 2100, by=10)
    rcp.yearly.tmp <- data.frame(sapply(data.frame(t(rcp.carbon.yearly)[-c(1:4),],
                                                   stringsAsFactors=F),
                                        as.numeric))
    names(rcp.yearly.tmp) <- rcp.carbon.yearly$Scenario
    rcp.yearly.tmp$Year <- years.tmp
    library(reshape2)
    rcp.yearly.long <- melt(rcp.yearly.tmp, id.vars=c("Year"))
    plot.obj <- plot.obj + geom_line(data=rcp.yearly.long,
                                     aes(x=Year, y=value, group=variable,
                                         color=variable), size=0.7, alpha=0.9,
                                     linetype=2) +
      scale_color_manual(limits=c("RCP8.5", "RCP6.0", "RCP4.5", "RCP2.6"),
                         values=c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00"))
  }
  
  # Plot historical data
  plot.obj <- plot.obj + geom_line(data=co2.sums.na.rm,
                                   aes(x=Year,y=Carbon))
  
  plot.obj
}

plot.co2.projections.cum <- function(data.medium,
                                     cum.emissions.quant,
                                     ylabel="Cumulative CO2 Emissions (gt CO2)") {
  # plot.co2.projections.cum plots projections of cumulative CO2 emissions worldwide
  
  library(doBy)
  library(ggplot2)
  data.restriction <- data.medium
  
  year.start <- 2010
  year.end <- 2100
  x.vals <- seq(year.start, year.end, by=5)
  
  # Create plot
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=22))
  
  color.us <- "#e41a1c"
  plot.obj <- plot.obj +
    ylab(ylabel) +
    geom_ribbon(aes(x=Year, ymin=Quant0.05, ymax=Quant0.95), data=cum.emissions.quant,
                alpha=0.3, fill=color.us) +
    geom_ribbon(aes(x=Year, ymin=Quant0.025, ymax=Quant0.975), data=cum.emissions.quant,
                alpha=0.2, fill=color.us) +
    geom_line(data=cum.emissions.quant, aes(x=Year, y=Quant0.5),
              alpha=1, size=1.3, color=color.us) +
    scale_x_continuous(breaks=c(2010, 2025, 2050, 2075, 2100), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))
  
  # Plot RCP numbers
  years.tmp <- seq(2010, 2100, by=10)
  rcp.cum.tmp <- data.frame(sapply(data.frame(t(rcp.carbon.cum)[-c(1:2),],
                                              stringsAsFactors=F),
                                   as.numeric))
  names(rcp.cum.tmp) <- rcp.carbon.cum$Scenario
  rcp.cum.tmp$Year <- years.tmp
  library(reshape2)
  rcp.cum.long <- melt(rcp.cum.tmp, id.vars=c("Year"))
  plot.obj <- plot.obj + geom_line(data=rcp.cum.long,
                                   aes(x=Year, y=value, group=variable,
                                       color=variable), size=0.7, alpha=0.9,
                                   linetype=2) +
    scale_color_manual(limits=c("RCP8.5", "RCP6.0", "RCP4.5", "RCP2.6"),
                       values=c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00"))
  
  plot.obj
}

plot.co2 <- plot.co2.projections(data.medium,
                                 proj.evals.2010.ar1.const$trajs.quants)
pdf(paste0(plot.location, "carbon_yearly.pdf"))
print(plot.co2)
dev.off()


plot.co2.cum <- plot.co2.projections.cum(data.medium,
                                         proj.evals.2010.ar1.const$cum.emissions.quantile)
pdf(paste0(plot.location, "carbon_cum.pdf"))
print(plot.co2.cum)
dev.off()



plot.co2.projections.cum.ssa <- function(data.medium,
                                         cum.emissions.quant,
                                         ylabel="Cumulative CO2 Emissions (gt CO2)") {
  # plot.co2.projections.cum.ssa plots projections of yearly cumulative
  # CO2 emissions for Sub-Saharan Africa
  
  library(doBy)
  library(ggplot2)
  data.restriction <- subset(data.medium, Isocode %in% ssa.isos)
  
  year.start <- 2010
  year.end <- 2100
  x.vals <- seq(year.start, year.end, by=5)
  
  # Create plot
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=22))
  
  color.us <- "#e41a1c"
  plot.obj <- plot.obj +
    ylab(ylabel) +
    geom_ribbon(aes(x=Year, ymin=Quant0.05, ymax=Quant0.95), data=cum.emissions.quant,
                alpha=0.3, fill=color.us) +
    geom_ribbon(aes(x=Year, ymin=Quant0.025, ymax=Quant0.975), data=cum.emissions.quant,
                alpha=0.2, fill=color.us) +
    geom_line(data=cum.emissions.quant, aes(x=Year, y=Quant0.5),
              alpha=1, size=1.3, color=color.us) +
    scale_x_continuous(breaks=c(2010, 2025, 2050, 2075, 2100), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))
  
  plot.obj
}

plot.co2.cum.ssa <- plot.co2.projections.cum.ssa(data.medium,
                                                 proj.evals.2010.ar1.const$cum.emissions.ssa.quantile)

pdf(paste0(plot.location, "ssa_carbon_cum.pdf"))
print(plot.co2.cum.ssa)
dev.off()

plot.co2.yearly <- plot.co2.projections(data.medium,
                                        proj.evals.2010.ar1.const$trajs.quants)
print(plot.co2.yearly)

pdf(paste0(plot.location, "carbon_yearly.pdf"))
print(plot.co2.yearly)
dev.off()

data.ssa <- subset(data.medium, Isocode %in% ssa.isos)
plot.co2.ssa.yearly <- plot.co2.projections(data.ssa,
                                            proj.evals.2010.ar1.const$co2.ssa.quants,
                                            ylab="Sub-Saharan African emissions (gt CO2)",
                                            plot.rcp=F,
                                            ybreaks=seq(0,16,2))
print(plot.co2.ssa.yearly)

pdf(paste0(plot.location, "ssa_carbon_yearly.pdf"))
print(plot.co2.ssa.yearly)
dev.off()


plot.outofsample.validate <- function(data.medium,
                                      trajs.quants,
                                      year.start) {
  # plot.outofsample.validate makes plots of projections of worldwide yearly CO2
  # emissions when performing out of sample validation, along with what actually happened.
  
  library(ggplot2)
  library(doBy)
  trajs.quants.transp <- data.frame(t(trajs.quants)[-1,])
  names(trajs.quants.transp) <- paste0("Quant", c(0.025, 0.05, 0.5, 0.95, 0.975))
  year.seq <- seq(year.start, 2010, by=5)
  trajs.quants.transp$Year <- year.seq
  
  if (year.start >= 1990) {
    data.restriction <- data.medium
  } else {
    data.medium.1980 <- subset(data.medium, Year == 1980)
    isos.include <- data.medium.1980$Isocode[! is.na(data.medium.1980$Tech)]
    data.restriction <- subset(data.medium, Isocode %in% isos.include)
  }
  
  co2.sums.na.rm <- summaryBy(CO2*PopTotal*1000 ~ Year, data.restriction, FUN=function(x) sum(x, na.rm=T))
  names(co2.sums.na.rm)[2] <- "Carbon"
  co2.sums.na.rm$Carbon <- co2.sums.na.rm$Carbon / 10^9
  
  # Create plot
  color.us <- "#e41a1c"
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=22)) +
    scale_x_continuous(breaks=seq(1960, 2010, by=10))
  
  plot.obj <- plot.obj +
    ylab("CO2 emissions (gt CO2)") +
    geom_ribbon(aes(x=Year, ymin=Quant0.05, ymax=Quant0.95), data=trajs.quants.transp,
                alpha=0.3, fill=color.us) +
    geom_ribbon(aes(x=Year, ymin=Quant0.025, ymax=Quant0.975), data=trajs.quants.transp,
                alpha=0.2, fill=color.us) +
    geom_line(data=trajs.quants.transp, aes(x=Year, y=Quant0.5,
                                            color="Projections"),
              size=1, color=color.us)
  
  # Plot historical data
  plot.obj <- plot.obj + geom_line(data=co2.sums.na.rm,
                                   aes(x=Year,y=Carbon))
}

plot.co2.1980 <- plot.outofsample.validate(data.medium,
                                           proj.evals.1980.ar1.const$trajs.quants,
                                           1980)
pdf(paste0(plot.location, "validate1980.pdf"))
print(plot.co2.1980)
dev.off()
plot.co2.1990 <- plot.outofsample.validate(data.medium,
                                           proj.evals.1990.ar1.const$trajs.quants,
                                           1990)
pdf(paste0(plot.location, "validate1990.pdf"))
print(plot.co2.1990)
dev.off()
plot.co2.2000 <- plot.outofsample.validate(data.medium,
                                           proj.evals.2000.ar1.const$trajs.quants,
                                           2000)
pdf(paste0(plot.location, "validate2000.pdf"))
print(plot.co2.2000)
dev.off()

ipat.olddata <- get.ipat.pastdata(data.medium, seq(1960, 2010, 5), c(names.countries, "USA"))
ipat.olddata.rm <- get.ipat.pastdata(data.medium, seq(1960, 2010, 5), c(names.countries, "USA"),
                                     na.remove=T)
ipat.olddata.ssa.rm <- get.ipat.pastdata(data.medium, seq(1960, 2010, 5), ssa.isos,
                                         na.remove=T)
ipat.olddata.ssa <- get.ipat.pastdata(data.medium, seq(1960, 2010, 5), ssa.isos,
                                      na.remove=F)

plot.ipat.projections <- function(ipat.quants,
                                  ipat.olddata.rm) {
  # plot.ipat.projections makes plots of projections of worldwide IPAT components
  
  library(ggplot2)
  var.array <- c("Pop", "GDPpercapita", "Tech")
  quantiles <- as.character(c(0.025, 0.05, 0.5, 0.95, 0.975))
  
  # First deal with historical data
  ipat.olddata.rm.log <- ipat.olddata.rm
  base.shift <- list()
  for (var in var.array) {
    var.log <- paste0("log", var)
    base.shift[[var.log]] <- log(ipat.olddata.rm.log[, var])[1]
    ipat.olddata.rm.log[, var.log] <- log(ipat.olddata.rm.log[, var]) -
      base.shift[[var.log]]
  }
  
  # Now look at our projections
  ipat.quants.reformed <- list()
  for (var in var.array) {
    var.log <- paste0("log", var)
    ipat.quants.reformed[[var.log]] <- data.frame(Year=ipat.quants[[1]]$Year)
    for (q in quantiles) {
      ipat.quants.reformed[[var.log]][, q] <- log(ipat.quants[[q]][, var]) -
        base.shift[[var.log]]
    }
    names(ipat.quants.reformed[[var.log]]) <- c("Year", paste0("Quant", quantiles))
  }
  
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    ylab("IPAT Components, log scale") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=22)) +
    scale_x_continuous(limits=c(1960, 2120), breaks=seq(1960,2100,by=20)) +
    scale_y_continuous(breaks=seq(-3,3))
  
  color.i <- 0
  color.vars <- c("#7fc97f", "#beaed4", "#fdc086")
  for (var.log in paste0("log", var.array)) {
    color.i <- color.i + 1
    color.tmp <- color.vars[color.i]
    
    # Plot our projections
    plot.obj <- plot.obj +
      geom_ribbon(aes(x=Year, ymin=Quant0.05, ymax=Quant0.95), data=ipat.quants.reformed[[var.log]],
                  alpha=0.3, fill=color.tmp) +
      geom_ribbon(aes(x=Year, ymin=Quant0.025, ymax=Quant0.975), data=ipat.quants.reformed[[var.log]],
                  alpha=0.2, fill=color.tmp) +
      geom_line(data=ipat.quants.reformed[[var.log]], aes(x=Year, y=Quant0.5),
                size=1, color=color.tmp)
    # Plot the variable name
    if (var.log == "logGDPpercapita") {
      var.print <- "GDP per
      capita"
    } else if (var.log == "logTech") {
      var.print <- "Intensity"
    } else {
      var.print <- "Population"
    }
    data.last.row <- ipat.quants.reformed[[var.log]][dim(ipat.quants.reformed[[var.log]])[1],]
    plot.obj <- plot.obj + annotate("text",
                                    x=data.last.row$Year, y=data.last.row$Quant0.5,
                                    label=var.print, hjust=0, color="grey20", family=font)
    # Plot historical ipat data
    plot.obj <- plot.obj + geom_line(data=ipat.olddata.rm.log,
                                     aes_string(x="Year", y=var.log),
                                     color=color.tmp)
  }
  
  plot.obj
}

plot.ipat.paper <- plot.ipat.projections(proj.evals.2010.ar1.const$ipat.quantiles,
                                         ipat.olddata.rm)
pdf(paste0(plot.location, "ipat_components.pdf"))
print(plot.ipat.paper)
dev.off()



plot.ipat.ssa.projections <- function(ipat.quants,
                                      ipat.olddata.rm, ipat.olddata) {
  # plot.ipat.ssa.projections plots projections of IPAT components aggregated
  # over Sub-Saharan Africa
  
  library(ggplot2)
  var.array <- c("Pop", "GDPpercapita", "Tech")
  quantiles <- as.character(c(0.025, 0.05, 0.5, 0.95, 0.975))
  
  # First deal with historical data
  ipat.olddata.rm.log <- ipat.olddata.rm
  base.shift <- list()
  for (var in var.array) {
    var.log <- paste0("log", var)
    base.shift[[var.log]] <- log(ipat.olddata.rm.log[, var])[1]
    ipat.olddata.rm.log[, var.log] <- log(ipat.olddata.rm.log[, var]) -
      base.shift[[var.log]]
  }
  
  # Now look at our projections
  ipat.quants.reformed <- list()
  for (var in var.array) {
    var.log <- paste0("log", var)
    ipat.quants.reformed[[var.log]] <- data.frame(Year=ipat.quants[[1]]$Year)
    for (q in quantiles) {
      ipat.quants.reformed[[var.log]][, q] <- log(ipat.quants[[q]][, var]) -
        base.shift[[var.log]]
    }
    names(ipat.quants.reformed[[var.log]]) <- c("Year", paste0("Quant", quantiles))
  }
  
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    ylab("IPAT Components, log scale") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=22)) +
    scale_x_continuous(limits=c(1960, 2120), breaks=seq(1960,2100,by=20)) +
    scale_y_continuous(breaks=seq(-3,3))
  
  color.i <- 0
  color.vars <- c("#7fc97f", "#beaed4", "#fdc086")
  for (var.log in paste0("log", var.array)) {
    color.i <- color.i + 1
    color.tmp <- color.vars[color.i]
    
    # Plot our projections
    plot.obj <- plot.obj +
      geom_ribbon(aes(x=Year, ymin=Quant0.05, ymax=Quant0.95), data=ipat.quants.reformed[[var.log]],
                  alpha=0.3, fill=color.tmp) +
      geom_ribbon(aes(x=Year, ymin=Quant0.025, ymax=Quant0.975), data=ipat.quants.reformed[[var.log]],
                  alpha=0.2, fill=color.tmp) +
      geom_line(data=ipat.quants.reformed[[var.log]], aes(x=Year, y=Quant0.5),
                size=1, color=color.tmp)
    # Plot the variable name
    if (var.log == "logGDPpercapita") {
      var.print <- "GDP per
      capita"
    } else if (var.log == "logTech") {
      var.print <- "Intensity"
    } else {
      var.print <- "Population"
    }
    data.last.row <- ipat.quants.reformed[[var.log]][dim(ipat.quants.reformed[[var.log]])[1],]
    plot.obj <- plot.obj + annotate("text",
                                    x=data.last.row$Year, y=data.last.row$Quant0.5,
                                    label=var.print, hjust=0, color="grey20", family=font)
    # Plot historical ipat data
    plot.obj <- plot.obj + geom_line(data=ipat.olddata.rm.log,
                                     aes_string(x="Year", y=var.log),
                                     color=color.tmp)
  }
  
  plot.obj
}


plot.ipat.ssa.paper <- plot.ipat.projections(proj.evals.2010.ar1.const$ipat.ssa.quantiles,
                                             ipat.olddata.ssa.rm[1:10,])

pdf(paste0(plot.location, "ssa_ipat_components.pdf"))
print(plot.ipat.ssa.paper)
dev.off()


# =========================================

# For getting numbers we use in the paper

getCountry(as.character(subset(data.medium, Year == 2010 & GDP > 30000)$Isocode))
subset(data.medium, Year == 2010 & GDP > 30000)[, c("Isocode", "Country", "Tech")]

subset(data.medium, Year == 2008 & Tech < 1.3)[, c("Isocode", "Country", "Tech", "GDP")]
dim(subset(data.medium, Year == 2008 & Tech < 1.3)[, c("Isocode", "Country", "Tech", "GDP")])
sum(subset(data.medium, Year == 2008 & Tech < 1.3)$GDP < 4000)

subset(data.medium, Year == 2008 & Tech > 5)[, c("Isocode", "Country", "Tech", "GDP")]
dim(subset(data.medium, Year == 2008 & Tech > 5)[, c("Isocode", "Country", "Tech", "GDP")])


#=========================================================

getCountry(as.character(subset(data.medium, Year == 2010 & GDP > 30000)$Isocode))
subset(data.medium, Year == 2010 & GDP > 30000)[, c("Isocode", "Country", "Tech")]

subset(data.medium, Year == 2008 & Tech < 1.3)[, c("Isocode", "Country", "Tech", "GDP")]
dim(subset(data.medium, Year == 2008 & Tech < 1.3)[, c("Isocode", "Country", "Tech", "GDP")])
sum(subset(data.medium, Year == 2008 & Tech < 1.3)$GDP < 4000)

subset(data.medium, Year == 2008 & Tech > 5)[, c("Isocode", "Country", "Tech", "GDP")]
dim(subset(data.medium, Year == 2008 & Tech > 5)[, c("Isocode", "Country", "Tech", "GDP")])



# =========================================


# Making the table of variances
find.ipat.vars <- function(trajs.worldwide) {
  # find.ipat.vars finds the variances of the different IPAT components
  # worldwide on the log scale
  
  index.2050 <- which(trajs.worldwide[[1]]$Year == 2050)
  index.2100 <- which(trajs.worldwide[[1]]$Year == 2100)
  logpop.2050 <- loggdp.2050 <- logtech.2050 <- rep(NA, length(trajs.worldwide))
  logpop.2100 <- loggdp.2100 <- logtech.2100 <- rep(NA, length(trajs.worldwide))
  logco2.2100 <- co2.2100 <- logco2.2050<- co2.2050 <- rep(NA, length(trajs.worldwide))
  for (i in 1:length(trajs.worldwide)) {
    logpop.2050[i] <- log(trajs.worldwide[[i]]$Pop[index.2050])
    loggdp.2050[i] <- log(trajs.worldwide[[i]]$GDPpercapita[index.2050])
    logtech.2050[i] <- log(trajs.worldwide[[i]]$Tech[index.2050])
    logpop.2100[i] <- log(trajs.worldwide[[i]]$Pop[index.2100])
    loggdp.2100[i] <- log(trajs.worldwide[[i]]$GDPpercapita[index.2100])
    logtech.2100[i] <- log(trajs.worldwide[[i]]$Tech[index.2100])
    
    logco2.2050[i] <- log(trajs.worldwide[[i]]$CO2[index.2050])
    logco2.2100[i] <- log(trajs.worldwide[[i]]$CO2[index.2100])
    co2.2050[i] <- trajs.worldwide[[i]]$CO2[index.2050]
    co2.2100[i] <- trajs.worldwide[[i]]$CO2[index.2100]
  }
  var.decomp <- data.frame(matrix(nrow=4, ncol=3))
  names(var.decomp) <- c("Component", "VarLog2050", "VarLog2100")
  var.decomp[,1] <- c("Population", "GDP", "Intensity", "CO2")
  var.decomp[,2] <- c(var(logpop.2050), var(loggdp.2050), var(logtech.2050), var(logco2.2050))
  var.decomp[,3] <- c(var(logpop.2100), var(loggdp.2100), var(logtech.2100), var(logco2.2100))
  print("In 2050:")
  print(paste("Var(log pop):", var(logpop.2050)))
  print(paste("Var(log gdp):", var(loggdp.2050)))
  print(paste("Var(log tech):", var(logtech.2050)))
  print("In 2100:")
  print(paste("Var(log pop):", var(logpop.2100)))
  print(paste("Var(log gdp):", var(loggdp.2100)))
  print(paste("Var(log tech):", var(logtech.2100)))
  
  print("In 2050:")
  print(paste("Var(log CO2):", var(logco2.2050)))
  print(paste("Var(CO2):", var(co2.2050)))
  print("In 2100:")
  print(paste("Var(log CO2):", var(logco2.2100)))
  print(paste("Var(CO2):", var(co2.2100)))
  var.decomp
}

#find.ipat.vars(model.ar1.const.results.2010$trajs.worldwide)


summarize.ipat.country <- function(data.medium, ipat.components.bycountry,
                                   ipat.quantiles.bycountry, iso, var.name,
                                   plot.log=F,
                                   n.trajs.plot=100) {
  # summarize.ipat.country makes plots of projections for a specific IPAT
  # component for a specific country (given by iso)
  
  data.tmp <- subset(data.medium, Isocode == iso)
  
  if (var.name == "GDP") {
    y.vals <- data.tmp[, "GDP"] * (10^3 * data.tmp[, "PopTotal"])
  } else if (var.name == "GDPpercapita") {
    y.vals <- data.tmp[, "GDP"] / 10^3
  } else if (var.name == "Pop") {
    y.vals <- data.tmp[, "PopTotal"]
  } else if (var.name == "CO2") {
    y.vals <- data.tmp[, "CO2"] * (10^3 * data.tmp[, "PopTotal"]) / 10^6
  } else {
    y.vals <- data.tmp[, var.name]
  }
  ylim.tmp <- c(0, max(y.vals, na.rm=T))
  quants.names <- names(ipat.quantiles.bycountry[[iso]])
  ylim.tmp[2] <- max(c(ylim.tmp[2],
                       ipat.quantiles.bycountry[[iso]][[quants.names[length(quants.names)]]][, var.name]),
                     na.rm=T)
  if (var.name == "Pop") {
    ylab.name.tmp <- "Population (1000s)"
    main.name.tmp <- "Population"
  } else if (var.name == "Tech") {
    ylab.name.tmp <- "Carbon Intensity"
    main.name.tmp <- "Carbon Intensity"
  } else if (var.name == "GDP") {
    ylab.name.tmp <- "GDP (2010 USD)"
    main.name.tmp <- "GDP"
  } else if (var.name == "GDPpercapita") {
    ylab.name.tmp <- "GDP per capita"
    main.name.tmp <- "GDP per capita"
  } else if (var.name == "CO2") {
    ylab.name.tmp <- "CO2 Emissions per year"
    main.name.tmp <- "CO2 Emissions"
  } else {
    stop("Variable unrecognized by summarize.ipat.country")
  }
  if (plot.log) {
    if (var.name == "CO2") {
      # Divide by 10^6 to get megatonnes
      ylim.tmp <- c(min(c(y.vals, ipat.quantiles.bycountry[[iso]][[quants.names[1]]][, var.name] / 10^6), na.rm=T),
                    max(c(y.vals, ipat.quantiles.bycountry[[iso]][[quants.names[length(quants.names)]]][, var.name] / 10^6), na.rm=T))
    } else if (var.name == "GDPpercapita") {
      # Divide by 10^3 to get in 1000s
      ylim.tmp <- c(min(c(y.vals, ipat.quantiles.bycountry[[iso]][[quants.names[1]]][, var.name] / 10^3), na.rm=T),
                    max(c(y.vals, ipat.quantiles.bycountry[[iso]][[quants.names[length(quants.names)]]][, var.name] / 10^3), na.rm=T))
    } else {
      ylim.tmp <- c(min(c(y.vals, ipat.quantiles.bycountry[[iso]][[quants.names[1]]][, var.name]), na.rm=T),
                    max(c(y.vals, ipat.quantiles.bycountry[[iso]][[quants.names[length(quants.names)]]][, var.name]), na.rm=T))
    }
  }
  ylab.name <- ylab.name.tmp
  if (iso == "KOR") {
    main.tmp <- "South Korea"
  } else {
    main.tmp <- getCountry(iso)
  }
  
  if (plot.log) {
    plot(data.tmp[, "Year"], y.vals, type="l",
         yaxs="i", xlim=c(1960, 2100),
         ylim=ylim.tmp,
         log="y",
         main=main.tmp,
         xlab="Year", ylab=ylab.name,
         xaxs="i")
  } else {
    plot(data.tmp[, "Year"], y.vals, type="l",
         yaxs="i", xlim=c(1960, 2100),
         ylim=ylim.tmp,
         main=main.tmp,
         xlab="Year", ylab=ylab.name,
         xaxs="i")
  }
  
  for (q in quants.names) {
    if (q == "0.5") {
      lty=1
    } else {
      lty=2
    }
    if (var.name == "CO2") {
      if (plot.log) {
        points(seq(2010, 2100, by=5), ipat.quantiles.bycountry[[iso]][[q]][, var.name] / 10^6,
               type="l", lwd=2, lty=lty, col="red")
      } else {
        points(seq(2010, 2100, by=5), ipat.quantiles.bycountry[[iso]][[q]][, var.name] / 10^6,
               type="l", lwd=2, lty=lty, col="red")
      }
    } else if (var.name == "GDPpercapita") {
      if (plot.log) {
        points(seq(2010, 2100, by=5), ipat.quantiles.bycountry[[iso]][[q]][, var.name] / 10^3,
               type="l", lwd=2, lty=lty, col="red")
      } else {
        points(seq(2010, 2100, by=5), ipat.quantiles.bycountry[[iso]][[q]][, var.name] / 10^3,
               type="l", lwd=2, lty=lty, col="red")
      }
    } else {
      if (plot.log) {
        points(seq(2010, 2100, by=5), ipat.quantiles.bycountry[[iso]][[q]][, var.name],
               type="l", lwd=2, lty=lty, col="red")
      } else {
        points(seq(2010, 2100, by=5), ipat.quantiles.bycountry[[iso]][[q]][, var.name],
               type="l", lwd=2, lty=lty, col="red")
      }
    }
  }
}


plot.5region.co2 <- function(co2.5regions, co2.5regions.quants,
                             save.plots=F) {
  # plot.5region.co2 plots projections of CO2 emissions for the 5 IPCC regions
  if (save.plots) pdf(paste0(plot.location, "co2_regions.pdf"),
                      width=7.20, height=9.72)
  ggplot.list <- list()
  num_screen <- 1
  for (reg in rcp.5region.names.altorder) {
    rcp.5region.tmp <- subset(rcp.5region, Region == reg)
    y.min <- min(0, min(rcp.5region.tmp[, -(1:3)])*1.04)
    y.max <- max(max(co2.5regions.quants[[reg]][,-1] / 10^9), 
                 max(rcp.5region.tmp[, -(1:3)])) * 1.04
    plot.obj <- ggplot() +
      labs(color="Scenario") +
      theme(panel.background=element_rect(fill="white"),
            panel.grid.minor=element_line(color="white"),
            panel.grid.major=element_line(color="grey90"),
            text=element_text(family=font, size=15)) +
      scale_x_continuous(breaks=seq(1960,2100,by=20))
    row.ind <- which(olddata.5region$Region == reg)
    # Plot historical data
    region.toplot <- t(as.matrix(co2.5regions.quants[[reg]][,-1])/1e9)
    region.toplot <- as.data.frame(region.toplot)
    names(region.toplot) <- paste0('Quantiles', c(0.025, 0.05, 0.5, 0.95, 0.975))
    region.toplot$Year <- seq(2010, 2100, by=5)
    #browser()
    plot.obj <- plot.obj + 
      ylab('Total CO2 Emissions (gt)') + 
      geom_ribbon(aes(x=Year, ymin=Quantiles0.05, ymax=Quantiles0.95), data=region.toplot,
                  alpha=0.3, fill='blue') +
      geom_ribbon(aes(x=Year, ymin=Quantiles0.025, ymax=Quantiles0.975), data=region.toplot,
                  alpha=0.2, fill='blue') +
      geom_line(data=region.toplot, aes(x=Year, y=Quantiles0.5,
                                        color='Projection'),
                alpha=1, size=1.3, color='blue')
    # browser()
    historical.data <- data.frame(Year = 1960:2010, CO2 = as.numeric(olddata.5region[row.ind, -1] / 10^9))
    plot.obj <- plot.obj + geom_line(aes(Year, CO2), data = historical.data, size = 1)
    plot.obj <- plot.obj + ggtitle(paste('Region', reg))
    # Plot IPCC RCPs
    rcp.col.ind <- 0
    for (scenario in c("RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5")) {
      rcp.col.ind <- rcp.col.ind + 1
      rcp.col <- rcp.colors[rcp.col.ind]
      row.ind <- which(rcp.5region.tmp$Scenario == scenario)
      years.tmp <- c(2000, 2005, seq(2010, 2100, by=10))
      vars.tmp <- paste0("Carbon", years.tmp)
      data.scenario <- data.frame(Year = years.tmp, CO2 = as.numeric(rcp.5region.tmp[row.ind, vars.tmp]))
      plot.obj <- plot.obj + geom_line(data = data.scenario, aes(Year, CO2), color = rcp.col, size = 1.3)
    }
    
    ggplot.list[[num_screen]] <- plot.obj
    num_screen <- num_screen + 1
  }
  if (save.plots) 
  {
    grid.arrange(ggplot.list[[1]], ggplot.list[[2]], ggplot.list[[3]], 
                 ggplot.list[[4]], ggplot.list[[5]], ncol = 2)
    dev.off()
  }
  return(ggplot.list)
}

plot.5region.co2.back <- function(co2.5regions, co2.5regions.quants,
                                  year.start, year.end,
                                  save.plots=F, screen.index) {
  # plot.5region.co2.back plots projections of CO2 emissions in the 5 IPCC regions
  # for the out of sample validation work
  
  if (save.plots) {
    pdf(paste0(plot.location, "co2_regions_", year.start, ".pdf"),
        height=10, width=7)
    par(mfrow=c(3,2))
  }
  num_screen <- 0
  ggplot.list <- list()
  for (reg in rcp.5region.names.altorder) {
    if (year.start == 1980 & reg == "REF") next
    y.min <- 0
    year.seq <- seq(year.start, year.end, by=5)
    row.ind <- which(olddata.5region$Region == reg)
    col.inds <- 1:length(year.seq)
    y.max <- max(olddata.5region[row.ind,-1], 
                 max(co2.5regions.quants[[reg]][, 1+col.inds])) / 10^9 * 1.04
    screen(screen.index[num_screen])
    num_screen <- num_screen + 1
    plot.obj <- ggplot() +
      labs(color="Scenario") +
      theme(panel.background=element_rect(fill="white"),
            panel.grid.minor=element_line(color="white"),
            panel.grid.major=element_line(color="grey90"),
            text=element_text(family=font, size=15)) +
      scale_x_continuous(breaks=seq(1960,2100,by=15))
    
    region.toplot <- t(as.matrix(co2.5regions.quants[[reg]][,-1])/1e9)
    region.toplot <- as.data.frame(region.toplot)
    names(region.toplot) <- paste0('Quantiles', c(0.025, 0.05, 0.5, 0.95, 0.975))
    region.toplot$Year <- seq(year.start, year.end, by=5)
    #browser()
    plot.obj <- plot.obj + 
      ylab('Total CO2 Emissions (gt)') + 
      geom_ribbon(aes(x=Year, ymin=Quantiles0.05, ymax=Quantiles0.95), data=region.toplot,
                  alpha=0.3, fill='blue') +
      geom_ribbon(aes(x=Year, ymin=Quantiles0.025, ymax=Quantiles0.975), data=region.toplot,
                  alpha=0.2, fill='blue') +
      geom_line(data=region.toplot, aes(x=Year, y=Quantiles0.5,
                                        color='Projection'),
                alpha=1, size=1.3, color='blue')
    #browser()
    # Plot historical data
    historical.data <- data.frame(Year = 1960:2010, CO2 = as.numeric(olddata.5region[row.ind, -1] / 10^9))
    plot.obj <- plot.obj + geom_line(aes(Year, CO2), data = historical.data)
    plot.obj <- plot.obj + ggtitle(paste('Region', reg))
    ggplot.list[[num_screen]] <- plot.obj
  }
  if (save.plots) {
    dev.off()
    # set par back to normal
    par(mfrow=c(1,1))
  }
  return(ggplot.list)
}



# =================================

# for supplementary information

get.posterior.table <- function(corr.model.out.params,
                                names.countries.intensity,
                                names.countries.gdp) {
  # get.posterior.table outputs a table of posterior distributions of different
  # model parameters
  
  table.formatted.tmp <- summary(corr.model.out.params,
                                 quantiles=c(0.05, 0.5, 0.95))[[2]]
  stopifnot(! "USA" %in% names.countries.intensity)
  stopifnot(! "USA" %in% names.countries.gdp)
  names.countries.intensity <- getCountry(names.countries.intensity)
  names.countries.gdp <- getCountry(names.countries.gdp)
  rownames.tmp <- rownames(table.formatted.tmp)
  table.formatted <- table.formatted.tmp
  for (var.name in c("delta[", "sig.eps[")) {
    rows.tmp <- which(substr(rownames.tmp, 1, nchar(var.name)) == var.name)
    rownames(table.formatted)[rows.tmp] <- paste0(var.name, names.countries.intensity, "]")
  }
  for (var.name in c("phi[", "sig.gap[")) {
    rows.tmp <- which(substr(rownames.tmp, 1, nchar(var.name)) == var.name)
    rownames(table.formatted)[rows.tmp] <- paste0(var.name, names.countries.gdp, "]")
  }
  
  return(table.formatted)
}

params.posterior.formatted <- get.posterior.table(model.ar1.const.results.2010$corr.model.out[[2]],
                                                  model.ar1.const.results.2010$names.countries.intensity,
                                                  model.ar1.const.results.2010$corr.model.out[[1]]$names.countries.gdp)
xtable(params.posterior.formatted, digits=4)



params.posterior.table <- summary(model.ar1.const.results.2010$corr.model.out[[2]],
                                  quantiles=c(0.05, 0.5, 0.95))
# parameter statistics
params.posterior.table[[1]]

# parameter quantiles
params.posterior.table[[2]]

library(xtable)
xtable(params.posterior.table[[2]], digits=4)




# =========================================



make.table.true.tests <- function(prop.true.tests.ipat, which.quant) {
  # make.table.true.tests creates tables of out of sample validation checks
  # by IPAT component
  library(xtable)
  dim.wide <- length(names(prop.true.tests.ipat)) + 1
  dim.long <- length(names(prop.true.tests.ipat[[1]]))
  table.true.tests <- data.frame(matrix(rep(NA, dim.wide*dim.long), ncol=dim.wide))
  names(table.true.tests) <- c("Year", names(prop.true.tests.ipat))
  table.true.tests$Year <- names(prop.true.tests.ipat[[1]])
  for (var.name in names(prop.true.tests.ipat)) {
    row.index <- 0
    for (year in names(prop.true.tests.ipat[[var.name]])) {
      row.index <- row.index + 1
      table.true.tests[row.index, var.name] <- 
        prop.true.tests.ipat[[var.name]][[year]][which.quant]
    }
  }
  
  print(xtable(table.true.tests, digits=4), include.rownames=F)
  return(table.true.tests)
}



# =========================================


# Plotting ipat components
# Based on similar code in estimate_model_corr.R
for (var.name in c("CO2", "GDP", "GDPpercapita", "Tech", "Pop")) {
  for (plot.log in c(TRUE, FALSE)) {
    if (plot.log & (var.name == "Pop")) next
    i.iso <- 0
    for (iso in unique(data.medium$Isocode)) {
      if (i.iso %% 8 == 0) {
        if (i.iso != 0) dev.off()
        pdf(paste0(plot.location, var.name, "_log", plot.log, "_", i.iso/8, ".pdf"),
            height=800, width=500)
        par(mfrow=c(4,2))
      }
      i.iso <- i.iso + 1
      summarize.ipat.country(data.medium,
                             proj.evals.2010.ar1.const$ipat.components.bycountry,
                             proj.evals.2010.ar1.const$ipat.quantiles.bycountry,
                             iso, var.name, plot.log=plot.log)
    }
    dev.off()
  }
}


# Out of sample validation tables
make.table.true.tests(proj.evals.2000.ar1.const$prop.belowmedian, 1)
make.table.true.tests(proj.evals.1990.ar1.const$prop.belowmedian, 1)
make.table.true.tests(proj.evals.1980.ar1.const$prop.belowmedian, 1)


make.table.true.tests(proj.evals.2000.ar1.const$prop.true.tests.ipat, 1)
make.table.true.tests(proj.evals.1990.ar1.const$prop.true.tests.ipat, 1)
make.table.true.tests(proj.evals.1980.ar1.const$prop.true.tests.ipat, 1)


make.table.true.tests
# Plotting for 5 IPCC regions
pdf(paste0(plot.location, "co2_regions", ".pdf"),
    height=15, width=10)

gglist.5region <- plot.5region.co2(proj.evals.2010.ar1.const$co2.5regions,
                                   proj.evals.2010.ar1.const$co2.5regions.quants, save.plots=F)
grid.arrange(gglist.5region[[1]], gglist.5region[[2]], gglist.5region[[3]],
             gglist.5region[[4]], gglist.5region[[5]], ncol = 2)
dev.off()

pdf(paste0(plot.location, "co2_regions_all", ".pdf"),
    height=15, width=9)
gglist1 <- plot.5region.co2.back(proj.evals.2000.ar1.const$co2.5regions,
                                 proj.evals.2000.ar1.const$co2.5regions.quants,
                                 year.start=2000, year.end=2010,
                                 save.plots=F, screen.index = seq(1,15,3))

gglist2 <- plot.5region.co2.back(proj.evals.1990.ar1.const$co2.5regions,
                                 proj.evals.1990.ar1.const$co2.5regions.quants,
                                 year.start=1990, year.end=2010,
                                 save.plots=F, screen.index = seq(1,15,3))

gglist3 <- plot.5region.co2.back(proj.evals.1980.ar1.const$co2.5regions,
                                 proj.evals.1980.ar1.const$co2.5regions.quants,
                                 year.start=1980, year.end=2010,
                                 save.plots=F, screen.index = seq(1,15,3))
grid.arrange(gglist3[[1]], gglist2[[1]], gglist1[[1]],
             gglist3[[2]], gglist2[[2]], gglist1[[2]],
             gglist3[[3]], gglist2[[3]], gglist1[[3]],
             gglist3[[4]], gglist2[[4]], gglist1[[4]],
             blankPanel, gglist2[[5]], gglist1[[5]],
             ncol=3)
dev.off()

plot.5region.co2.back(proj.evals.2000.ar1.const$co2.5regions,
                      proj.evals.2000.ar1.const$co2.5regions.quants,
                      year.start=2000, year.end=2010,
                      save.plots=T)

plot.5region.co2.back(proj.evals.1990.ar1.const$co2.5regions,
                      proj.evals.1990.ar1.const$co2.5regions.quants,
                      year.start=1990, year.end=2010,
                      save.plots=T)

plot.5region.co2.back(proj.evals.1980.ar1.const$co2.5regions,
                      proj.evals.1980.ar1.const$co2.5regions.quants,
                      year.start=1980, year.end=2010,
                      save.plots=T)


##################################################################
#  Analysis based on Figure SPM.5 of AR5 Summary Report	 #
##################################################################
# see calc_temp_changes.R

## 3/9/2016: Analysis without upper bound
# Find spline approximations to upper and lower bounds of temperature curve
# cum.emissions.2100 is an array of cumulative emissions in 2100, for the different trajectories
cum.emissions.2100 <- proj.evals.2010.ar1.const$cum.emissions[,19]
#load ("cum_emissions_2100.Rda")		# Posterior distribution of 2100 cumCO2, gives cum.emissions.2100
spm5 <- read.csv(paste0(data.location, "ipcc_ar5_temp_FigSPM5.csv"))
spm5[10,1] <- spm5[10,1] * 1000		# That one value was in different units
cumco2.spm5 <- spm5[,1]
templower.spm5 <- spm5[,2]
tempupper.spm5 <- spm5[1:7,3]
templower.splinefun <- splinefun (x=cumco2.spm5, y=templower.spm5, 
                                  method="hyman")
tempupper.splinefun <- splinefun (x=cumco2.spm5[1:7], y=tempupper.spm5,
                                  method="hyman")

# Plot the spline functions
co2 <- seq (1010,8950,by=100)
templower <- templower.splinefun (co2)
tempupper <- tempupper.splinefun (co2)
tempmean <- (templower+tempupper)/2
pdf(file=paste0(plot.location, "SPM5rangeExtrapolated.pdf"))
#pdf (file="SPM5rangeExtrapolated.pdf")
matplot (co2,cbind(templower,tempmean,tempupper),lty=c(2,1,2),type="l",
         col=c("red","black","red"), main="Figure SPM.5 extrapolated",
         xlab="Cumulative CO2 emissions", ylab="Temperature increase to 2100" )
dev.off ()

# Find posterior distribution of temperature increase
# Note: posterior distribution is in gt 2010-2100
#  Figure SPM.5 is in cumulative CO2 emissions since 1861-1880
# From AR5 summary, p.10, emissions to 2011 were about 1900.

nsim <- length(cum.emissions.2100)
cumco2.post <- cum.emissions.2100 * (10^(-9)) + 1900 # Convert to gt 1870-2100.
tempinc.spm5 <- rep (NA, nsim)
n.sims.pertraj <- 100
for (k in 1:nsim) {
  co2 <- cumco2.post[k]
  co2 <- min (co2, 8950)  # 8950 is the upper bound of cum C02 in Figure SPM.5.
  templower <- templower.splinefun(co2)
  tempupper <- tempupper.splinefun(co2)
  tempmean <- (templower + tempupper)/2
  tempsd <- (tempupper - templower)/3.28
  temp <- rnorm (n.sims.pertraj, tempmean, tempsd)
  ind.init <- (k-1)*n.sims.pertraj + 1
  tempinc.spm5[ind.init:(ind.init+n.sims.pertraj-1)] <- temp
}
hist(tempinc.spm5)


round (quantile (tempinc.spm5, c(.025, .05, .1, .5, .9, .95, .975)), 1)
length (tempinc.spm5[tempinc.spm5 <= 2.0]) / nsim
length (tempinc.spm5[tempinc.spm5 <= 1.5]) / nsim
round (quantile (cumco2.post, c(.025, .05, .1, .5, .9, .95, .975)))
round (quantile (cumco2.post-1900, c(.025, .05, .1, .5, .9, .95, .975)))

pdf(file = paste0(plot.location, "tempinc.SPM5.pdf"))
#pdf (file = "tempinc.SPM5.pdf")
hist (tempinc.spm5, xlab="Temperature Increase to 2100 (deg C)", 
      ylab="Frequency (out of 1000)", main="")
dev.off ()

color.hist <- "#e41a1c"
library(scales) # My older version of ggplot2 doesn't have alpha()
temp.inc.plot <- ggplot(tempinc.spm5.data, aes(x=tempinc.spm5)) + 
  geom_histogram(alpha=0.5, binwidth=0.25, fill=color.hist, color=alpha(color.hist, 1)) +
  ylab("Count (out of 100,000)") +
  xlab(expression(paste("Temperature change relative to 1861-1880 (", ~degree~C, ")"))) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.minor=element_line(color="white"),
        panel.grid.major=element_line(color="grey90")) +
  scale_x_continuous(limits=c(0, 8), breaks=0:8)
pdf(paste0(plot.location, "tempinc_SPM5.pdf"))
print(temp.inc.plot)
dev.off()



tempinc.spm5.data <- data.frame(tempinc.spm5)
ggplot(tempinc.spm5.data, aes(x=tempinc.spm5)) + geom_bar()

color.hist <- "#e41a1c"
library(scales) # My older version of ggplot2 doesn't have alpha()
library(ggplot2)
temp.inc.plot <- ggplot(tempinc.spm5.data, aes(x=tempinc.spm5)) + 
  geom_histogram(alpha=0.5, binwidth=0.25, fill=color.hist, color=alpha(color.hist, 1)) +
  ylab("Count (out of 100,000)") +
  xlab(expression(paste('Temperature change relative to 1861-1880 (', degree ~ C, ')'))) +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.minor=element_line(color="white"),
        panel.grid.major=element_line(color="grey90")) +
  scale_x_continuous(limits=c(0, 8), breaks=0:8)
pdf(paste0(plot.location, "tempinc_SPM5.pdf"))
print(temp.inc.plot)
dev.off()



#=======================================
# Probability of INDC targets being achieved under our model
#=======================================
#========================
# China:
# Emissions peaking by 2030

china.co2.max.years <- sapply(1:n.trajectories, function(i.tmp) {
  which.max(proj.evals.2010.ar1.const$ipat.components.bycountry[[i.tmp]][["CHN"]][, "CO2"])
})
ind.2030.tmp <- which(proj.evals.2010.ar1.const$ipat.components.bycountry[[1]][["CHN"]][, "Year"] == 2030)
# 25% of the trajectories have China's peak CO2 emissions by 2030
mean(china.co2.max.years <= ind.2030.tmp)
#[1] 0.251
# Note that 5% of these peak emissions for China occur in 2030
# (our projections have 5-year granularity)
mean(china.co2.max.years == ind.2030.tmp)
#[1] 0.049



# Carbon intensity being reduced by 60-65% by 2030 relative to 2005 levels
tech.china.2005.tmp <- subset(data.medium, Year == 2005 & Isocode == "CHN")$Tech
is.china.meeting.reduction <- sapply(1:n.trajectories, function(i.tmp) {
  china.2030.tech.val.tmp <- proj.evals.2010.ar1.const$ipat.components.bycountry[[i.tmp]][["CHN"]][ind.2030.tmp, "Tech"]
  return(china.2030.tech.val.tmp < (1 - 0.6) * tech.china.2005.tmp)
})
# China meets its weak intensity goal of reducing by 60% in 12% of cases
mean(is.china.meeting.reduction)
#[1] 0.116
is.china.meeting.reduction.strong <- sapply(1:n.trajectories, function(i.tmp) {
  china.2030.tech.val.tmp <- proj.evals.2010.ar1.const$ipat.components.bycountry[[i.tmp]][["CHN"]][ind.2030.tmp, "Tech"]
  return(china.2030.tech.val.tmp < (1 - 0.65) * tech.china.2005.tmp)
})
# China meets its strong intensity goal of reducing by 65% in 5% of cases
mean(is.china.meeting.reduction.strong)
#[1] 0.054

plot.paristarget <- function(Iso, objective, indeces, intensity = F, newdata = NULL, title)
{
  data.historical <- subset(data.medium, Isocode == Iso)
  
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=15)) +
    scale_x_continuous(breaks=seq(2010,2100,by=20))
  
  data.plot <- data.frame(Year = proj.evals.2010.ar1.const$ipat.quantiles.bycountry[[Iso]][["0.025"]]$Year[indeces])
  quantiles.need <- c(0.025, 0.05, 0.5, 0.95, 0.975)
  obj.data <- data.frame(Year = objective[1], obj = objective[2])
  if(intensity)
  {
    for(quant in quantiles.need)
    {
      data.plot[[paste0('quantiles', quant)]] <- proj.evals.2010.ar1.const$ipat.quantiles.bycountry[[Iso]][[as.character(quant)]]$Tech[indeces]
    }
    plot.obj <- plot.obj + 
      ylab('tonnes of CO2 per $10,000') + 
      geom_ribbon(aes(x=Year, ymin=quantiles0.05, ymax=quantiles0.95), data=data.plot,
                  alpha=0.3, fill=color.us) +
      geom_ribbon(aes(x=Year, ymin=quantiles0.025, ymax=quantiles0.975), data=data.plot,
                  alpha=0.2, fill=color.us) +
      geom_line(data=data.plot, aes(x=Year, y=quantiles0.5,
                                    color="Projections"),
                alpha=1, size=1.3, color=color.us)
    plot.obj <- plot.obj + geom_line(data = data.historical, aes(x=Year, y = Tech), size = 1.2) + 
      geom_point(aes(x = Year, y = obj), data = obj.data, color="Blue", size = 4)
  }
  else
  {
    data.historical$CO2Total <- data.historical$CO2 * (1000 * data.historical$PopTotal)/1e9
    data.2015 <- data.frame(Year = 2015, CO2 = newdata)
    for(quant in quantiles.need)
    {
      data.plot[[paste0('quantiles', quant)]] <- proj.evals.2010.ar1.const$ipat.quantiles.bycountry[[Iso]][[as.character(quant)]]$CO2[indeces]/1e9
    }
    plot.obj <- plot.obj + 
      ylab('Total CO2 Emissions (gt)') + 
      geom_ribbon(aes(x=Year, ymin=quantiles0.05, ymax=quantiles0.95), data=data.plot,
                  alpha=0.3, fill=color.us) +
      geom_ribbon(aes(x=Year, ymin=quantiles0.025, ymax=quantiles0.975), data=data.plot,
                  alpha=0.2, fill=color.us) +
      geom_line(data=data.plot, aes(x=Year, y=quantiles0.5,
                                    color="Projections"),
                alpha=1, size=1.3, color=color.us)
    plot.obj <- plot.obj + geom_line(data = data.historical, aes(x=Year, y = CO2Total), size = 1.2) + 
      geom_point(aes(x = Year, y = obj), data = obj.data, color="Blue", size = 4) + 
      geom_point(aes(x = Year, y = CO2), data = data.2015,  color="black", size = 4)
  }
  plot.obj <- plot.obj + ggtitle(title)
  return (plot.obj)
}


pdf(paste0(plot.location, "Responses/ParisTarget_China.pdf"), height = 5, width = 5)
plot.china <- plot.paristarget('CHN', objective = c(2030, (1 - 0.6) * tech.china.2005.tmp), 
                               indeces = 1:5, intensity = T, title = 'China Intensity')
print(plot.china)
dev.off()


#========================
# United States:
# Reduction of emissions by 26-28% by 2025 over 2005 levels

ind.2025.tmp <- which(proj.evals.2010.ar1.const$ipat.components.bycountry[[1]][["USA"]][, "Year"] == 2025)

base.usa.2005.tmp <- subset(data.medium, Year == 2005 & Isocode == "USA")$CO2 *
  (1000 * subset(data.medium, Year == 2005 & Isocode == "USA")$PopTotal)
is.usa.meeting.reduction <- sapply(1:n.trajectories, function(i.tmp) {
  usa.2025.val.tmp <- proj.evals.2010.ar1.const$ipat.components.bycountry[[i.tmp]][["USA"]][ind.2025.tmp, "CO2"]
  return(usa.2025.val.tmp < (1-.26) * base.usa.2005.tmp)
  #  return(usa.2025.val.tmp <  base.usa.2005.tmp)
})
# The USA meets its commitment of 26% in 0.7% of our trajectories
mean(is.usa.meeting.reduction)
#[1] 0.007
is.usa.meeting.reduction.strong <- sapply(1:n.trajectories, function(i.tmp) {
  usa.2025.val.tmp <- proj.evals.2010.ar1.const$ipat.components.bycountry[[i.tmp]][["USA"]][ind.2025.tmp, "CO2"]
  return(usa.2025.val.tmp < (1-.28) * base.usa.2005.tmp)
  #  return(usa.2025.val.tmp <  base.usa.2005.tmp)
})
# The USA meets its commitment of 28% in 0.7% of our trajectories
mean(is.usa.meeting.reduction.strong)
#[1] 0.003

data.us.historical <- subset(data.medium, Isocode == "USA")
data.us.historical$CO2Total <- data.us.historical$CO2 * (1000 * data.us.historical$PopTotal)

pdf(paste0(plot.location, "Responses/ParisTarget_USA.pdf"), height = 5, width = 5)
plot.USA <- plot.paristarget('USA', c(2025, (1-.28) * base.usa.2005.tmp / 10^9), indeces = 1:4, 
                             intensity = F, newdata = 1.38997 * (11 / 3), title = 'United States Total CO2 Emissions')
print(plot.USA)
dev.off()


#========================
# EU (European Union)

#CHANGE
#n.trajectories <- 1000
n.trajectories <- length(proj.evals.2010.ar1.const$ipat.components.bycountry)
# First, find emissions of the EU in 1990
data.eu.tmp <- subset(data.medium, Isocode %in% isos.eu & Year == 1990)
eu.1990.co2.tmp <- sum(data.eu.tmp$CO2 * (1000 * data.eu.tmp$PopTotal))

# Now, look at projections for the EU

# we don't have Cyprus, Luxembourg, or Malta in our dataset
isos.eu.indata <- isos.eu[isos.eu %in% unique(data.medium$Isocode)]
is.eu.meeting.reduction <- sapply(1:n.trajectories, function(i.tmp) {
  co2.eu.total.2030 <- 0
  for (iso.tmp in isos.eu.indata) {
    co2.eu.total.2030 <- co2.eu.total.2030 +
      proj.evals.2010.ar1.const$ipat.components.bycountry[[i.tmp]][[iso.tmp]][ind.2030.tmp, "CO2"]
  }
  return(co2.eu.total.2030 < (1 - .4) * eu.1990.co2.tmp)
})
# The EU is meeting its reduction goals in none of our scenarios
mean(is.eu.meeting.reduction)
#[1] 0

# Plot these projections
get.eu.projections.2030 <- function(n.trajectories,
                                    proj.evals.2010.ar1.const,
                                    quantiles=c(0.025,0.05,0.5,0.95,0.975)) {
  #sapply(1:n.trajectories, function(i.tmp) {
  years.tmp <- seq(2010, 2030, by=5)
  co2.eu.total <- matrix(0, nrow=n.trajectories, ncol=length(years.tmp),
                         dimnames=list(c(),paste0("CO2", years.tmp)))
  
  year.ind <- 0
  for (year.tmp in years.tmp) {
    year.ind <- year.ind + 1
    # Arbitrary country, GBR, used below
    ind.year.tmp <- which(proj.evals.2010.ar1.const$ipat.components.bycountry[[1]][["GBR"]][, "Year"] == year.tmp)
    year.tmp.chr <- as.character(year.tmp)
    for (i.tmp in 1:n.trajectories) {
      #co2.eu.total[[year.tmp.chr]] <- 0
      for (iso.tmp in isos.eu.indata) {
        co2.eu.total[i.tmp, year.ind] <- co2.eu.total[i.tmp, year.ind] +
          proj.evals.2010.ar1.const$ipat.components.bycountry[[i.tmp]][[iso.tmp]][ind.year.tmp, "CO2"]
      }
    }
  }
  #return(co2.eu.total.2030 < (1 - .4) * eu.1990.co2.tmp)
  
  co2.eu.total.quant <- matrix(NA, nrow=length(quantiles), ncol=length(years.tmp),
                               dimnames=list(quantiles, paste0("CO2", years.tmp)))
  for (i.tmp in 1:length(years.tmp)) {
    co2.eu.total.quant[,i.tmp] <- quantile(co2.eu.total[,i.tmp], probs=quantiles)
  }
  
  #  return(co2.eu.total)
  return(co2.eu.total.quant)
}

eu.co2.project.quant <- get.eu.projections.2030(n.trajectories, proj.evals.2010.ar1.const)
# TODO: isos.eu, get historic emissions data
data.eu.historical <- subset(data.medium, Isocode %in% isos.eu)
data.eu.historical$CO2Total <- data.eu.historical$CO2 * (1000 * data.eu.historical$PopTotal)
library(doBy)
#summaryBy
data.eu.comb.historical <- summaryBy(CO2Total ~ Year, data=data.eu.historical, FUN=sum)



# multiply the below by (11/3)
# Austria 16.7983776507337
# Belgium 25.1137773007619
# Bulgaria 12.1354873237362
# Croatia 4.87836248509944

# Cyprus 1.64287078513157
# Czech Republic 26.2551319120317
# Denmark 8.8025108485609
# Estonia 5.36021545311619

# Finland 10.5338276676075
# France 83.2041991152356
# Germany 196.198319733223
# Greece 17.2949265206269
# Hungary 11.7966485607259

# Ireland 9.61241258402785
# Italy 91.6478739196349
# Latvia 1.95672503810325
# Lithuania 3.31001835842373

# Luxembourg 2.81076692352811
# Malta 0.613791910302059
# Netherlands 44.7276131711787
# Poland 79.538115826858

# Portugal 13.3771820605989
# Romania 19.7011736114519
# Slovakia 8.56222094873605
# Slovenia 3.94916932957372

# Spain 67.5899261282673
# Sweden 11.7108871218644
# United Kingdom 108.467020521798

data.eu.2015 <- (11/3) * (16.7983776507337 + 25.1137773007619 + 12.1354873237362 + 
                            4.87836248509944 + 1.64287078513157 + 26.2551319120317 +
                            8.8025108485609 + 5.36021545311619 + 10.5338276676075 + 
                            83.2041991152356 + 196.198319733223 + 17.2949265206269 + 
                            11.7966485607259 + 9.61241258402785 + 91.6478739196349 + 
                            1.95672503810325 + 3.31001835842373 + 2.81076692352811 +
                            0.613791910302059 + 44.7276131711787 + 79.538115826858 +
                            13.3771820605989 + 19.7011736114519 + 8.56222094873605 +
                            3.94916932957372 + 67.5899261282673 + 11.7108871218644 + 
                            108.467020521798) / 1000
pdf(paste0(plot.location, "Responses/ParisTarget_EU.pdf"), height = 5, width = 5)

plot.eu <- ggplot() +
  labs(color="Scenario") +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.minor=element_line(color="white"),
        panel.grid.major=element_line(color="grey90"),
        text=element_text(family=font, size=15)) +
  scale_x_continuous(breaks=seq(2010,2100,by=20))
data.plot <- data.frame(t(eu.co2.project.quant))/1e9
quantiles.need <- c(0.025, 0.05, 0.5, 0.95, 0.975)
data.plot$Year <- seq(2010, 2030, 5)
names(data.plot)[1:5] <- paste0('quantiles', quantiles.need)
data.2015 <- data.frame(Year = 2015, CO2 = data.eu.2015)
obj.data <- data.frame(Year = 2030, obj = (1 - .4) * eu.1990.co2.tmp / 10^9)
plot.eu <- plot.eu + 
  ylab('Total CO2 Emissions (gt)') + 
  geom_ribbon(aes(x=Year, ymin=quantiles0.05, ymax=quantiles0.95), data=data.plot,
              alpha=0.3, fill=color.us) +
  geom_ribbon(aes(x=Year, ymin=quantiles0.025, ymax=quantiles0.975), data=data.plot,
              alpha=0.2, fill=color.us) +
  geom_line(data=data.plot, aes(x=Year, y=quantiles0.5,
                                color="Projections"),
            alpha=1, size=1.3, color=color.us)
plot.eu <- plot.eu + geom_line(data = data.eu.comb.historical, aes(x=Year, y = CO2Total.sum/1e9), size = 1.2) + 
  geom_point(aes(x = Year, y = obj), data = obj.data, color="Blue", size = 4) + 
  geom_point(aes(x = Year, y = CO2), data = data.2015,  color="black", size = 4) + 
  ggtitle("EU Total CO2 Emissions")
print(plot.eu)

dev.off()

#========================
# India:
# Reduction of carbon intensity by 33-35% over 2005 levels by 2030

names(proj.evals.2010.ar1.const)
names(proj.evals.2010.ar1.const$ipat.components.bycountry[[1]][["CHN"]])

# Carbon intensity being reduced by 33-35% by 2030 relative to 2005 levels
tech.india.2005.tmp <- subset(data.medium, Year == 2005 & Isocode == "IND")$Tech
is.india.meeting.reduction <- sapply(1:n.trajectories, function(i.tmp) {
  india.2030.tech.val.tmp <- proj.evals.2010.ar1.const$ipat.components.bycountry[[i.tmp]][["IND"]][ind.2030.tmp, "Tech"]
  return(india.2030.tech.val.tmp < (1 - 0.33) * tech.india.2005.tmp)
})
# India meets its weak intensity goal of reducing by 33% in 54% of cases
mean(is.india.meeting.reduction)
#[1] 0.539
is.india.meeting.reduction <- sapply(1:n.trajectories, function(i.tmp) {
  india.2030.tech.val.tmp <- proj.evals.2010.ar1.const$ipat.components.bycountry[[i.tmp]][["IND"]][ind.2030.tmp, "Tech"]
  return(india.2030.tech.val.tmp < (1 - 0.35) * tech.india.2005.tmp)
})
# India meets its strong intensity goal of reducing by 35% in 46% of cases
mean(is.india.meeting.reduction)
#[1] 0.458

data.india.historical <- subset(data.medium, Isocode == "IND")
pdf(paste0(plot.location, "Responses/ParisTarget_India.pdf"), height = 5, width = 5)
plot.India <- plot.paristarget('IND', c(2030, (1 - 0.35) * tech.india.2005.tmp), indeces = 1:5, 
                               intensity = T, title = 'India Intensity')
print(plot.India)
dev.off()


#========================
# Russia:
# Reduction of emissions by 25-30% over 1990 levels by 2030


base.russia.1990.tmp <- subset(data.medium, Year == 1990 & Isocode == "RUS")$CO2 *
  (1000 * subset(data.medium, Year == 1990 & Isocode == "RUS")$PopTotal)
is.russia.meeting.reduction <- sapply(1:n.trajectories, function(i.tmp) {
  russia.2030.val.tmp <- proj.evals.2010.ar1.const$ipat.components.bycountry[[i.tmp]][["RUS"]][ind.2030.tmp, "CO2"]
  return(russia.2030.val.tmp < (1-.25) * base.russia.1990.tmp)
})
# Russia meets its commitment of 25% in 78% of our trajectories
mean(is.russia.meeting.reduction)
#[1] 0.777

is.russia.meeting.reduction.strong <- sapply(1:n.trajectories, function(i.tmp) {
  russia.2030.val.tmp <- proj.evals.2010.ar1.const$ipat.components.bycountry[[i.tmp]][["RUS"]][ind.2030.tmp, "CO2"]
  return(russia.2030.val.tmp < (1-.3) * base.russia.1990.tmp)
})
# Russia meets its commitment of 30% in 73% of our trajectories
mean(is.russia.meeting.reduction.strong)
#[1] 0.728

data.russia.historical <- subset(data.medium, Isocode == "RUS")
data.russia.historical$CO2Total <- data.russia.historical$CO2 * (1000 * data.russia.historical$PopTotal)

pdf(paste0(plot.location, "Responses/ParisTarget_Russia.pdf"), height = 5, width = 5)
plot.rus <- plot.paristarget(Iso = 'RUS', c(2030, (1-.3) * base.russia.1990.tmp / 10^9), indeces = 1:5, 
                             intensity = F, newdata = 471.432845976477 * (11 / 3000), title = "Russia Total CO2 Emissions")
print(plot.rus)
dev.off()



#========================
# Japan:
# Reduction of emissions by 25.4% by 2030 over 2005 levels


base.japan.2005.tmp <- subset(data.medium, Year == 2005 & Isocode == "JPN")$CO2 *
  (1000 * subset(data.medium, Year == 2005 & Isocode == "JPN")$PopTotal)
is.japan.meeting.reduction <- sapply(1:n.trajectories, function(i.tmp) {
  japan.2030.val.tmp <- proj.evals.2010.ar1.const$ipat.components.bycountry[[i.tmp]][["JPN"]][ind.2030.tmp, "CO2"]
  return(japan.2030.val.tmp < (1-.254) * base.japan.2005.tmp)
})
# Japan meets its commitment of 25.4% in 0.7% of our trajectories
mean(is.japan.meeting.reduction)
#[1] 0.058


data.japan.historical <- subset(data.medium, Isocode == "JPN")
data.japan.historical$CO2Total <- data.japan.historical$CO2 * (1000 * data.japan.historical$PopTotal)

pdf(paste0(plot.location, "Responses/ParisTarget_Japan.pdf"), height = 5, width = 5)
plot.jpn <- plot.paristarget(Iso = 'JPN', c(2030, (1-.254) * base.japan.2005.tmp / 10^9), indeces = 1:5, 
                             intensity = F, newdata =  1.184861, title = "Japan Total CO2 Emissions")
print(plot.jpn)
dev.off()


#========================================================
pdf(paste0(plot.location, "Responses/ParisTarget_All.pdf"), height = 15, width = 10)
grid.arrange(plot.USA, plot.eu, plot.jpn, plot.rus, plot.china, plot.India, ncol = 2)
dev.off()


#========================================================

# Plots for individual countries.
countries.toplot <- c('Brazil', 'China', 'Egypt', 'France', 'India', 'Indonesia', 'Japan', 'South Korea',
                      'Nigeria', 'Russia', 'Saudi Arabia', 'South Africa', 'Thailand', 'United Kingdom', 
                      'United States')

pdf(paste0(plot.location, "GDP_per_capita_by_country", ".pdf"),
    height=14, width=10)
par(mfrow = c(5,3), mar = c(4.5, 4.5, 1.5, 1) + 0.5)

isos.toplot <- sapply(countries.toplot, getIso)
for(iso in isos.toplot)
{
  data.toplot.tmp <- data.medium[data.medium$Isocode == iso, c('Year', 'CO2', 'Tech', 'PopTotal', 'GDP')]
  data.toplot.proj.list <- proj.evals.2010.ar1.const$ipat.quantiles.bycountry[[iso]]
  data.toplot.proj <- data.frame(Year = seq(2010, 2100, 5), q025 = data.toplot.proj.list$`0.025`$GDPpercapita, 
                                 q05 = data.toplot.proj.list$`0.05`$GDPpercapita, q5 = data.toplot.proj.list$`0.5`$GDPpercapita, 
                                 q95 = data.toplot.proj.list$`0.95`$GDPpercapita, q975 = data.toplot.proj.list$`0.975`$GDPpercapita)
  data.toplot.proj[,-1] <- data.toplot.proj[,-1]/1000
  ymin <- min(min(data.toplot.tmp$GDP/1000), min(data.toplot.proj[,-1]), na.rm = T)
  ymax <- max(max(data.toplot.tmp$GDP/1000), max(data.toplot.proj[,-1]), na.rm = T)
  plot(c(1960,2100), c(ymin, ymax), type = 'n', xlab = '', 
       ylab = 'GDP per capita', main = names(isos.toplot)[which(isos.toplot == iso)], log = 'y')
  lines(data.toplot.tmp$Year, data.toplot.tmp$GDP/1000, lwd= 1)
  lines(data.toplot.proj$Year, data.toplot.proj$q5, col = 'red', lty =1, lwd=2)
  lines(data.toplot.proj$Year, data.toplot.proj$q025, col = 'red', lty =2, lwd=2)
  lines(data.toplot.proj$Year, data.toplot.proj$q05, col = 'red', lty =2, lwd=2)
  lines(data.toplot.proj$Year, data.toplot.proj$q95, col = 'red', lty =2, lwd=2)
  lines(data.toplot.proj$Year, data.toplot.proj$q975, col = 'red', lty =2, lwd=2)
  
  
}
dev.off()


temp.inc.toplot <- data.frame(Year = seq(2010, 2100, 5))
q025 <- c()
q050 <- c()
q500 <- c()
q950 <- c()
q975 <- c()
for(year in seq(2010, 2100, 5))
{
  quantiles_tempinc <- quantile(tempinc.spm5.yearly[[as.character(year)]], probs = c(0.025, 0.05, 0.5, 0.95, 0.975))
  q025 <- c(q025, quantiles_tempinc[1])
  q050 <- c(q050, quantiles_tempinc[2])
  q500 <- c(q500, quantiles_tempinc[3])
  q950 <- c(q950, quantiles_tempinc[4])
  q975 <- c(q975, quantiles_tempinc[5])
}
temp.inc.toplot$q025 <- q025
temp.inc.toplot$q050 <- q050
temp.inc.toplot$q500 <- q500
temp.inc.toplot$q950 <- q950
temp.inc.toplot$q975 <- q975

pdf(paste0(plot.location, "Temperature_Change.pdf"))

plot.obj <- ggplot() +
  labs(color="Scenario") +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.minor=element_line(color="white"),
        panel.grid.major=element_line(color="grey90"),
        text=element_text(family=font, size=22)) +
  scale_x_continuous(breaks=seq(2010,2100,by=20))

plot.obj <- plot.obj + 
  ylab(expression(paste('Temperature increase relative to 1861-1880 (', degree ~ C, ')'))) + 
  geom_ribbon(aes(x=Year, ymin=q050, ymax=q950), data=temp.inc.toplot,
              alpha=0.3, fill=color.us) +
  geom_ribbon(aes(x=Year, ymin=q025, ymax=q975), data=temp.inc.toplot,
              alpha=0.2, fill=color.us) +
  geom_line(data=temp.inc.toplot, aes(x=Year, y=q500,
                                      color="Projections"),
            alpha=1, size=1.3, color=color.us)
print(plot.obj)
dev.off()

GBR_emission <- c()
for(i in 1:1000)
{
  GBR_emission <- c(GBR_emission, proj.evals.2010.ar1.const$ipat.components.bycountry[[i]]$GBR$CO2[9]/1e6)
}

year.sequence <- seq(2010, 2100, 5)
trajs.worldwide <- lapply(1:length(model.ar1.const.results.2010$data.proj), function(i) {
  x <- model.ar1.const.results.2010$data.proj[[i]]
  get.ipat.components.total(x$TechData, x$GDPData, x$GDPFrontier, x$TechUSA,
                            model.ar1.const.results.2010$preds.countries.trajs[[i]], 
                            model.ar1.const.results.2010$co2.projections[[i]], 
                            year.sequence=seq(2010, 2100, 5))
})

tech2010 <- c()
tech2050 <- c()
tech2100 <- c()
CO2percapita2010 <- c()
CO2percapita2050 <- c()
CO2percapita2100 <- c()
for(i in 1:1000)
{
  tech2010 <- c(tech2010, trajs.worldwide[[i]]$Tech[1])
  tech2050 <- c(tech2050, trajs.worldwide[[i]]$Tech[9])
  tech2100 <- c(tech2100, trajs.worldwide[[i]]$Tech[19])
  CO2percapita2010 <- c(CO2percapita2010, trajs.worldwide[[i]]$CO2percapita[1])
  CO2percapita2050 <- c(CO2percapita2050, trajs.worldwide[[i]]$CO2percapita[9])
  CO2percapita2100 <- c(CO2percapita2100, trajs.worldwide[[i]]$CO2percapita[19])
}