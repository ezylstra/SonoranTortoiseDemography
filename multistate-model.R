################################################################################
# Multi-state model to estimate demographic rates of Sonoran desert tortoises in 
# Arizona, 1987-2020

# ER Zylstra
# Last updated: 4 April 2024
################################################################################

library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(nimble)
library(MCMCvis)
library(ggplot2)
library(cowplot)
library(popbio)

rm(list=ls())

#------------------------------------------------------------------------------# 
# Load data
#------------------------------------------------------------------------------#

cr <- read.csv("data/CapRecapData.csv", header = TRUE)
plots <- read.csv("data/Plots_NoCoords.csv", header = TRUE)
surveys <- read.csv("data/Surveys.csv", header = TRUE)
pdsi <- read.csv("data/PDSI.csv", header = TRUE)
precip <- read.csv("data/Precip_Monthly.csv", header = TRUE)

#------------------------------------------------------------------------------# 
# Survey data
#------------------------------------------------------------------------------#

# Calculate survey intervals (number of years between consecutive surveys)
  surveys <- surveys %>%
    arrange(code, yr) %>%
    mutate(interval = NA)
  for (i in 2:nrow(surveys)) {
    surveys$interval[i] <- ifelse(surveys$code[i] != surveys$code[i-1], NA,
                                  surveys$yr[i] - surveys$yr[i-1])
  }  
  
  # Summarizing survey intervals
  summary(surveys$interval)

# Identify years when no surveys were done
  missingyrs <- setdiff(1987:2020, sort(unique(surveys$yr)))
  
# Create wide version of the persondays data
  persondays.w <- surveys %>%
    select(-c(interval, plot, area.sqmi)) %>%
    add_row(code = surveys$code[1],
            yr = missingyrs,
            persondays = 0) %>%
    pivot_wider(id_cols = code,
                names_from = yr,
                names_prefix = "y",
                names_sort = TRUE,
                values_from = persondays,
                values_fill = 0) %>%
    data.frame()
  
# Create wide version of survey data, where 1/0 indicates when surveys were done
  surv.w <- surveys %>%
    select(-c(interval, plot, area.sqmi, persondays)) %>%
    mutate(survey = 1) %>%
    add_row(code = surveys$code[1],
            yr = missingyrs,
            survey = 0) %>%    
    pivot_wider(id_cols = code,
                names_from = yr,
                names_prefix = "y",
                names_sort = TRUE,
                values_from = survey,
                values_fill = 0) %>%
    data.frame()  

# Summarize survey effort by plot
  surv.p <- surveys %>%
    group_by(code, plot) %>%
    summarize(area.sqmi = mean(area.sqmi),
              n.surveys = length(persondays),
              first.surv = min(yr),
              last.surv = max(yr),
              .groups = "keep") %>%
    data.frame()
  
  # Average number of surveys per plot
  summary(surv.p$n.surveys)

# Summarize survey effort by year
  surv.y <- surveys %>%
    select(code, yr, persondays) %>%
    add_row(code = surveys$code[1],
            yr = missingyrs,
            persondays = 0) %>%   
    group_by(yr) %>%
    summarize(n.plots = sum(persondays > 0),
              mean.eff = ifelse(n.plots == 0, 
                                NA, 
                                round(sum(persondays)/n.plots, 1))) %>%
    # Identify two "periods" (1987-2008, 2009-2020)
    mutate(period = ifelse(yr %in% 1987:2008, 1, 2)) %>%
    data.frame()
  
  # Average number of plots per year
  summary(surv.y$n.plots)
  # Average number of plots per year, 1987-2008
  summary(surv.y$n.plots[surv.y$yr %in% 1987:2008])
  # Total number of plots surveyed between 2009-2004
  sum(surv.y$n.plots[surv.y$yr %in% 2009:2014])
  # Average number of plots per year, 2015-2020
  summary(surv.y$n.plots[surv.y$yr %in% 2015:2020])

#------------------------------------------------------------------------------# 
# Create capture histories for each tortoise
#------------------------------------------------------------------------------#

# MCL: Used to differentiate adults (>= 180 mm), juveniles (< 180 mm)

# Sex: 1 = female, 2 = male, 3 = juvenile/unknown 
  # When there were occasional discrepancies, we used sex at last capture, since 
  # we're more likely to correctly classify sex in larger individuals. Only
  # seven individuals that were captured as adults and not assigned sex
  
# Retain only a single capture each year (using MCL at first capture of year) 
# and identify stage (juvenile = 1, adult = 2)
  cr.yr <- cr %>%
    mutate(obsdate = ymd(obsdate)) %>%
    group_by(plot, tort, sex, yr) %>%
    summarize(mcl = MCL[1],
              .groups = "keep") %>%
    mutate(stage = ifelse(mcl < 180, 1, 2)) %>%
    data.frame()

# Identify how many unique tortoises were detected at least once during
# each year a plot was surveyed
  plot.yr <- cr.yr %>%
    group_by(plot, yr) %>%
    summarize(n.torts = length(tort),
              n.adultm = sum(sex == 2 & mcl > 179),
              n.adultf = sum(sex == 1 & mcl > 179),
              n.adultu = sum(sex == 3 & mcl > 179),
              n.juv = sum(mcl < 180),
              .groups = "keep") %>%
    data.frame()
  
# Identify how many unique tortoises were ever marked on each plot
  plot.tort <- cr.yr %>%
    group_by(plot) %>%
    summarize(n.torts = length(unique(tort))) %>%
    data.frame()

# Create matrix with capture histories 
  # 1 = captured as juvenile
  # 2 = captured as adult 
  # 3 = plot surveyed but tortoise not captured
  # NA = plot not surveyed
  
  ch <- cr.yr %>%
    select(-mcl) %>%
    add_row(plot = cr.yr$plot[1],
            tort = cr.yr$tort[1],
            sex = cr.yr$sex[1],
            yr = missingyrs,
            stage = NA) %>%
    pivot_wider(id_cols = c(plot, tort, sex),
                names_from = yr,
                names_prefix = "y",
                names_sort = TRUE,
                values_from = stage,
                values_fill = NA) %>%
    data.frame()

  # Change NAs to 3 when the plot was surveyed, but the tortoise wasn't captured
  surv.w.mat <- as.matrix(surv.w[, -1])
  ch.mat <- ch[ ,4:ncol(ch)]
  for (i in 1:nrow(ch)){
    years <- which(surv.w.mat[surv.w$code == ch$plot[i],] == 1)
    for(j in years){
      ch.mat[i,j] <- ifelse(is.na(ch.mat[i,j]), 3, ch.mat[i,j])
    }
  }
  ch <- cbind(ch[ ,1:3], ch.mat)
  
  # Check: number of tortoises captured each year at each plot 
  # (same totals from cr.yr and ch dataframes?)
    # chcheck <- ch %>%
    #   pivot_longer(cols = starts_with("y"),
    #                names_to = "yr",
    #                names_prefix = "y",
    #                names_transform = list(yr = as.integer),
    #                values_to = "stage") %>%
    #   data.frame()
    # plot.yr2 <- chcheck %>%
    #   filter(!is.na(stage) & stage < 3) %>%
    #   group_by(plot, yr) %>%
    #   summarize(n.torts = length(tort)) %>%
    #   data.frame()
    # all.equal(plot.yr$n.torts, plot.yr2$n.torts)
    # rm(chcheck, plot.yr2)

# Create vector indicating the first year each tortoise was caught:
  first <- rep(NA, nrow(ch.mat))
  for (i in 1:nrow(ch.mat)) {
    first[i] <- which(!is.na(ch.mat[i,]) & ch.mat[i,] < 3)[1]
  }
  
#------------------------------------------------------------------------------#  
# Format covariates
#------------------------------------------------------------------------------#

# Format individual covariates
  sex <- ch$sex
  sex[sex == 3] <- NA
  male.ind <- sex - 1

# Before formatting site covariate, check that plots are in the same order
# they're found in capture histories
  all.equal(plots$plot, unique(ch$plot))

# Format site covariate: distance to nearest major city (pop > 10,000)
  city.mn <- mean(plots$city.km)
  city.sd <- sd(plots$city.km)
  city <- (plots$city.km - city.mn) / city.sd
  
# Format site covariate: precipitation normals
  precipnorm.mn <- mean(plots$pptnorms.mm)
  precipnorm.sd <- sd(plots$pptnorms.mm)
  precip.norm <- (plots$pptnorms.mm - precipnorm.mn) / precipnorm.sd

# Format site*year covariate: drought
  # Calculate mean PDSI averaged over previous 24 months in each climate
  # division (index from 2013 paper)
  pdsi <- pdsi %>%
    arrange(div, yr, mon) %>%
    mutate(pdsi.24 = NA)
  jul.index <- which(pdsi$yr %in% 1988:2020 & pdsi$mon == 7)
  for (i in jul.index) {
    pdsi$pdsi.24[i] <- mean(pdsi$pdsi[(i-23):i])
  }
  pdsi.mn <- pdsi %>%
    filter(!is.na(pdsi.24)) %>%
    select(div, yr, pdsi.24) %>%
    rename(climate = div)
  
  # Link climate division data to plots
  pdsi.plot <- expand.grid(plot = plots$plot, 
                           yr = 1988:2020, 
                           stringsAsFactors = FALSE,
                           KEEP.OUT.ATTRS = FALSE) %>%
    left_join(plots[, c("plot", "climate")], by = "plot") %>%
    left_join(pdsi.mn, by = c("climate", "yr")) %>%
    mutate(pdsi.24.z = (pdsi.24 - mean(pdsi.24)) / sd(pdsi.24))
  # Save mean, SD for use later
  pdsi.24.mn <- mean(pdsi.plot$pdsi.24)
  pdsi.24.sd <- sd(pdsi.plot$pdsi.24)
  
  # Convert to wide format
  pdsi.24 <- pdsi.plot %>%
    pivot_wider(id_cols = plot,
                names_from = yr,
                names_prefix = "y",
                names_sort = TRUE,
                values_from = pdsi.24.z) %>%
    data.frame()
  
  # Check that plots are in the same order they're found in capture histories
  all.equal(pdsi.24$plot, unique(ch$plot))
  
  # Convert to matrix
  pdsi.24.mat <- as.matrix(pdsi.24[, -1])
  
  # Summarize how drought conditions at each plot in 1996-2019
  drought19962019 <- pdsi.plot %>%
    filter(yr %in% 1996:2019) %>%
    group_by(plot) %>%
    summarize(pdsi.24_mn = mean(pdsi.24),
              pdsi.24_min = min(pdsi.24),
              pdsi.24_max = max(pdsi.24),
              yrs_neg = sum(pdsi.24 < 0)) %>%
    data.frame()
  drought19962019

# Format site*year covariate: precipitation	
  # Calculate cumulative precipitation (mm) at each plot from Aug-Jul
  precip <- precip %>%
    mutate(yr = as.numeric(str_sub(yr.m, 1, 4)),
           mon = as.numeric(str_sub(yr.m, 6, 7)),
           season = ifelse(mon < 8, yr, yr + 1))
  precip.aj.plot <- precip %>%
    filter(season %in% 1988:2020) %>%
    group_by(plot, season) %>%
    summarize(ppt = sum(ppt),
              .groups = "keep") %>%
    data.frame() %>%
    mutate(ppt.z = (ppt - mean(ppt)) / sd(ppt))
  # Save mean, SD for use later
  ppt.mn <- mean(precip.aj.plot$ppt)
  ppt.sd <- sd(precip.aj.plot$ppt)
  
  # Convert to wide format
  precip.aj <- precip.aj.plot %>%
    pivot_wider(id_cols = plot,
                names_from = season,
                names_prefix = "y",
                names_sort = TRUE,
                values_from = ppt.z) %>%
    data.frame()
    
  # Check that plots are in the same order they're found in capture histories
  all.equal(precip.aj$plot, unique(ch$plot))
  
  # Convert to matrix
  precip.aj.mat <- as.matrix(precip.aj[, -1])  

# Format site*year covariate: survey effort 
# (persondays or persondays scaled by plot area)
  surveys <- surveys %>%
    select(-plot) %>%
    rename(plot = code) %>%
    mutate(area.sqkm = area.sqmi * 2.59,
           effort.sc = persondays / area.sqkm,
           persondays.z = (persondays - mean(persondays)) / sd(persondays),
           effort.sc.z = (effort.sc - mean(effort.sc)) / sd(effort.sc)) %>%
    arrange(plot)
 
  # Convert to wide format
  effort <- surveys %>%
    filter(yr %in% 1988:2020) %>%
    select(plot, yr, persondays.z) %>%
    add_row(plot = surveys$plot[1],
            yr = missingyrs,
            persondays.z = NA) %>%
    pivot_wider(id_cols = plot,
                names_from = yr,
                names_prefix = "y",
                names_sort = TRUE,
                values_from = persondays.z) %>%
    data.frame()
  effort.sc <- surveys %>%
    filter(yr %in% 1988:2020) %>%
    select(plot, yr, effort.sc.z) %>%
    add_row(plot = surveys$plot[1],
            yr = missingyrs,
            effort.sc.z = NA) %>%
    pivot_wider(id_cols = plot,
                names_from = yr,
                names_prefix = "y",
                names_sort = TRUE,
                values_from = effort.sc.z) %>%
    data.frame()
  
  # Check that plots are in the same order they're found in capture histories
  all.equal(effort$plot, unique(ch$plot))
  
  # Convert to matrix and replace NAs with -9999 (need some value since JAGS
  # does not take NAs for covariates, but values for years when surveys weren't 
  # conducted won't be used)
  effort.mat <- as.matrix(effort[, -1])  
  effort.mat[is.na(effort.mat)] <- -9999
  effort.sc.mat <- as.matrix(effort.sc[, -1])
  effort.sc.mat[is.na(effort.sc.mat)] <- -9999
  
# Create plot index for each tortoise
  plots$plot.index <- 1:nrow(plots)
  plot.index <- plots$plot.index[match(ch$plot, plots$plot)]	
  
#------------------------------------------------------------------------------#  
# Functions to create initial values
#------------------------------------------------------------------------------#

# Create vector indicating the first year each tortoise was caught as juvenile:
  first1 <- rep(NA, nrow(ch.mat))
  for (i in 1:nrow(ch.mat)) {
    first1[i] <- which(!is.na(ch.mat[i,]) & ch.mat[i,] == 1)[1]
  }
  
# Create vector indicating the first year each tortoise was caught as adult:  
  first2 <- rep(NA, nrow(ch.mat))
  for (i in 1:nrow(ch.mat)) {
    first2[i] <- which(!is.na(ch.mat[i,]) & ch.mat[i,] == 2)[1]
  }
  
# Create a matrix of initial values for latent states (z)
  # NAs up to and including the first occasion, then 1/2 through the remainder
  ch.init <- function(y, f1, f2){
    for (i in 1:length(f1)) {
      if(!is.na(f1[i]) & f1[i] == ncol(y)) {y[i,] <- NA} else
        if (is.na(f1[i]) & !is.na(f2[i]) & f2[i] == ncol(y)) {y[i,] <- NA} else
          if (is.na(f1[i]) & !is.na(f2[i]) & f2[i] != ncol(y)) {y[i, 1:f2[i]] <- NA; y[i, (f2[i] + 1):ncol(y)] <- 2} else
            if (!is.na(f1[i]) & f1[i] != ncol(y) & is.na(f2[i])) {y[i, 1:f1[i]] <- NA; y[i, (f1[i] +1 ):ncol(y)] <- 1} else
              if (!is.na(f1[i]) & !is.na(f2[i]) & (f2[i] - f1[i] == 1)) {y[i, 1:f1[i]] <- NA; y[i, f2[i]:ncol(y)] <- 2} else
                {y[i, 1:f1[i]] <- NA; y[i,(f1[i] + 1):(f2[i] - 1)] <- 1; y[i, f2[i]:ncol(y)] <- 2}}
    return(y)
  }
  
#------------------------------------------------------------------------------#  
# Run multistate model
#------------------------------------------------------------------------------#  

# Prep data objects
  ntorts <- nrow(ch.mat)             # Number of tortoises
  nyears <- ncol(ch.mat)             # Number of occasions
  nplots <- length(unique(ch$plot))  # Number of plots
  
  # Sequence to evaluate a (logit) linear trend in adult survival
  trend <- seq(0, nyears - 2)        
  trend.z <- (trend - mean(trend)) / sd(trend)
  trend.z2 <- trend.z * trend.z
  
  tortdata <- list(y = as.matrix(ch.mat),
                   ntorts = ntorts,
                   nyears = nyears,
                   nplots = nplots,
                   first = first,
                   male = male.ind,
                   plot = plot.index,
                   city = city,
                   mn.precip = precip.norm,
                   drought = pdsi.24.mat,
                   trend = trend.z,
                   trend2= trend.z2,
                   precip = precip.aj.mat,
                   effort = effort.mat)
  
  # Model
  tortcode <- nimbleCode({
    
    #-- Priors and constraints
    
    alpha.p1 ~ dlogis(0, 1)
    alpha.p2 ~ dlogis(0, 1)
    beta.phi1 ~ dlogis(0, 1)
    beta.phi2 ~ dlogis(0, 1)
    gamma.psi ~ dlogis(0, 1)
    
    a1.precip ~ dnorm(0, 0.1)
    a1.effort ~ dnorm(0, 0.1)
    a2.male ~ dnorm(0, 0.1)
    a2.precip ~ dnorm(0, 0.1)
    a2.effort ~ dnorm(0, 0.1)
    
    b1.city ~ dnorm(0, 0.1)
    b1.mnprecip ~ dnorm(0, 0.1)
    b1.drought ~ dnorm(0, 0.1)
    b1.int ~ dnorm(0, 0.1)
    b2.male ~ dnorm(0, 0.1)
    b2.city ~ dnorm(0, 0.1)
    b2.mnprecip ~ dnorm(0, 0.1)
    b2.drought ~ dnorm(0, 0.1)
    b2.int ~ dnorm(0, 0.1)
    b2.trend ~ dnorm(0, 0.1)
    b2.trend2 ~ dnorm(0, 0.1)
    
    c.mnprecip ~ dnorm(0, 0.1)
    
    omega ~ dunif(0, 1)
    
    sigma.site.1 ~ T(dt(0, pow(2.5, -2), 1), 0, )
    tau.site.1 <- 1 / (sigma.site.1 * sigma.site.1)
    
    sigma.site.2 ~ T(dt(0, pow(2.5, -2), 1), 0, )
    tau.site.2 <- 1 / (sigma.site.2 * sigma.site.2)
    
    sigma.site.t ~ T(dt(0, pow(2.5, -2), 1), 0, )
    tau.site.t <- 1 / (sigma.site.t * sigma.site.t)
    
    for (p in 1:nplots) {
      e.site.1[p] ~ dnorm(0, tau.site.1)
      e.site.2[p] ~ dnorm(0, tau.site.2)
      e.site.t[p] ~ dnorm(0, tau.site.t)
    }
    
    for (i in 1:ntorts) {
      for (t in first[i]:(nyears - 1)) {
        
        # Juvenile recapture probability
        logit(p1[i, t]) <- alpha.p1 + a1.precip * precip[plot[i], t] + 
          a1.effort * effort[plot[i], t]
        # Juvenile survival probability
        logit(phi1[i, t]) <- beta.phi1 + b1.city * city[plot[i]] +
          b1.mnprecip * mn.precip[plot[i]] + 
          b1.drought * drought[plot[i], t] +
          b1.int * mn.precip[plot[i]] * drought[plot[i], t] +
          e.site.1[plot[i]]
        
        # Transition probability
        logit(psi12[i, t]) <- gamma.psi + c.mnprecip * mn.precip[plot[i]] +
          e.site.t[plot[i]]
        
        # Adult recapture probability
        logit(p2[i, t]) <- alpha.p2 + a2.male * male[i] + 
          a2.precip * precip[plot[i],t] + 
          a2.effort * effort[plot[i],t]
        # Adult survival probability
        logit(phi2[i, t]) <- beta.phi2 + b2.male * male[i] + 
          b2.city * city[plot[i]] + 
          b2.mnprecip * mn.precip[plot[i]] +
          b2.drought * drought[plot[i], t] + 
          b2.int * mn.precip[plot[i]] * drought[plot[i], t] +
          b2.trend * trend[t] + b2.trend2 * trend2[t] + 
          e.site.2[plot[i]]
        
        # Define state transition probabilities
        # First index = state at time t-1, last index = state at time t
        ps[1, i ,t, 1] <- phi1[i, t] * (1 - psi12[i, t])
        ps[1, i, t, 2] <- phi1[i, t] * psi12[i, t]
        ps[1, i, t, 3] <- 1-phi1[i, t]
        ps[2, i, t, 1] <- 0
        ps[2, i, t, 2] <- phi2[i, t]
        ps[2, i, t, 3] <- 1-phi2[i, t]
        ps[3, i, t, 1] <- 0
        ps[3, i, t, 2] <- 0
        ps[3, i, t, 3] <- 1
        
        # Define stage-dependent detection probabilities
        # First index = state at time t, last index = detection type at time t
        po[1, i, t, 1] <- p1[i, t]
        po[1, i, t, 2] <- 0
        po[1, i, t, 3] <- 1-p1[i, t]
        po[2, i, t, 1] <- 0
        po[2, i, t, 2] <- p2[i, t]
        po[2, i, t, 3] <- 1-p2[i, t]
        po[3, i, t, 1] <- 0
        po[3, i, t, 2] <- 0
        po[3, i, t, 3] <- 1
        
      } # t
    } # i
    
    #-- Likelihood
    
    for (i in 1:ntorts) {
      z[i, first[i]] <- y[i, first[i]]
      male[i] ~ dbern(omega)
      
      for (t in (first[i] + 1):nyears) {
        
        # State process: draw State[t] given State[t-1]
        z[i, t] ~ dcat(ps[z[i, t-1], i, t-1, 1:3])
        
        # Observation process: draw Obs[t] given State[t]
        y[i, t] ~ dcat(po[z[i, t], i, t-1, 1:3])
        
      } # t
    } # i
  })
  
  # Parameters to monitor
  params <- c("alpha.p1", "a1.precip", "a1.effort",
              "beta.phi1", "b1.city", "b1.mnprecip", "b1.drought", "b1.int",
              "gamma.psi", "c.mnprecip",
              "alpha.p2", "a2.male", "a2.precip", "a2.effort",
              "beta.phi2", "b2.male", "b2.city", "b2.mnprecip", 
              "b2.drought", "b2.int", "b2.trend", "b2.trend2",
              "omega", "sigma.site.1", "e.site.1", "sigma.site.2", "e.site.2",
              "sigma.site.t", "e.site.t")
  
  # Initial values
  inits <- function() {list(alpha.p1 = runif(1, -1, 1),
                            alpha.p2 = runif(1, -1, 2),
                            beta.phi1 = runif(1, -1, 1),
                            beta.phi2 = runif(1, 1, 3),
                            gamma.psi = runif(1, -2, 0),
                            a1.precip = runif(1, -0.5, 0.5),
                            a1.effort = runif(1, -0.5, 0.5),
                            a2.male = runif(1, -0.5, 0.5),
                            a2.precip = runif(1, -0.5, 0.5),
                            a2.effort = runif(1, -0.5, 0.5),
                            b1.city = runif(1, -0.5, 0.5),
                            b1.mnprecip = runif(1, -0.5, 0.5),
                            b1.drought = runif(1, -0.5, 0.5),
                            b1.int = runif(1, -0.5, 0.5),
                            b2.male = runif(1, -0.5, 0.5),
                            b2.city = runif(1, -0.5, 0.5),
                            b2.mnprecip = runif(1, -0.5, 0.5),
                            b2.drought = runif(1, -0.5, 0.5),
                            b2.int = runif(1, -0.5, 0.5),
                            b2.trend = runif(1, -0.5, 0.5),
                            b2.trend2 = runif(1, -0.5, 0.5),
                            c.mnprecip = runif(1, -0.5, 0.5),
                            omega = runif(1, 0, 1),
                            sigma.site.1 = runif(1, 0, 3),
                            sigma.site.2 = runif(1, 0, 3),
                            sigma.site.t = runif(1, 0, 3),
                            male = ifelse(is.na(male.ind), 1, NA),
                            z = ch.init(as.matrix(ch.mat), first1, first2))}
  
# Separate constants and data for NIMBLE
  tortconstants <- tortdata
  tortconstants[c("male", "y")] <- NULL

# Create model object (~ 15 min)
  # tortmodel <- nimbleModel(code = tortcode, constants = tortconstants, 
  #                          calculate = FALSE)

# Set data and inits
  # tortmodel$setData(list(y = as.matrix(ch.mat),
  #                        male = male.ind))
  # set.seed(123)
  # tortmodel$setInits(inits())

# Build MCMC (~ 69 min)
  # tortmcmc <- buildMCMC(tortmodel,
  #                       monitors = params)

# Compile the model and MCMC (~ 21 min)
  # Ctortmodel <- compileNimble(tortmodel)
  # Ctortmcmc <- compileNimble(tortmcmc, project = tortmodel)

# MCMC settings
  # n.chains <- 3
  # n.iter <- 30000
  # n.burn <- 5000
  # n.thin <- 15
  # ni.tot <- n.iter + n.burn
    
# Run the MCMC and extract the samples (~ 13 hrs)
  # samples <- runMCMC(
  #   Ctortmcmc,
  #   nchains = n.chains,
  #   niter = ni.tot,
  #   nburnin = n.burn,
  #   thin = n.thin,
  #   samplesAsCodaMCMC = TRUE
  # )

# Save samples
  # saveRDS(samples, "MS-samples-6000.rds")

#------------------------------------------------------------------------------#  
# Post-processing
#------------------------------------------------------------------------------#    
  
# Load samples from previous run (if needed)
  samples <- readRDS("MS-samples-6000.rds")

# Produce summary table, look at trace & density plots
  MCMCsummary(samples, 
              round = 2, 
              params = "all", 
              probs = c(0.025, 0.975))
  # MCMCtrace(samples,
  #           params = "all",
  #           pdf = TRUE,
  #           open_pdf = FALSE)
  # MCMCplot(samples,
  #          params = "all", # excl = ""
  #          ci = c(50, 90))  

# Create matrix, dataframe with samples  
  samples_mat <- MCMCchains(samples, 
                            params = "all",
                            mcmc.list = FALSE)
  samples_df <- data.frame(samples_mat)
  
# Create a table with parameter estimates for the manuscript
  params <- c("beta.phi1", "b1.city", "b1.mnprecip", "b1.drought", "b1.int", 
              "sigma.site.1", "beta.phi2", "b2.male", "b2.city", "b2.mnprecip",
              "b2.drought", "b2.int", "b2.trend", "b2.trend2", "sigma.site.2",
              "gamma.psi", "c.mnprecip", "sigma.site.t", "alpha.p1", "a1.precip",
              "a1.effort", "alpha.p2", "a2.male", "a2.precip", "a2.effort")
  param_table <- data.frame(parameter = params) %>%
    mutate(Mean = apply(samples_mat[, params], 2, mean),
           lcl = apply(samples_mat[, params], 2, quantile, probs = 0.025),
           ucl = apply(samples_mat[, params], 2, quantile, probs = 0.975),
           greater0 = apply(samples_mat[, params], 2, function(x) sum(x > 0)/length(x)),
           f = if_else(greater0 < 0.5, 1 - greater0, greater0)) %>%
    select(-greater0)

# Write to file  
# write.csv(param_table, "output/parameter-estimates.csv", row.names = FALSE)

#------------------------------------------------------------------------------# 
# Table: plot characteristics, capture data
#------------------------------------------------------------------------------#  

# Identify first and last known state of each tortoise 
# (state1 = first; state2 = last)
  firstlast <- ch[, 1:3]
  firstlast$state1 <- unlist(apply(ch.mat, 1, 
                                   function(x) head(x[!is.na(x) & x < 3], 1)))
  firstlast$state2 <- unlist(apply(ch.mat, 1, 
                                   function(x) tail(x[!is.na(x) & x < 3], 1)))
  # Check that there are no backward transitions (ad -> juv)
  count(firstlast, state1, state2)
  
# Create new sex column with all tortoises that were only captured as juveniles
# (ie, never captured as an adult) listed as "J"
  firstlast <- firstlast %>%
    mutate(sex_new = if_else(state2 == 1, "J",
                             if_else(sex == 1, "F",
                                     if_else(sex == 2, "M", "U")))) %>%
    select(-sex) %>%
    rename(sex = sex_new)
  
# Calculate for each plot:
  # n_torts = total number of marked tortoises
  # n_juv = total number of marked tortoises first captured as juveniles
  # prop_juv = proportion of marked tortoises first captured as juveniles
  # n_transition = number of marked tortoises first captured as juveniles that
    # were subsequently captured as an adult
  # prop_transition = proportion of juveniles captured as adults
  # n_adults = number of tortoises ever captured as an adult
  # prop_female = proportion of adults that were female
  
caps_byplot <- firstlast %>%
  group_by(plot) %>%
  summarize(n_torts = length(tort),
            n_juv = sum(state1 == 1),
            prop_juv = n_juv / n_torts,
            n_transition = sum(state1 == 1 & state2 == 2),
            prop_transition = n_transition / n_juv,
            n_adults = sum(sex != "J"),
            prop_female = sum(sex == "F") / n_adults) %>%
  data.frame()

plot_summaries <- plots %>%
  select(-plot.index) %>%
  left_join(caps_byplot, by = "plot") %>%
  arrange(plot)
to_add <- data.frame(plot = "Mean/Total", climate = NA, 
                     city.km = mean(plot_summaries$city.km),
                     pptnorms.mm = mean(plot_summaries$pptnorms.mm),
                     n_torts = sum(plot_summaries$n_torts),
                     n_juv = sum(plot_summaries$n_juv),
                     prop_juv = sum(firstlast$state1 == 1) / nrow(firstlast),
                     n_transition = sum(plot_summaries$n_transition),
                     prop_transition = sum(firstlast$state1 == 1 & firstlast$state2 == 2) / sum(firstlast$state1 == 1),
                     n_adults = sum(plot_summaries$n_adults),
                     prop_female = sum(firstlast$sex == "F") / sum(firstlast$sex != "J"))
plot_summaries <- rbind(plot_summaries, to_add)

# Write to file
# write.csv(plot_summaries, "output/plot-summaries.csv", row.names = FALSE)

#------------------------------------------------------------------------------# 
# Settings for figures
#------------------------------------------------------------------------------#	

# Use mean or median for measure of central tendency
  ctend <- mean
  # ctend <- median

# 90% or 95% credible intervals
  # qprobs <- c(0.05,0.95)
  qprobs <- c(0.025,0.975) 

#------------------------------------------------------------------------------# 
# Effect of drought on survival
#------------------------------------------------------------------------------#	
# Note: not taking site-level random effects into account in CI calculations
  
# Adults
  # Use survival estimates for last year (2019-2020)
  # Assume average distance from city (standardized distance = 0)
  phi2 <- samples_df %>% 
    select(beta.phi2, b2.male, b2.mnprecip, b2.drought, b2.int, b2.trend, 
           b2.trend2) %>%
    as.matrix()

  # Generate covariate values for figure
  # Use 3 values of mnprecip (values range from 177-409; mean = 280)
  mnprecip3 <- c(180, 280, 380)
  mnprecip3.z <- (mnprecip3 - precipnorm.mn) / precipnorm.sd
  xmin <- min(pdsi.24.mat)
  xmax <- max(pdsi.24.mat)
  predx <- cbind(int = 1, male = rep(0:1, each = 300),
                 mnprecip = rep(rep(mnprecip3.z, 2), each = 100),
                 drought = rep(seq(xmin, xmax, length = 100), 6))
  predx <- cbind(predx, interact = predx[,3] * predx[,4],
                 trend = tail(trend.z, 1), trend2 = tail(trend.z2, 1))
  predl <- predx %*% t(phi2)  # [600 * 7] %*% [7 * 6000] = [600 * 6000]
  pred <- exp(predl) / (1 + exp(predl))  
  
  # Calculate mean/median and CRIs:
  adsurv_df <- as.data.frame(predx) %>%
    mutate(central = apply(pred, 1, ctend),
           lcl = apply(pred, 1, quantile, probs = qprobs[1]),
           ucl = apply(pred, 1, quantile, probs = qprobs[2]),
           sex = ifelse(male == 1, "M", "F"),
           regime = ifelse(mnprecip == mnprecip3.z[1], "Arid",
                              ifelse(mnprecip == mnprecip3.z[2], 
                                     "Arid/Semiarid", "Semiarid")),
           group = paste(regime, sex, sep = ":"),
           pdsi = drought * pdsi.24.sd + pdsi.24.mn)

  # Set colors, linetypes for figure (types: 1 = solid, 2 = dashed, 3 = dotted) 
  groups <- data.frame(unique(adsurv_df[, c("group", "sex", "regime")])) %>%
    mutate(col = ifelse(sex == "F", "salmon3", "steelblue4"),
           linetype = ifelse(regime == "Arid", 1,
                             ifelse(regime == "Semiarid", 2, 3)))
  linewidth <- 0.3
  
  # Just using estimates for M/F at dries and wettest sites
  adsurv_df4 <- filter(adsurv_df, regime != "Arid/Semiarid") %>%
    mutate(group = as.factor(group))
  groups4 <- filter(groups, regime != "Arid/Semiarid") %>%
    mutate(group = as.factor(group)) %>%
    arrange(group)

  adsurv_plot <- ggplot(data = adsurv_df4) +
    geom_vline(xintercept = 0, col = "gray", linetype = 3, linewidth = linewidth) +
    geom_line(aes(x = pdsi, y = central, group = group,
                  color = group, linetype = group), linewidth = linewidth) +     
    labs(x = "PDSI (24-month)", y = "Estimated adult survival") +
    scale_x_continuous(breaks = seq(-4, 5, by = 2)) +
    scale_color_manual(values = groups4$col) +
    scale_linetype_manual(values = groups4$linetype) +
    theme_classic() +
    theme(text = element_text(size = 9),
          axis.text = element_text(size = 9),
          legend.title = element_blank(),
          legend.position = c(1, 0.02),
          legend.justification = c("right", "bottom"),
          legend.key.height = unit(0.15, "in"))

  # ggsave("output/adult-survival-drought.jpg",
  #        adsurv_plot,
  #        device = "jpeg",
  #        dpi = 600,
  #        width = 3,
  #        height = 3,
  #        units = "in")
  # Can easily change device to pdf (and remove dpi argument)

# Juveniles
  # Use survival estimates for last year (2019-2020)
  # Assume average distance from city (standardized distance = 0)
  phi1 <- samples_df %>% 
    select(beta.phi1, b1.mnprecip, b1.drought, b1.int) %>%
    as.matrix()
  
  # Generate covariate values for figure
  pred1x <- cbind(int = 1, mnprecip = rep(mnprecip3.z, each=100),
                 drought = rep(seq(xmin, xmax, length = 100), 3))
  pred1x <- cbind(pred1x, interact = pred1x[,2] * pred1x[,3])
  pred1l <- pred1x %*% t(phi1)
  pred1 <- exp(pred1l) / (1 + exp(pred1l))
  
  # Calculate mean/median and CRIs:
  juvsurv_df <- as.data.frame(pred1x) %>%
    mutate(central = apply(pred1, 1, ctend),
           lcl = apply(pred1, 1, quantile, probs = qprobs[1]),
           ucl = apply(pred1, 1, quantile, probs = qprobs[2]),
           regime = ifelse(mnprecip == mnprecip3.z[1], "Arid",
                           ifelse(mnprecip == mnprecip3.z[2], 
                                  "Arid/Semiarid", "Semiarid")),
           pdsi = drought * pdsi.24.sd + pdsi.24.mn)  
  
  # Set colors, linetypes for figure (types: 1 = solid, 2 = dashed, 3 = dotted) 
  groups_juv <- data.frame(regime = unique(juvsurv_df$regime)) %>%
    mutate(col = "black",
           linetype = ifelse(regime == "Arid", 1,
                             ifelse(regime == "Semiarid", 2, 3)))
  linewidth <- 0.3
  
  # Just using estimates  at driest and wettest sites
  juvsurv_df2 <- filter(juvsurv_df, regime != "Arid/Semiarid") %>%
    mutate(regime = as.factor(regime))
  groups_juv2 <- filter(groups_juv, regime != "Arid/Semiarid") %>%
    mutate(regime = as.factor(regime)) %>%
    arrange(regime)
  
  juvsurv_plot <- ggplot(data = juvsurv_df2) +
    geom_vline(xintercept = 0, col = "gray", linetype = 3, linewidth = linewidth) +
    geom_line(aes(x = pdsi, y = central, group = regime,
                  color = regime, linetype = regime), linewidth = linewidth) +     
    labs(x = "PDSI (24-month)", y = "Estimated juvenile survival") +
    scale_x_continuous(breaks = seq(-4, 5, by = 2)) +
    scale_y_continuous(limits = c(0.47, 1), breaks = seq(0.5, 1, by = 0.1), 
                       labels = c("0.50", "0.60", "0.70", "0.80", "0.90", "1.00")) +
    scale_color_manual(values = groups_juv2$col) +
    scale_linetype_manual(values = groups_juv2$linetype) +
    theme_classic() +
    theme(text = element_text(size = 9),
          axis.text = element_text(size = 9),
          legend.title = element_blank(),
          legend.position = c(1, 0.02),
          legend.justification = c("right", "bottom"),
          legend.key.height = unit(0.15, "in"))

  # ggsave("output/juv-survival-drought.jpg",
  #        juvsurv_plot,
  #        device = "jpeg",
  #        dpi = 600,
  #        width = 3,
  #        height = 3,
  #        units = "in")
 
# Adults and juveniles, stacked
  # Remove x-axis text from juveniles (top) figure
  juvsurv_stack <- ggplot(data = juvsurv_df2) +
    geom_vline(xintercept = 0, col = "gray", linetype = 3, linewidth = linewidth) +
    geom_line(aes(x = pdsi, y = central, group = regime,
                  color = regime, linetype = regime), linewidth = linewidth) +     
    labs(x = "PDSI (24-month)", y = "Estimated juvenile survival") +
    scale_x_continuous(breaks = seq(-4, 5, by = 2)) +
    scale_y_continuous(limits = c(0.47, 1), breaks = seq(0.5, 1, by = 0.1), 
                       labels = c("0.50", "0.60", "0.70", "0.80", "0.90", "1.00")) +
    scale_color_manual(values = groups_juv2$col) +
    scale_linetype_manual(values = groups_juv2$linetype) +
    theme_classic() +
    theme(text = element_text(size = 9),
          axis.text = element_text(size = 9),
          legend.title = element_blank(),
          legend.position = c(1, 0.02),
          legend.justification = c("right", "bottom"),
          legend.key.height = unit(0.15, "in"),
          axis.title.x = element_blank(), 
          axis.text.x = element_blank())
  surv_stack <- plot_grid(juvsurv_stack, adsurv_plot, ncol = 1)

  # ggsave("output/both-survival-drought.jpg",
  #        surv_stack,
  #        device = "jpeg",
  #        dpi = 600,
  #        width = 3,
  #        height = 5.5,
  #        units = "in")

# May want to simplify legends...
  
#------------------------------------------------------------------------------# 
# Temporal trends in adult survival
#------------------------------------------------------------------------------#	  
# Again, not taking site-level random effects into account for CI calculations

# For an average site (mean distance from city, annual precip) with PDSI = mean
# value during study period (-1.08; drought.z = 0)
phi2t <- samples_df %>% 
  select(beta.phi2, b2.male, b2.trend, b2.trend2) %>%
  as.matrix()
  
predtx <- cbind(int = 1, male = rep(0:1, each = length(trend)), 
                trend = rep(trend.z, 2), trend2 = rep(trend.z2, 2))
predtl <- predtx %*% t(phi2t)
predt <- exp(predtl) / (1 + exp(predtl))  

# Calculate mean/median and CRIs:
adtrend_df <- as.data.frame(predtx) %>%
  mutate(endyr = rep(1988:2020, 2),
         interval = rep(paste(1987:2019, 1988:2020, sep = "-"), 2),
         central = apply(predt, 1, ctend),
         lcl = apply(predt, 1, quantile, probs = qprobs[1]),
         ucl = apply(predt, 1, quantile, probs = qprobs[2]),
         sex = ifelse(male == 1, "M", "F"))

# Set colors, linetypes for figure (types: 1 = solid, 2 = dashed, 3 = dotted) 
groups_tr <- data.frame(sex = unique(adtrend_df$sex)) %>%
  mutate(col = ifelse(sex == "F", "salmon3", "steelblue4"),
         fill = col, 
         linetype = 1)
linewidth <- 0.3

trend_plot <- ggplot(data = adtrend_df, 
                     aes(x = endyr, group = sex)) +
  geom_line(aes(y = central, color = sex, linetype = sex), linewidth = linewidth) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = sex), alpha = 0.2) +
  labs(x = "Year", y = "Estimated adult survival (95% CI)") +
  scale_color_manual(values = groups_tr$col) +
  scale_linetype_manual(values = groups_tr$linetype) +
  scale_fill_manual(values = groups_tr$fill) +
  scale_y_continuous(limits = c(0.75, 1), breaks = seq(0.75, 1, by = 0.05)) +
  theme_classic() +
  theme(text = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.title = element_blank(),
        legend.position = c(1, 0.02),
        legend.justification = c("right", "bottom"),
        legend.key.height = unit(0.15, "in"))

# ggsave("output/adult-survival-trends.jpg",
#        trend_plot,
#        device = "jpeg",
#        dpi = 600,
#        width = 3,
#        height = 3,
#        units = "in")

#------------------------------------------------------------------------------# 
# Temporal trend in PDSI values
#------------------------------------------------------------------------------#	

pdsi24t <- pdsi.plot %>%
  rename(div = climate) %>%
  select(yr, div, pdsi.24) %>%
  distinct() %>%
  mutate(yr0 = yr - min(yr))

# Quick run of ML linear regression models
summary(lm.trend1 <- lm(pdsi.24 ~ yr0, data = pdsi24t))  
summary(lm.trend5 <- lm(pdsi.24 ~ yr0 * factor(div), data = pdsi24t))
AIC(lm.trend1); AIC(lm.trend5) # AIC 10 points lower for the simpler model

# Bayesian linear regression
regcode <- nimbleCode({
  b0 ~ dnorm(0, sd = 100)
  b1 ~ dnorm(0, sd = 100)
  sigma ~ dunif(0, 100)

  for(i in 1:nobs){
    y[i] ~ dnorm(b0 + b1*x[i], sd = sigma)
  }
})

regdata <- list(y = pdsi24t$pdsi.24)
regconstants <- list(x = pdsi24t$yr0,
                     nobs = nrow(pdsi24t))

params <- c("b0", "b1", "sigma")

inits <- function() {list(b0 = runif(1, -10, 10),
                          b1 = runif(1, -2, 2),
                          sigma = runif(1, 1, 10))} 

# Create model obejct
set.seed(123)
regmodel <- nimbleModel(code = regcode, 
                        data = regdata,
                        constants = regconstants,
                        inits = inits())

# Build MCMC
regmcmc <- buildMCMC(regmodel, monitors = params)

# Compile the model and MCMC
Cregmodel <- compileNimble(regmodel)
Cregmcmc <- compileNimble(regmcmc, project = regmodel)

# MCMC settings
n.chains <- 3
n.iter <- 10000
n.burn <- 30000
n.thin <- 1
ni.tot <- n.iter + n.burn

# Run the MCMC and extract the samples
regsamples <- runMCMC(
  Cregmcmc,
  nchains = n.chains,
  niter = ni.tot,
  nburnin = n.burn,
  thin = n.thin,
  samplesAsCodaMCMC = TRUE
)
MCMCsummary(regsamples, 
            round = 2, 
            params = "all", 
            probs = c(0.025, 0.975))
MCMCtrace(regsamples,
          params = "all",
          pdf = FALSE)

# Create matrix with samples, and thin so left with 6000 samples
regsamples_mat <- MCMCchains(regsamples, 
                             params = "all",
                             mcmc.list = FALSE)
regsamples_mat <- regsamples_mat[seq(1, nrow(regsamples_mat), by = 5),]
betas <- regsamples_mat[, c("b0", "b1")]

pdsi_df <- data.frame(int = 1,yr0 = seq(0, 32, length=100))
pdsipreds <- as.matrix(pdsi_df) %*% t(betas)
pdsi_df <- pdsi_df %>%
  mutate(central = apply(pdsipreds, 1, ctend),
         lcl = apply(pdsipreds, 1, quantile, probs = qprobs[1]),
         ucl = apply(pdsipreds, 1, quantile, probs = qprobs[2]),
         yr = yr0 + 1988)

# Set colors, linetypes for figure
pdsi24t$div <- as.factor(pdsi24t$div)
divs <- data.frame(div = sort(unique(pdsi24t$div))) %>%
  mutate(col = c('mediumpurple4','steelblue4','darkseagreen4',
                 'goldenrod4','salmon4'),
         shape = 21:25)
linetype <- 1
linewidth <- 0.3

pdsi_plot <- ggplot() +
  geom_line(data = pdsi24t, aes(x = yr, y = pdsi.24, group = div, col = div),
            linewidth = linewidth, alpha = 0.6, show.legend = NA) +
  geom_point(data = pdsi24t, aes(x = yr, y = pdsi.24, col = div, fill = div, shape = div), 
             size = 1.2, alpha = 0.6) +
  geom_line(data = pdsi_df, aes(x = yr, y = central), show.legend = NA) +
  geom_ribbon(data = pdsi_df, aes(x = yr, ymin = lcl, ymax = ucl), alpha = 0.2,
              show.legend = NA) +
  labs(x = "Year", y = "PDSI (24-month)") +
  scale_x_continuous(limits = c(1988, 2020), expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(-4, 4.7), breaks = seq(-4, 4, by = 2)) +
  scale_color_manual(values = divs$col, name = "Division") +
  scale_fill_manual(values = divs$col, name = "Division") +
  scale_shape_manual(values = divs$shape, name = "Division") +
  theme_classic() +
  theme(text = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.position = c(0.92, 1),
        legend.justification = c("right", "top"),
        legend.key.height = unit(0.12, "in"))

# ggsave("output/pdsi-trend.jpg",
#        pdsi_plot,
#        device = "jpeg",
#        dpi = 600,
#        width = 6.5,
#        height = 3,
#        units = "in")

#------------------------------------------------------------------------------# 
# Plot-specific estimates of demographic rates, lambdas
#------------------------------------------------------------------------------#
# For adult survival, generating values for 2019-2020
# For survival, generating values for 4 levels of drought: 
  # mean PDSI = -3, -1.08 (mean value across plots and years), 0, +3
# Will account for site-level random effects in all demographic parameters
# For lambda, will assume recruitment = 0.32 f/f/yr 

drought4 <- c(-3, pdsi.24.mn, 0, 3)
drought4.z <- (drought4 - pdsi.24.mn) / pdsi.24.sd

# Juvenile survival
  # Grab random effects for each plot
  phi1_RE <- select(samples_df, contains("e.site.1"))
  # Grab posterior samples for fixed parameters
  phi1_samples <- samples_df %>% 
    select(beta.phi1, b1.city, b1.mnprecip, b1.drought, b1.int)
  # Gather covariate values
  phi1_X <- data.frame(plot = rep(plots$plot, 4),
                       intercept = 1,
                       city = rep(city, 4),
                       mnprecip = rep(precip.norm, 4),
                       drought = rep(drought4.z, each = nrow(plots))) %>%
    mutate(int = mnprecip * drought)
    # check: 
    # head(phi1_X, 18)
  # Matrix math
  phi1_l <- as.matrix(select(phi1_X, -plot)) %*% t(phi1_samples)
  # Add in random effects
  phi1_l <- phi1_l + rbind(t(phi1_RE), t(phi1_RE), t(phi1_RE), t(phi1_RE))
  # Put on probability scale 
  phi1_mat <- exp(phi1_l) / (1 + exp(phi1_l))
  
# Adult male survival
  # Grab random effects for each plot
  phi2_RE <- select(samples_df, contains("e.site.2"))
  # Grab posterior samples for fixed parameters
  phi2_samples <- samples_df %>% 
    select(beta.phi2, b2.male, b2.city, b2.mnprecip, b2.drought, b2.int, 
           b2.trend, b2.trend2)
  # Gather covariate values
  phi2m_X <- data.frame(plot = rep(plots$plot, 4),
                        intercept = 1,
                        male = 1,
                        city = rep(city, 4),
                        mnprecip = rep(precip.norm, 4),
                        drought = rep(drought4.z, each = nrow(plots))) %>%
    mutate(int = mnprecip * drought) %>%
    mutate(trend = tail(trend.z, 1),
           trend2 = trend * trend)
    # check: 
    # head(phi2m_X, 18)
  # Matrix math
  phi2m_l <- as.matrix(select(phi2m_X, -plot)) %*% t(phi2_samples)
  # Add in random effects
  phi2m_l <- phi2m_l + rbind(t(phi2_RE), t(phi2_RE), t(phi2_RE), t(phi2_RE))
  # Put on probability scale 
  phi2m_mat <- exp(phi2m_l) / (1 + exp(phi2m_l))  
  
# Adult female survival
  # Gather covariate values
  phi2f_X <- phi2m_X %>% mutate(male = 0)
  # Matrix math
  phi2f_l <- as.matrix(select(phi2f_X, -plot)) %*% t(phi2_samples)
  # Add in random effects
  phi2f_l <- phi2f_l + rbind(t(phi2_RE), t(phi2_RE), t(phi2_RE), t(phi2_RE))
  # Put on probability scale 
  phi2f_mat <- exp(phi2f_l) / (1 + exp(phi2f_l))  
  
# Transition rate
  # Grab random effects for each plot
  psi_RE <- select(samples_df, contains("e.site.t"))
  # Grab posterior samples for fixed parameters
  psi_samples <- samples_df %>% 
    select(gamma.psi, c.mnprecip)
  # Gather covariate values
  psi_X <- data.frame(plot = plots$plot,
                      intercept = 1,
                      mnprecip = precip.norm)
    # check: 
    # psi_X
  # Matrix math
  psi_l <- as.matrix(select(psi_X, -plot)) %*% t(psi_samples)
  # Add in random effects
  psi_l <- psi_l + t(psi_RE)
  # Put on probability scale 
  psi_mat <- exp(psi_l) / (1 + exp(psi_l))
  
# Calculate lambdas for each plot, level of drought
# (Note that loop takes several minutes to run)
  
  # Will create a list with 4 lambda matrices, one for each level of drought
  lambda_plots <- list() 
  # Will also extract elasticity values for adult female survival
  fem_surv_elas <- list()
  
  for (d in 1:4) {
    
    lambda_plots[[d]] <- matrix(NA, nrow = nrow(psi_mat), ncol = ncol(psi_mat))
    fem_surv_elas[[d]] <- matrix(NA, nrow = nrow(psi_mat), ncol = ncol(psi_mat))
    # Extract matrix of juvenile survival values for d level of drought
    phi1 <- phi1_mat[which(phi1_X$drought == drought4.z[d]),]
    # Extract matrix of female survival values for d level of drought
    phi2 <- phi2f_mat[which(phi2f_X$drought == drought4.z[d]),]
    
    for (i in 1:nrow(psi_mat)) {
      for (j in 1:ncol(psi_mat)) {
        
        proj.mat <- matrix(c(phi1[i, j] * (1 - psi_mat[i, j]), 0.32,
                             phi1[i, j] * psi_mat[i, j], phi2[i, j]),
                           nrow = 2, ncol = 2, byrow = TRUE)
        lambda_plots[[d]][i, j] <- eigen(proj.mat)$values[1]
        fem_surv_elas[[d]][i, j] <- elasticity(proj.mat)[2, 2]
        
      }
    }
  } 

# Summarize estimates:
  lambda_plots <- do.call(rbind, lambda_plots)
  ests_by_plot <- data.frame(plot = rep(plots$plot, 4),
                             drought = rep(drought4, each = nrow(plots))) %>%
    mutate(juv = apply(phi1_mat, 1, ctend),
           juv_lcl = apply(phi1_mat, 1, quantile, qprobs[1]),
           juv_ucl = apply(phi1_mat, 1, quantile, qprobs[2]),
           ad_fem = apply(phi2f_mat, 1, ctend),
           ad_fem_lcl = apply(phi2f_mat, 1, quantile, qprobs[1]),
           ad_fem_ucl = apply(phi2f_mat, 1, quantile, qprobs[2]),
           ad_male = apply(phi2m_mat, 1, ctend),
           ad_male_lcl = apply(phi2m_mat, 1, quantile, qprobs[1]),
           ad_male_ucl = apply(phi2m_mat, 1, quantile, qprobs[2]),
           trans = rep(apply(psi_mat, 1, ctend), 4),
           trans_lcl = rep(apply(psi_mat, 1, quantile, qprobs[1]), 4),
           trans_ucl = rep(apply(psi_mat, 1, quantile, qprobs[2]), 4),
           lambda = apply(lambda_plots, 1, ctend),
           lambda_lcl = apply(lambda_plots, 1, quantile, qprobs[1]),
           lambda_ucl = apply(lambda_plots, 1, quantile, qprobs[2]),
           probdecline = apply(lambda_plots, 1, function(x) sum(x < 1) / length(x)))
  
# Write to file
# write.csv(ests_by_plot, "output/plot-level-estimates.csv", row.names = FALSE)

# Summarize elasticity values for adult female survival (at PDSI = -1.08)
  mean(fem_surv_elas[[2]]) # mean = 0.77
  quantile(fem_surv_elas[[2]], probs = qprobs) # 95% CI = 0.58-0.94

#------------------------------------------------------------------------------# 
# Comparing demographic rates and lambda at each plot
#------------------------------------------------------------------------------#	  
# For now, focusing comparisons/correlations for PDSI = long-term mean (-1.08)
  
  ests_droughtmean <- ests_by_plot %>%
    filter(round(drought) == -1)

# 3-panel figure with plot-specific demographic rates vs lambda  
  
  barc <- "gray"
  barw <- 0.3
  
  juvlambda <- ggplot(ests_droughtmean, aes(x = juv, y = lambda)) +
    geom_errorbar(aes(xmin = juv_lcl, xmax = juv_ucl), 
                  col = barc, linewidth = barw) +
    geom_errorbar(aes(ymin = lambda_lcl, ymax = lambda_ucl), 
                  col = barc, linewidth = barw) +
    geom_point() +
    labs(x = "Juvenile survival", y = "Rate of popuation change") +
    theme_classic() +
    theme(text = element_text(size = 9),
          axis.text = element_text(size = 9))
  adlambda <- ggplot(ests_droughtmean, aes(x = ad_fem, y = lambda)) +
    geom_errorbar(aes(xmin = ad_fem_lcl, xmax = ad_fem_ucl), 
                  col = barc, linewidth = barw) +
    geom_errorbar(aes(ymin = lambda_lcl, ymax = lambda_ucl), 
                  col = barc, linewidth = barw) +
    geom_point() +
    labs(x = "Adult female survival", y = "Rate of popuation change") +
    scale_x_continuous(limits = c(0.89, 0.99), breaks = seq(0.9, 0.98, 0.04)) +
    theme_classic() +
    theme(text = element_text(size = 9),
          axis.text = element_text(size = 9),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
  translambda <- ggplot(ests_droughtmean, aes(x = trans, y = lambda)) +
    geom_errorbar(aes(xmin = trans_lcl, xmax = trans_ucl), 
                  col = barc, linewidth = barw) +
    geom_errorbar(aes(ymin = lambda_lcl, ymax = lambda_ucl), 
                  col = barc, linewidth = barw) +
    geom_point() +
    labs(x = "Transition rate", y = "Rate of popuation change") +
    theme_classic() +
    theme(text = element_text(size = 9),
          axis.text = element_text(size = 9),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
  demog_lambda <- plot_grid(juvlambda, adlambda, translambda, nrow = 1,
                            rel_widths = c(1.23, 1, 1))
  
  # ggsave("output/demog-v-lambda.jpg",
  #        demog_lambda,
  #        device = "jpeg",
  #        dpi = 600,
  #        width = 6.5,
  #        height = 3,
  #        units = "in")
  
# Correlations among demographic parameter estimates
  cor(ests_droughtmean[, c("juv", "ad_fem", "ad_male", "trans")])
  plot(ests_droughtmean[, c("juv", "ad_fem", "ad_male", "trans")])
  
# Correlations between vital rates and lambda
  # Adult survival (female)
  plot(lambda ~ ad_fem, data = ests_droughtmean)
  cor.test(ests_droughtmean$lambda, ests_droughtmean$ad_fem)
  # r = 0.05 (P = 0.84)
  
  # Juvenile survival
  plot(lambda ~ juv, data = ests_droughtmean)
  cor.test(ests_droughtmean$lambda, ests_droughtmean$juv)
  # r = 0.88 (P < 0.001)
  
  # Transition rate
  plot(lambda ~ trans, data = ests_droughtmean)
  cor.test(ests_droughtmean$lambda, ests_droughtmean$trans)
  # r = 0.53 (P = 0.03), though one plot (WM) seems to have a lot of leverage
  # with very high transition rate and lambda
  cor.test(ests_droughtmean$lambda[c(1:15, 17)], 
           ests_droughtmean$trans[c(1:15, 17)])
  # r = 0.32 (P = 0.23)
  
# Correlations between lambda values and latitude/longitude
# (Commented out because plot coordinates not publicly available)
  # coords <- read.csv(".../Plots_Coords.csv")
  # 
  # ests_by_plot <- ests_by_plot %>%
  #   left_join(coords[, c("plot", "lat", "long")], by = "plot")
  # 
  # cor.test(ests_by_plot$lambda[round(ests_by_plot$drought) == -1], 
  #          ests_by_plot$lat[round(ests_by_plot$drought) == -1])
  # # r = -0.05 (P = 0.86)
  # cor.test(ests_by_plot$lambda[round(ests_by_plot$drought) == -1], 
  #          ests_by_plot$long[round(ests_by_plot$drought) == -1])
  # # r = 0.26 (P = 0.32)
  
#------------------------------------------------------------------------------# 
# Estimates of mean demographic, detection rates across sampled populations
#------------------------------------------------------------------------------#	
# For adult survival, generating values for 2019-2020
# For survival, generating values at PDSI = -1.08 (mean across plots and years)
  
# Note: we could do this 2 ways: 1) use the intercept from the linear model for
# each rate, or, 2) calculate the rate for each plot (city and mean precip = 
# observed values, drought.z = 0, including random effect), average across plots 
# for each iteration, and then summarize the resulting distribution. 

# Option 1 would reflect a hypothetical plot that was at the mean distance from 
# a city and had mean annual precipitation.
# Option 2 would reflect the mean across sampled populations.  
# Going with option 2 (commented out approach using option 1 below)

# Option 1: use intercept from model for each rate
  # samples_df <- samples_df %>%
  #     mutate(phi1.mn = exp(beta.phi1) / (1 + exp(beta.phi1)),
  #            phi2.f = exp(beta.phi2) / (1 + exp(beta.phi2)),
  #            phi2.m = exp(beta.phi2 + b2.male) / (1 + exp(beta.phi2 + b2.male)),
  #            psi12.mn = exp(gamma.psi) / (1 + exp(gamma.psi)),
  #            p1.mn = exp(alpha.p1) / (1 + exp(alpha.p1)),
  #            p2.f = exp(alpha.p2) / (1 + exp(alpha.p2)),
  #            p2.m = exp(alpha.p2 + a2.male) / (1 + exp(alpha.p2 + a2.male)))
  # samples_means_int <- samples_df %>%
  #   select(phi1.mn, phi2.f, phi2.m, psi12.mn, p1.mn, p2.f, p2.m)
  # rate_means_int <- data.frame(mean = apply(samples_means_int, 2, mean),
  #                              median = apply(samples_means_int, 2, median),
  #                              lcl = apply(samples_means_int, 2, quantile, qprobs[1]),
  #                              ucl = apply(samples_means_int, 2, quantile, qprobs[2]))
  # rate_means_int$param <- rownames(rate_means_int)
  # rate_means_int <- rate_means_int %>%
  #   relocate(param, .before = mean)
  # row.names(rate_means_int) <- NULL
  # rate_means_int

# Option 2: calculate rate for each plot, average across plots for each 
# iteration, then summarize distribution
  # Grab samples for each demographic rate & plot, with PDSI = long-term mean (-1)
  phi1_mat_0 <- phi1_mat[18:34, ]
  phi2m_mat_0 <- phi2m_mat[18:34, ]
  phi2f_mat_0 <- phi2f_mat[18:34, ]
  psi_mat_0 <- psi_mat
  # For each iteration, average demographic rates across plots
  phi1_0 <- apply(phi1_mat_0, 2, mean)
  phi2m_0 <- apply(phi2m_mat_0, 2, mean)
  phi2f_0 <- apply(phi2f_mat_0, 2, mean)
  psi_0 <- apply(psi_mat, 2, mean)
  
  # Generate plot-specific detection probs for juveniles, with effort = mean
  p1_samples <- samples_df %>%
    select(alpha.p1, a1.precip)
  p1_X <- data.frame(intercept = 1,
                     mnprecip = precip.norm)
  p1_plots_l <- as.matrix(p1_X) %*% t(p1_samples)
  p1_plots <- exp(p1_plots_l) / (1 + exp(p1_plots_l))
  p1 <- apply(p1_plots, 2, mean)
  
  # Generate plot-specific detection probs for adult males, with effort = mean
  p2_samples <- samples_df %>%
    select(alpha.p2, a2.male, a2.precip)
  p2_X <- data.frame(intercept = 1,
                     male = 1,
                     mnprecip = precip.norm)
  p2m_plots_l <- as.matrix(p2_X) %*% t(p2_samples)
  p2m_plots <- exp(p2m_plots_l) / (1 + exp(p2m_plots_l))
  p2m <- apply(p2m_plots, 2, mean)
  
  # Generate plot-specific detection probs for adult females, with effort = mean
  p2_samples <- samples_df %>%
    select(alpha.p2, a2.male, a2.precip)
  p2_X <- p2_X %>% mutate(male = 0)
  p2f_plots_l <- as.matrix(p2_X) %*% t(p2_samples)
  p2f_plots <- exp(p2f_plots_l) / (1 + exp(p2f_plots_l))
  p2f <- apply(p2f_plots, 2, mean)
  
  samples_means <- cbind(phi1_0, phi2f_0, phi2m_0, psi_0, p1, p2f, p2m)
  
  rate_means <- data.frame(param = c("phi1.mn", "phi2.f", "phi2.m", 
                                     "psi12.mn", "p1.mn", "p2.f", "p2.m"),
                           mean = apply(samples_means, 2, mean),
                           median = apply(samples_means, 2, median),
                           lcl = apply(samples_means, 2, quantile, qprobs[1]),
                           ucl = apply(samples_means, 2, quantile, qprobs[2]))
  row.names(rate_means) <- NULL
  rate_means

# Calculate mean lambda (with mean PDSI) across sampled populations
# (assuming recruitment = 0.32 f/f/yr) 

  lambda_0 <- rep(NA, length = length(phi1_0))
  
  for (i in 1:length(phi1_0)) {
      proj.mat <- matrix(c(phi1_0[i] * (1 - psi_0[i]), 0.32,
                           phi1_0[i] * psi_0[i], phi2f_0[i]),
                         nrow = 2, ncol = 2, byrow = TRUE)
      lambda_0[i] <- eigen(proj.mat)$values[1]
  }
  
  rate_means_add <- data.frame(param = "lambda",
                               mean = mean(lambda_0),
                               median = median(lambda_0),
                               lcl = quantile(lambda_0, qprobs[1]),
                               ucl = quantile(lambda_0, qprobs[2]))
  rate_means <- rbind(rate_means, rate_means_add)
  
# Write to file
# write.csv(rate_means, "output/rate-means.csv", row.names = FALSE)

# Calculate mean probability of decline (mean across sampled plots)
  sum(lambda_0 < 1) / length(lambda_0)
  
#------------------------------------------------------------------------------# 
# Projected estimates of lambda under a range of drought conditions
#------------------------------------------------------------------------------#	
# For adult survival, generating values for 2019-2020
# For survival, generating values for 4 levels of drought: 
  # mean PDSI = -3, -1.08 (mean value across plots and years), 0, +3
# Will assume city and meanprecip = mean (0, when standardized)
# Will account for spatial random effects in all demographic parameters
# Will assume recruitment = 0.32 f/f/yr   

# Set a seed, since we'll be generating random values of spatial random effects
# for each iteration
set.seed(123)

# Juvenile survival
  # Generate a random effect from normal distribution with iteration-specific SDs
  phi1_e <- rnorm(nrow(samples_df), mean = 0, sd = samples_df$sigma.site.1)
  phi1_e <- matrix(rep(phi1_e, 4), nrow = 4, byrow = TRUE)
  # Calculate fixed effect of drought
  phi1_samples <- samples_df %>% select(beta.phi1, b1.drought)
  phi1_X <- data.frame(int = 1, drought.z = drought4.z)
  phi1_fixed <- as.matrix(phi1_X) %*% t(phi1_samples)
  # Combine fixed and random components
  phi1_l <- phi1_fixed + phi1_e
  # Put on probability scale
  phi1_mat <- exp(phi1_l) / (1 + exp(phi1_l))

# Adult male survival
  # Generate a random effect from normal distribution with iteration-specific SDs
  phi2_e <- rnorm(nrow(samples_df), mean = 0, sd = samples_df$sigma.site.2)
  phi2_e <- matrix(rep(phi2_e, 4), nrow = 4, byrow = TRUE)
  # Calculate fixed effect of drought
  phi2_samples <- samples_df %>%
    select(c(beta.phi2, b2.male, b2.drought, b2.trend, b2.trend2))
  phi2m_X <- data.frame(int = 1, male = 1, drought.z = drought4.z,
                        trend = tail(trend.z, 1), 
                        trend2 = tail(trend.z, 1) * tail(trend.z, 1))
  phi2m_fixed <- as.matrix(phi2m_X) %*% t(phi2_samples)
  # Combine fixed and random components
  phi2m_l <- phi2m_fixed + phi2_e
  # Put on probability scale
  phi2m_mat <- exp(phi2m_l) / (1 + exp(phi2m_l))
  
# Adult female survival
  # Calculate fixed effect of drought
  phi2f_X <- phi2m_X %>% mutate(male = 0)
  phi2f_fixed <- as.matrix(phi2f_X) %*% t(phi2_samples)
  # Combine fixed and random components
  phi2f_l <- phi2f_fixed + phi2_e
  # Put on probability scale
  phi2f_mat <- exp(phi2f_l) / (1 + exp(phi2f_l))
 
# Transition rate
  # Generate a random effect from normal distribution with iteration-specific SDs
  psi_e <- rnorm(nrow(samples_df), mean = 0, sd = samples_df$sigma.site.t)
  psi_e <- matrix(rep(psi_e, 4), nrow = 4, byrow = TRUE)
  # Calculate fixed effect of drought
  psi_samples <- samples_df %>% select(gamma.psi)
  psi_X <- data.frame(int = rep(1, 4))
  psi_fixed <- as.matrix(psi_X) %*% t(psi_samples)
  # Combine fixed and random components
  psi_l <- psi_fixed + psi_e
  # Put on probability scale
  psi_mat <- exp(psi_l) / (1 + exp(psi_l)) 
  
# Calculate lambdas for each level of drought

  lambda_overall <- matrix(NA, nrow = nrow(phi1_mat), ncol = ncol(phi1_mat))

  for (i in 1:nrow(phi1_mat)) {
    for (j in 1:ncol(phi1_mat)) {
      
      proj.mat <- matrix(c(phi1_mat[i, j] * (1 - psi_mat[i, j]), 0.32,
                           phi1_mat[i, j] * psi_mat[i, j], phi2f_mat[i, j]),
                         nrow = 2, ncol = 2, byrow = TRUE)
      lambda_overall[i, j] <- eigen(proj.mat)$values[1]
      
    }
  }

# Generate figure with overlapping histograms for lambda at PDSI = -3, 0, 3
  # Convert matrix to dataframe
  lambdas_df <- data.frame(PDSI = rep(drought4, each = ncol(lambda_overall)),
                           lambda = c(lambda_overall[1, ], lambda_overall[2, ],
                                      lambda_overall[3, ], lambda_overall[4,])) %>%
    mutate(PDSI = as.factor(PDSI))
  
  lambdas_df %>%
    group_by(PDSI) %>%
    summarize(mn = mean(lambda),
              md = median(lambda),
              min = min(lambda),
              max = max(lambda),
              sd = sd(lambda)) %>%
    data.frame()
  
  lambdas_fig <- ggplot(filter(lambdas_df, PDSI %in% c(-3, 0, 3))) + 
    geom_density(aes(x = lambda, group = PDSI, fill = PDSI), alpha = 0.25) +
    geom_vline(xintercept = 1, linetype = 2, col = "gray40") +
    labs(x = "Estimated rate of population change", y = "Density") +
    theme_classic() +
    theme(text = element_text(size = 8),
          axis.text = element_text(size = 8),
          legend.position = c(0.02, 0.98), 
          legend.justification = c(0, 1))
    
  # ggsave("output/lambda-distributions.jpg",
  #        lambdas_fig,
  #        device = "jpeg",
  #        dpi = 600,
  #        width = 3,
  #        height = 3,
  #        units = "in")

# Summarize everything:
  ests_overall <- data.frame(drought = drought4) %>%
    mutate(juv = apply(phi1_mat, 1, ctend),
           juv_lcl = apply(phi1_mat, 1, quantile, qprobs[1]),
           juv_ucl = apply(phi1_mat, 1, quantile, qprobs[2]),
           ad_fem = apply(phi2f_mat, 1, ctend),
           ad_fem_lcl = apply(phi2f_mat, 1, quantile, qprobs[1]),
           ad_fem_ucl = apply(phi2f_mat, 1, quantile, qprobs[2]),
           ad_male = apply(phi2m_mat, 1, ctend),
           ad_male_lcl = apply(phi2m_mat, 1, quantile, qprobs[1]),
           ad_male_ucl = apply(phi2m_mat, 1, quantile, qprobs[2]),
           trans = apply(psi_mat, 1, ctend),
           trans_lcl = apply(psi_mat, 1, quantile, qprobs[1]),
           trans_ucl = apply(psi_mat, 1, quantile, qprobs[2]),
           lambda = apply(lambda_overall, 1, ctend),
           lambda_lcl = apply(lambda_overall, 1, quantile, qprobs[1]),
           lambda_ucl = apply(lambda_overall, 1, quantile, qprobs[2]),
           probdecline = apply(lambda_overall, 1, function(x) sum(x < 1) / length(x)))
  
# Write to file
# write.csv(ests_overall, "output/statewide-projections.csv", row.names = FALSE)
