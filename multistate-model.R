################################################################################
# Multi-state model to estimate demographic rates of Sonoran desert tortoises in 
# Arizona, 1987-2020

# ER Zylstra
# Last updated: 4 July 2023
################################################################################

library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(nimble)
library(MCMCvis)

rm(list=ls())

#------------------------------------------------------------------------------# 
# Load data
#------------------------------------------------------------------------------#

cr <- read.csv("data/CapRecapData.csv", header = TRUE)
plots <- read.csv("data/Plots.csv", header = TRUE)
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

# Identify how many unique tortoises were detected at least once during once 
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
  
  # First and last known state of each tortoise
  firstlast <- ch[, 1:3]
  firstlast$state1 <- unlist(apply(ch.mat, 1, 
                                   function(x) head(x[!is.na(x) & x < 3], 1)))
  firstlast$state2 <- unlist(apply(ch.mat, 1, 
                                   function(x) tail(x[!is.na(x) & x < 3], 1)))
  # Check that there are no backward transitions (ad -> juv)
  count(firstlast, state1, state2)
  # Proportion of individuals first captured as juvenile (n = 569)
  sum(firstlast$state1 == 1)/nrow(firstlast)
  # Proportion of individuals first captured as adult (n = 1466)
  sum(firstlast$state1 == 2)/nrow(firstlast)
  # Proportion of juveniles that were subsequently captured as adults (n = 91)
  sum(firstlast$state1 == 1 & firstlast$state2 == 2)/sum(firstlast$state1 == 1)
  
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
# Functions to create initial values for JAGS
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

      sigma.site.2 ~ T(dt(0, pow(2.5, -2), 1), 0, )
      tau.site.2 <- 1 / (sigma.site.2 * sigma.site.2)

      for (p in 1:nplots) {
        e.site.2[p] ~ dnorm(0, tau.site.2)
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
                               b1.int * mn.precip[plot[i]] * drought[plot[i], t]

          # Transition probability
          logit(psi12[i, t]) <- gamma.psi + c.mnprecip * mn.precip[plot[i]]

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
          # First index = states at time t-1, last index = states at time t
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
          # First index = states at time t, last index = detection type at time t
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
              "omega", "sigma.site.2", "e.site.2")
  
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
                            sigma.site.2 = runif(1, 0, 3),
                            male = ifelse(is.na(male.ind), 1, NA),
                            z = ch.init(as.matrix(ch.mat), first1, first2))}
  
# Separate constants and data for NIMBLE
  tortconstants <- tortdata
  tortconstants[c("male", "y")] <- NULL

# Create model object (~ 8 min)
  tortmodel <- nimbleModel(code = tortcode, constants = tortconstants, 
                           calculate = FALSE)

# Set data and inits
  tortmodel$setData(list(y = as.matrix(ch.mat),
                         male = male.ind))
  set.seed(123)
  tortmodel$setInits(inits())

# Build MCMC (~ 53 min)
  start_build <- Sys.time()
  tortmcmc <- buildMCMC(tortmodel,
                        monitors = params)
  end_build <- Sys.time()
  end_build - start_build

# Compile the model and MCMC (~ 26 min)
  Ctortmodel <-compileNimble(tortmodel)
  Ctortmcmc <- compileNimble(tortmcmc, project = tortmodel)

# MCMC settings, parameters, initial values  
  n.chains <- 3
  n.iter <- 30000
  n.burn <- 5000
  n.thin <- 15
  ni.tot <- n.iter + n.burn
    
# Run the MCMC and extract the samples (~13.4 hrs)
  start_run <- Sys.time()
  samples <- runMCMC(
    Ctortmcmc,
    nchains = n.chains,
    niter = ni.tot,
    nburnin = n.burn,
    thin = n.thin,
    samplesAsCodaMCMC = TRUE
  )
  end_run <- Sys.time()
  end_run - start_run

# Save samples
  saveRDS(samples, "MS-samples-6000.rds")
  
# Produce summary table, look at trace & density plots
  MCMCsummary(samples, 
              round = 2, 
              params = "all", 
              probs = c(0.025, 0.975))
  MCMCtrace(samples,
            params = "all",
            pdf = TRUE,
            open_pdf = FALSE)
  MCMCplot(samples,
           params = "all", # excl = ""
           ci = c(50, 90))

#------------------------------------------------------------------------------#  
# Post-processing
#------------------------------------------------------------------------------#    
  
# Load samples from previous run (if needed)
  samples <- readRDS("MS-samples-6000.rds")

# Create matrix with samples  
  samples_mat <- MCMCchains(samples, 
                            params = "all",
                            mcmc.list = FALSE)
  
# Calculate derived parameters
  samples_df <- data.frame(samples_mat) %>%
    mutate(psi12.mn = exp(gamma.psi) / (1 + exp(gamma.psi)),
           phi1.mn = exp(beta.phi1) / (1 + exp(beta.phi1)),
           p1.mn = exp(alpha.p1) / (1 + exp(alpha.p1)),
           phi2.f = exp(beta.phi2) / (1 + exp(beta.phi2)),
           phi2.m = exp(beta.phi2 + b2.male) / (1 + exp(beta.phi2 + b2.male)),
           p2.f = exp(alpha.p2) / (1 + exp(alpha.p2)),
           p2.m = exp(alpha.p2 + a2.male) / (1 + exp(alpha.p2 + a2.male)))
  
  
  
  
  #-- The following derived parameters had been in the original JAGS model
  #-- If we add them back in, need to add names to params
  # logit(psi12.mn) <- gamma.psi
  # logit(phi1.mn) <- beta.phi1
  # logit(p1.mn) <- alpha.p1
  # 
  # phi2.f <- exp(beta.phi2) / (1 + exp(beta.phi2))
  # phi2.m <- exp(beta.phi2 + b2.male) / (1 + exp(beta.phi2 + b2.male))
  # p2.f <- exp(alpha.p2) / (1 + exp(alpha.p2))
  # p2.m <- exp(alpha.p2 + a2.male) / (1 + exp(alpha.p2 + a2.male))
  
  # Run model with JAGS
    # fit.ms <- jags(data = tortdata, inits = inits, parameters.to.save = params,
    #                model.file="MS_siteRE_trend.txt",
    #                n.chains = n.chains, n.adapt = n.adapt, n.iter = ni.tot, 
    #                n.burnin = n.burn, n.thin = n.thin,
    #                parallel = TRUE, n.cores = 3, DIC = FALSE) 
    # 
    # # load(file.choose())
    # print(fit.ms, digits = 2)	
    # 
    # #Create a matrix of posterior samples
    # out <- fit.ms$samples
    # comb <- combine.mcmc(out)
    # phi1.s <- comb[ ,c("beta.phi1", colnames(comb)[grep("b1.", colnames(comb))])]
    # phi2.s <- comb[ ,c("beta.phi2", colnames(comb)[grep("b2.", colnames(comb))])]
    # phi2RE.s <- comb[ ,grep("e.site.2", colnames(comb))]
    # psi12.s <- comb[ ,c("gamma.psi", "c.mnprecip")]
  