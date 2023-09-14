################################################################################
# Multi-state model to estimate demographic rates of Sonoran desert tortoises in 
# Arizona, 1987-2020

# ER Zylstra
# Last updated: 14 September 2023
################################################################################

library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(nimble)
library(MCMCvis)
library(ggplot2)
library(cowplot)

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
  # tortmodel <- nimbleModel(code = tortcode, constants = tortconstants, 
  #                          calculate = FALSE)

# Set data and inits
  # tortmodel$setData(list(y = as.matrix(ch.mat),
  #                        male = male.ind))
  # set.seed(123)
  # tortmodel$setInits(inits())

# Build MCMC (~ 53 min)
  # tortmcmc <- buildMCMC(tortmodel,
  #                       monitors = params)

# Compile the model and MCMC (~ 26 min)
  # Ctortmodel <- compileNimble(tortmodel)
  # Ctortmcmc <- compileNimble(tortmcmc, project = tortmodel)

# MCMC settings
  n.chains <- 3
  n.iter <- 30000
  n.burn <- 5000
  n.thin <- 15
  ni.tot <- n.iter + n.burn
    
# Run the MCMC and extract the samples (~13.4 hrs)
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
                 mnprecip = rep(rep(mnprecip3.z, 2), each=100),
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
                                     "Semiarid", "Mesic")),
           group = paste(regime, sex, sep = ":"),
           pdsi = drought * pdsi.24.sd + pdsi.24.mn)

  # Set colors, linetypes for figure (types: 1 = solid, 2 = dashed, 3 = dotted) 
  groups <- data.frame(unique(adsurv_df[, c("group", "sex", "regime")])) %>%
    mutate(col = ifelse(sex == "F", "salmon3", "steelblue4"),
           linetype = ifelse(regime == "Arid", 1,
                             ifelse(regime == "Semiarid", 3, 2)))
  linewidth <- 0.3
  
  # Just using estimates for M/F at dry and wet sites
  adsurv_df4 <- filter(adsurv_df, regime != "Semiarid") %>%
    mutate(group = as.factor(group))
  groups4 <- filter(groups, regime != "Semiarid") %>%
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
    theme(text = element_text(size = 8),
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
                                  "Semiarid", "Mesic")),
           pdsi = drought * pdsi.24.sd + pdsi.24.mn)  
  
  # Set colors, linetypes for figure (types: 1 = solid, 2 = dashed, 3 = dotted) 
  groups_juv <- data.frame(regime = unique(juvsurv_df$regime)) %>%
    mutate(col = "black",
           linetype = ifelse(regime == "Arid", 1,
                             ifelse(regime == "Semiarid", 3, 2)))
  linewidth <- 0.3
  
  # Just using estimates  at dry and wet sites
  juvsurv_df2 <- filter(juvsurv_df, regime != "Semiarid") %>%
    mutate(regime = as.factor(regime))
  groups_juv2 <- filter(groups_juv, regime != "Semiarid") %>%
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
    theme(text = element_text(size = 8),
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
    theme(text = element_text(size = 8),
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

# For an average site (mean distance from city, annual precip) with PDSI = 0
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
  theme(text = element_text(size = 8),
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
                 'goldenrod4','salmon4'))
linetype <- 1
linewidth <- 0.3

pdsi_plot <- ggplot() +
  geom_line(data = pdsi24t, aes(x = yr, y = pdsi.24, group = div, col = div),
            linewidth = linewidth, alpha = 0.6, show.legend = NA) +
  geom_point(data = pdsi24t, aes(x = yr, y = pdsi.24, col = div, fill = div), 
             shape = 21, size = 1, alpha = 0.6) +
  geom_line(data = pdsi_df, aes(x = yr, y = central), show.legend = NA) +
  geom_ribbon(data = pdsi_df, aes(x = yr, ymin = lcl, ymax = ucl), alpha = 0.2,
              show.legend = NA) +
  labs(x = "Year", y = "PDSI (24-month)") +
  scale_x_continuous(limits = c(1988, 2020), expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(-4, 4.7), breaks = seq(-4, 4, by = 2)) +
  scale_color_manual(values = divs$col, name = "Division") +
  scale_fill_manual(values = divs$col, name = "Division") +
  theme_classic() +
  theme(text = element_text(size = 8),
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
# Plot-specific estimates of demographic rates
#------------------------------------------------------------------------------#
# Generating values for 2019-2020, with drought = 0

# Juvenile survival
  phi1p <- samples_df %>%
    select(beta.phi1, b1.city, b1.mnprecip)
  predp1x <- cbind(int = 1, city = city, mnprecip = precip.norm)
  predp1l <- predp1x %*% t(phi1p)
  predp1 <- exp(predp1l)/(1 + exp(predp1l))  
  plotests <- data.frame(plot = plots$plot) %>%
    mutate(juv = round(apply(predp1, 1, ctend), 2),
           juv_lcl = round(apply(predp1, 1, quantile, qprobs[1]), 2),
           juv_ucl = round(apply(predp1, 1, quantile, qprobs[2]), 2))

# Adult survival
  phi2p <- samples_df %>%
    select(beta.phi2, b2.male, b2.city, b2.mnprecip, b2.trend, b2.trend2)
  phi2RE <- select(samples_df, contains("e.site.2"))

  predpx <- cbind(int = 1, male = rep(0:1, each = 17), city = rep(city, 2),
                  mnprecip = rep(precip.norm, 2), trend = tail(trend.z, 1),
                  trend2 = tail(trend.z2, 1))
  predpl <- predpx %*% t(phi2p)
  RE_mat <- rbind(t(phi2RE), t(phi2RE))
  predpl <- predpl + RE_mat
  predp <- exp(predpl)/(1 + exp(predpl)) 
  
  plotests <- plotests %>%
    mutate(ad_fem = round(apply(predp[1:17, ], 1, ctend), 2),
           ad_fem_lcl = round(apply(predp[1:17, ], 1, quantile, qprobs[1]), 2),
           ad_fem_ucl = round(apply(predp[1:17, ], 1, quantile, qprobs[2]), 2),
           ad_male = round(apply(predp[18:34, ], 1, ctend), 2),
           ad_male_lcl = round(apply(predp[18:34, ], 1, quantile, qprobs[1]), 2),
           ad_male_ucl = round(apply(predp[18:34, ], 1, quantile, qprobs[2]), 2))
  
  # M/F survival in 2019-2020, across all plots
  phi2pA <- select(samples_df, c(beta.phi2, b2.male, b2.trend, b2.trend2))
  predpxA <- cbind(int = 1, male = 0:1, 
                   trend = tail(trend.z, 1), trend2 = tail(trend.z2, 1))
  predplA <- predpxA %*% t(phi2pA)
  predpA <- exp(predplA)/(1 + exp(predplA)) 
  adsurvivalests <- data.frame(sex = c("F", "M")) %>%
    mutate(mn = round(apply(predpA, 1, ctend), 2),
           lcl = round(apply(predpA, 1, quantile, qprobs[1]), 2),
           ucl = round(apply(predpA, 1, quantile, qprobs[2]), 2))
  adsurvivalests

# Transition rates
  psi12p <- select(samples_df, c(gamma.psi, c.mnprecip))
  predpsix <- cbind(int = 1, mnprecip = precip.norm)
  predpsil <- predpsix %*% t(psi12p)
  predpsi <- exp(predpsil) / (1 + exp(predpsil))
  plotests <- plotests %>%
    mutate(trans = round(apply(predpsi, 1, ctend), 2),
           trans_lcl = round(apply(predpsi, 1, quantile, qprobs[1]), 2),
           trans_ucl = round(apply(predpsi, 1, quantile, qprobs[2]), 2))

# Add overall estimates to bottom of table
plotests_add <- data.frame(plot = "Overall",
                           juv = ctend(samples_df$phi1.mn),
                           juv_lcl = quantile(samples_df$phi1.mn, qprobs[1]),
                           juv_ucl = quantile(samples_df$phi1.mn, qprobs[2]),
                           ad_fem = adsurvivalests$mn[adsurvivalests$sex == "F"],
                           ad_fem_lcl = adsurvivalests$lcl[adsurvivalests$sex == "F"],
                           ad_fem_ucl = adsurvivalests$ucl[adsurvivalests$sex == "F"],
                           ad_male = adsurvivalests$mn[adsurvivalests$sex == "M"],
                           ad_male_lcl = adsurvivalests$lcl[adsurvivalests$sex == "M"],
                           ad_male_ucl = adsurvivalests$ucl[adsurvivalests$sex == "M"],
                           trans = ctend(samples_df$psi12.mn),
                           trans_lcl = quantile(samples_df$psi12.mn, qprobs[1]),
                           trans_ucl = quantile(samples_df$psi12.mn, qprobs[2]),
                           row.names = NULL) 
plotests_add <- plotests_add %>%
  mutate(across(juv:trans_ucl, function(x) round(x, 2)))
plotests <- rbind(plotests, plotests_add)
plotests

#------------------------------------------------------------------------------# 
# Population growth rates
#------------------------------------------------------------------------------#	


