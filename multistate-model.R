################################################################################
# Multi-state model to estimate demographic rates of Sonoran desert tortoises in 
# Arizona, 1987-2020

# ER Zylstra
# Last updated: 16 July 2023
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
  # Ctortmodel <-compileNimble(tortmodel)
  # Ctortmcmc <- compileNimble(tortmcmc, project = tortmodel)

# MCMC settings, parameters, initial values  
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
  MCMCplot(samples,
           params = "all", # excl = ""
           ci = c(50, 90))  
  
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
  #ctend <- median

# 90% or 95% credible intervals
  #qprobs <- c(0.05,0.95)
  qprobs <- c(0.025,0.975) 
  
# Colors (n = 2)
  mycol <- col2rgb(c('salmon3','steelblue4'))
  col1 <- rgb(mycol[1,1],mycol[2,1],mycol[3,1],alpha=255,max=255)
  col1p <- rgb(mycol[1,1],mycol[2,1],mycol[3,1],alpha=0.2*255,max=255)
  col2 <- rgb(mycol[1,2],mycol[2,2],mycol[3,2],alpha=255,max=255)
  col2p <- rgb(mycol[1,2],mycol[2,2],mycol[3,2],alpha=0.2*255,max=255)
  
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
  mean.f.dry <- apply(pred[1:100,], 1, ctend)
  mean.f.avg <- apply(pred[101:200,], 1, ctend)
  mean.f.wet <- apply(pred[201:300,], 1, ctend)
  mean.m.dry <- apply(pred[301:400,], 1, ctend)
  mean.m.avg <- apply(pred[401:500,], 1, ctend)
  mean.m.wet <- apply(pred[501:600,], 1, ctend)
  ci.f.dry <- apply(pred[1:100,], 1, quantile, probs = qprobs)
  ci.f.avg <- apply(pred[101:200,], 1, quantile, probs = qprobs)
  ci.f.wet <- apply(pred[201:300,], 1, quantile, probs = qprobs)
  ci.m.dry <- apply(pred[301:400,], 1, quantile, probs = qprobs)
  ci.m.avg <- apply(pred[401:500,], 1, quantile, probs = qprobs)
  ci.m.wet <- apply(pred[501:600,], 1, quantile, probs = qprobs)
  
  plotx <- predx[1:100, "drought"] * pdsi.24.sd + pdsi.24.mn
  
  # Figure with M/F adult survival at extreme mnprecip values (wet/dry plots)
  ### Want to do this with ggplot (and save as pdf and jpg) ###
  
  par(mar = c(2.5, 3.5, 0.5, 0.6), cex = 0.8)
  plot(mean.f.dry ~ plotx, type = "l", lty = 1, xaxt= "n", yaxt= "n", xlab = "", 
       ylab = "",  ylim = c(0.78, 1), bty = "n", yaxs = "i", col = col1)
  axis(1, at = c(par("usr")[1], par("usr")[2]), tck = F, labels = F)
  axis(1, at = seq(-4, 4, by = 2), labels = seq(-4, 4, by = 2), tcl = -0.25,
       mgp = c(1.5, 0.4, 0))
  axis(2, at = c(par("usr")[3], par("usr")[4]), tck = F, labels = F)
  axis(2, at = seq(0.8, 1, by = 0.05), tcl = -0.25, las = 1, mgp = c(1.5, 0.5, 0),
       labels = c("0.80", "0.85", "0.90", "0.95", "1.00"))
  lines(mean.f.wet ~ plotx, type = "l", lty = 2, col = col1)
  lines(mean.m.dry ~ plotx, type = "l", lty = 1, col = col2)
  lines(mean.m.wet ~ plotx, type = "l", lty = 2, col = col2)
  arrows(x0 = 0, x1 = 0, y0 = 0.68, y1 = 1, length = 0, col = "gray50", lty = 3)
  mtext("Adult survival", side = 2, las = 0, line = 2.5, cex = 0.8)
  mtext("PDSI (24-month)", side = 1, line = 1.5, cex = 0.8)
  legend("bottomright", c("Arid:F", "Arid:M'", "Semiarid:F", "Semiarid:M"),
         lty = c(1, 1, 2, 2), col = c(col1, col2, col1, col2), bty = "n")

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
  mean.j.dry <- apply(pred1[1:100,], 1, ctend)
  mean.j.avg <- apply(pred1[101:200,], 1, ctend)
  mean.j.wet <- apply(pred1[201:300,], 1, ctend)
  ci.j.dry <- apply(pred1[1:100,], 1, quantile, probs = qprobs)
  ci.j.avg <- apply(pred1[101:200,], 1, quantile, probs = qprobs)
  ci.j.wet <- apply(pred1[201:300,], 1, quantile, probs = qprobs)
  
  # Figure with juvenile survival at extreme mnprecip values (wet/dry plots)
  ### Want to do this with ggplot (and save as pdf and jpg) ###
  
  par(mar = c(2.5, 3.5, 0.5, 0.6), cex = 0.8)
  plot(mean.j.dry ~ plotx, type = "l", lty = 1, xaxt= "n", yaxt= "n", xlab = "", 
       ylab = "",  ylim = c(0.45, 1), bty = "n", yaxs = "i", col = "black")
  axis(1, at = c(par("usr")[1], par("usr")[2]), tck = F, labels = F)
  axis(1, at = seq(-4, 4, by = 2), labels = seq(-4, 4, by = 2), tcl = -0.25,
       mgp = c(1.5, 0.4, 0))
  axis(2, at = c(par("usr")[3], par("usr")[4]), tck = F, labels = F)
  axis(2, at = seq(0.5, 1, by = 0.1), tcl = -0.25, las = 1, mgp = c(1.5, 0.5, 0),
       labels = c("0.50", "0.60", "0.70", "0.80", "0.90", "1.00"))
  lines(mean.j.wet ~ plotx, type = "l", lty = 2, col = "black")
  arrows(x0 = 0, x1 = 0, y0 = 0.68, y1 = 1, length = 0, col = "gray50", lty = 3)
  mtext("Juvenile survival", side = 2, las = 0, line = 2.5, cex = 0.8)
  mtext("PDSI (24-month)", side = 1, line = 1.5, cex = 0.8)
  legend("bottomright", c("Arid", "Semiarid"),
         lty = c(1, 2), col = "black", bty = "n")

# Adults and juveniles
  # Stacked figure with adult, juvenile survival at extreme mnprecip values
  ### Want to do this with ggplot (and save as pdf and jpg) ###
  
#------------------------------------------------------------------------------# 
# Temporal trends in adult survival?
#------------------------------------------------------------------------------#	  

# For an average site (mean distance from city, annual precip) with PDSI = 0
phi2t <- samples_df %>% 
  select(beta.phi2, b2.male, b2.trend, b2.trend2) %>%
  as.matrix()
  
predtx <- cbind(int = 1, male = rep(0:1, each = 33), trend = rep(trend.z, 2),
                trend2 = rep(trend.z2, 2))
predtl <- predtx %*% t(phi2t)
predt <- exp(predtl) / (1 + exp(predtl))  
  
predt.both <- data.frame(endyr = 1988:2020,
                         interval = rep(paste(1987:2019, 1988:2020, sep = "-")))
predt.both$female <- round(apply(predt[1:33,], 1, ctend), 3)
predt.both$female.lcl <- round(apply(predt[1:33,], 1, quantile, probs = qprobs[1]), 3)
predt.both$female.ucl <- round(apply(predt[1:33,], 1, quantile, probs = qprobs[2]), 3)
predt.both$male <- round(apply(predt[34:66,], 1, ctend), 3)
predt.both$male.lcl <- round(apply(predt[34:66,], 1, quantile, probs = qprobs[1]), 3)
predt.both$male.ucl <- round(apply(predt[34:66,], 1, quantile, probs = qprobs[2]), 3)   

#------------------------------------------------------------------------------# 
# Temporal trend in PDSI values
#------------------------------------------------------------------------------#	

#------------------------------------------------------------------------------# 
# Plot-specific estimates of demographic rates
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------# 
# Population growth rates
#------------------------------------------------------------------------------#	

