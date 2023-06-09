################################################################################
# Multi-state model for Sonoran desert tortoises in Arizona, 1987-2020

# ER Zylstra
# Last updated: 8 June 2023
################################################################################

library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)

rm(list=ls())

#------------------------------------------------------------------------------# 
# Load data
#------------------------------------------------------------------------------#

cr <- read.csv("data/CapRecapData.csv", header = TRUE)
surveys <- read.csv("data/Surveys.csv", header = TRUE)
pdsi <- read.csv("data/PDSI.csv", header = TRUE)
precip.norms <- read.csv("data/Precip_norms.csv", header = TRUE)
precip <- read.csv("data/Precip_Monthly.csv", header = TRUE)
city <- read.csv("data/DistToCity.csv", header = TRUE)

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
  
# Create a wide version of the persondays data
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
  
# Create a wide version of data, where 1/0 indicates when surveys were done
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
  # seven individual that were captured as adults but not assigned sex
  
# Retain only a single capture each year (using MCL at first capture) and 
# identify stage (juvenile = 1, adult = 2)
  cr.yr <- cr %>%
    mutate(obsdate = ymd(obsdate)) %>%
    group_by(plot, tort, sex, yr) %>%
    summarize(mcl = MCL[1],
              .groups = "keep") %>%
    mutate(stage = ifelse(mcl < 180, 1, 2)) %>%
    data.frame()

# Identify how many tortoises were caught at least once each year a plot was
# surveyed
  plot.yr <- cr.yr %>%
    group_by(plot, yr) %>%
    summarize(n.torts = length(tort),
              n.adultm = sum(sex == 2 & mcl > 179),
              n.adultf = sum(sex == 1 & mcl > 179),
              n.adultu = sum(sex == 3 & mcl > 179),
              n.juv = sum(mcl < 180),
              .groups = "keep") %>%
    data.frame()
  
# Identify how many tortoises were caught ever at each plot
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
  # Proportion of juveniles that were subsequently captured as adults (91)
  sum(firstlast$state1 == 1 & firstlast$state2 == 2)/sum(firstlast$state1 == 1)
  
#------------------------------------------------------------------------------#  
# Format covariates
#------------------------------------------------------------------------------#

# Formatting individual covariates
  sex <- ch$sex
  sex[sex == 3] <- NA
  male.ind <- sex - 1
  
# Format site covariate: distance to nearest "major" city (pop > 10,000)
  city <- city %>%
    rename(plot = code) %>%
    arrange(plot)
  # Check that plots are in the same order they're found in capture histories
  # all.equal(unique(ch$plot), unique(city$plot))
  dist.mn <- mean(city$dist.km)
  dist.sd <- sd(city$dist.km)
  distance <- (city$dist.km - dist.mn)/dist.sd

  