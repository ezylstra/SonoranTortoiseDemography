################################################################################
# Multi-state model for Sonoran desert tortoises in Arizona, 1987-2020

# ER Zylstra
# Last updated: 9 June 2023
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
  
  