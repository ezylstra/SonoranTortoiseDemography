### Analysis of capture-recapture data for Sonoran desert tortoises from long-term monitoring plots across Arizona. 

## Code
[multistate-model.R](multistate-model.R): loads capture-recapture and covariate data, runs a Bayesian multi-state model using the nimble package, and summarizes model results. 

## Data
1. [CapRecapData.csv](data/CapRecapData.csv): Capture-recapture data for all tortoises that were marked between 1987 and 2020 at each of 17 long-term monitoring plots.
    - plot: 2-letter code for each of the 17 plots
    - yr: capture year
    - obsdate: capture date (YYYY-MM-DD format)
    - tort: unique alphanumeric ID for each tortoise
    - sex: Sex of tortoise, where 1 = female, 2 = male, and 3 = unknown (most of which are juveniles). 
  - MCL: Maximum carapace length
2. [Plots_NoCoords.csv](data/Plots_NoCoords.csv): Covariate data for each plot.
    - plot: 2-letter code for each of the 17 plots
    - climate: climate division
    - city.km: distance (km) from each plot to the center of the nearest city with ???10,000 people (based on 2010 US Census data). 
    - pptnorms.mm: Mean annual precipitation (mm) at each plot (data from PRISM)
3. [PDSI.csv](data/PDSI.csv): Monthly measures of Palmer drought severity index (PDSI) for climate divisions in Arizona from 1985 to 2020.
    - div: climate division
    - yr: year
    - mon: month
    - pdsi: PDSI value (data from the National Climatic Data Center)
4. [Precip_Monthly.csv](data/Precip_Monthly.csv): Monthly precipitation totals for each plot from 1985 to 2020. 
    - plot: 2-letter code for each of the 17 plots
    - yr.m: year and month (YYYY-MM format)
    - ppt: Cumulative precipitation (mm; data from PRISM)
5. [Surveys.csv](data/Surveys.csv): Summaries of survey effort at each plot from 1987 to 2020. 
    - plot: Name of plot
    - code: 2-letter code for each of the 17 plots
    - area.sqmi: area (sq. km) of each plot
    - yr: Survey year
    - persondays: number of person days spent searching a plot for tortoises in a given year.