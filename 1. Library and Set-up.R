
# Library Set up File

setwd("C:/Users/akg9010/Dropbox/cvd_hurricane")

#Install INLA package
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

#Install libraries
library(spdep)
library(readxl)
library(sf)
library(pscl) # Zero inflated models
library(INLAspacetime)
library(INLA)
library(tidyverse)
library(kableExtra) #Better tables using Kable
library(SpatialEpi)
library(ggplot2)
library(rgeoda)
library(CARBayes)
library(CARBayesST)
library(fabricatr) # splitting quantiles
library(spdep)
library(sp)
library(FRK)
library(MatchIt)

# Reading files
zctas <- read_sf("Data files/US_ZCTAs/tl_2020_us_zcta510.shp") # US ZCTA shapefile
all_data_quarters <- read_excel("Data files/Q_CVD_DIED_COVR_adjCHRLSNADI.xlsx") # ZCTA-level CVD data from Medicare claims
nyczip <- read.csv("Data files/US_ZCTAs//nyc-zip-codes.csv") # NYC ZCTA indicator list

# Processing ZCTA data
colnames(zctas)[1] <- "indx_zip2" # Change name of variable for variable
nyczip$NYC <- "1"  # Adding NYC indicator 
colnames(nyczip)[3] <- "indx_zip2" # change name of variable to ensure merging complete
nyczip$indx_zip2 <- as.character(nyczip$indx_zip2) # change to character

# Merge the shape and the data files
all_data_quarters_wshape <- merge(zctas, all_data_quarters, by = "indx_zip2", all.y = TRUE)
all_data_quarters_wshape2 <- merge(nyczip, all_data_quarters_wshape, by = "indx_zip2", all.y = TRUE)

# Further wrangling in preparation of analysis

all_data_quarters_wshape2$state_flood <- ifelse(all_data_quarters_wshape2$NYC %in% "1", 
                                                "NYC", all_data_quarters_wshape2$state_flood) # Change state_flood variable to include NYC, other_NY, CT, and NJ
all_data_quarters_wshape2$hurricane <- 0 # Defining when hurricane occurred in the dataset by dichotomizing the hurricane variable (before and after)
all_data_quarters_wshape2$hurricane[all_data_quarters_wshape2$QUARTER > 10] <- 1
all_data_quarters_wshape2$tenperc <- ifelse(all_data_quarters_wshape2$stormsurgeextent < 0.1, 0, 1) # Change the definition of flooding as sensitivity analysis (exposure = >= 10%)
all_data_quarters_wshape2$median_hshld_incm <- as.numeric(all_data_quarters_wshape2$median_hshld_incm) # Change median income to numerical variable
all_data_quarters_wshape2 <- all_data_quarters_wshape2 %>% # Modify median household income variable by normalizing data
  mutate(median_hshld_incm2 = (median_hshld_incm - mean(median_hshld_incm))/sd(median_hshld_incm))
all_data_quarters_wshape2$idarea <- as.numeric(as.factor(all_data_quarters_wshape2$indx_zip2)) # Create an alternative variable for area and time - required to run the random effects components of the INLA models
all_data_quarters_wshape2$idarea1 <- all_data_quarters_wshape2$idarea
all_data_quarters_wshape2$idtime <- 1 + all_data_quarters_wshape2$year - min(all_data_quarters_wshape2$year)

# Saving
saveRDS(all_data_quarters_wshape2, "cvd_hurricane. .rds")