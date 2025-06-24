getwd()
######## Load library #########
library(tidyverse)
library(reshape2)
library(lubridate)
library(psych)
library(gganimate)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(vtable)
source("function_indirect.CO2.balance.method.R")

######## Import Data #########
#Load processed gas datasets
LUFA_FTIR_long <- read.csv("20250408-15_long_LUFA_FTIR.csv")
LUFA_FTIR_long$DATE.TIME <- as.POSIXct(LUFA_FTIR_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

ANECO_FTIR_long <- read.csv("20250408-15_long_ANECO_FTIR.csv")
ANECO_FTIR_long$DATE.TIME <- as.POSIXct(ANECO_FTIR_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

MBBM_FTIR_long <- read.csv("20250408-15_long_MBBM_FTIR.csv")
MBBM_FTIR_long$DATE.TIME <- as.POSIXct(MBBM_FTIR_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

ATB_FTIR_long <- read.csv("20250408-15_long_ATB_FTIR.1.csv")
ATB_FTIR_long$DATE.TIME <- as.POSIXct(ATB_FTIR_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

ATB_CRDS_long <- read.csv("20250408-15_long_ATB_CRDS.P8.csv")
ATB_CRDS_long$DATE.TIME <- as.POSIXct(ATB_CRDS_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

LUFA_CRDS_long <- read.csv("20250408-15_long_LUFA_CRDS.P8.csv")
LUFA_CRDS_long$DATE.TIME <- as.POSIXct(LUFA_CRDS_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

UB_CRDS_long <- read.csv("20250408-15_long_UB_CRDS.P8.csv")
UB_CRDS_long$DATE.TIME <- as.POSIXct(UB_CRDS_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

#load animal and temperature data
animal_temp <- read.csv("20250408-15_LVAT_Animal_Temp_data.csv")
animal_temp$DATE.TIME <- as.POSIXct(animal_temp$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")


# Join animal & temp data with gas datasets
emission_LUFA_FTIR  <- left_join(animal_temp, LUFA_FTIR_long,  by = "DATE.TIME")
emission_ANECO_FTIR <- left_join(animal_temp, ANECO_FTIR_long, by = "DATE.TIME")
emission_MBBM_FTIR  <- left_join(animal_temp, MBBM_FTIR_long,  by = "DATE.TIME")
emission_ATB_FTIR   <- left_join(animal_temp, ATB_FTIR_long,   by = "DATE.TIME")
emission_ATB_CRDS   <- left_join(animal_temp, ATB_CRDS_long,   by = "DATE.TIME")
emission_LUFA_CRDS  <- left_join(animal_temp, LUFA_CRDS_long,  by = "DATE.TIME")
emission_UB_CRDS    <- left_join(animal_temp, UB_CRDS_long,    by = "DATE.TIME")


########### Computation of ventilation rates and emissions #################
# Calculate emissions using the function
result_emission_LUFA_FTIR  <- indirect.CO2.balance.method(emission_LUFA_FTIR)
result_emission_ANECO_FTIR <- indirect.CO2.balance.method(emission_ANECO_FTIR)
result_emission_MBBM_FTIR  <- indirect.CO2.balance.method(emission_MBBM_FTIR)
result_emission_ATB_FTIR   <- indirect.CO2.balance.method(emission_ATB_FTIR)
result_emission_ATB_CRDS   <- indirect.CO2.balance.method(emission_ATB_CRDS)
result_emission_LUFA_CRDS  <- indirect.CO2.balance.method(emission_LUFA_CRDS)
result_emission_UB_CRDS    <- indirect.CO2.balance.method(emission_UB_CRDS)

# Write emissions results to CSV
write.csv(result_emission_LUFA_FTIR,  "result_emission_LUFA_FTIR.csv",  row.names = FALSE)
write.csv(result_emission_ANECO_FTIR, "result_emission_ANECO_FTIR.csv", row.names = FALSE)
write.csv(result_emission_MBBM_FTIR,  "result_emission_MBBM_FTIR.csv",  row.names = FALSE)
write.csv(result_emission_ATB_FTIR,   "result_emission_ATB_FTIR.csv",   row.names = FALSE)
write.csv(result_emission_ATB_CRDS,   "result_emission_ATB_CRDS.csv",   row.names = FALSE)
write.csv(result_emission_LUFA_CRDS,  "result_emission_LUFA_CRDS.csv",  row.names = FALSE)
write.csv(result_emission_UB_CRDS,    "result_emission_UB_CRDS.csv",    row.names = FALSE)

# Extract NH3 and CH4 emission columns from each result
result_emission_LUFA_FTIR <- result_emission_LUFA_FTIR %>%
        select(DATE.TIME, matches("^emission_NH3"), matches("^emission_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

result_emission_ANECO_FTIR <- result_emission_ANECO_FTIR %>%
        select(DATE.TIME, matches("^emission_NH3"), matches("^emission_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

result_emission_MBBM_FTIR <- result_emission_MBBM_FTIR %>%
        select(DATE.TIME, matches("^emission_NH3"), matches("^emission_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

result_emission_ATB_FTIR <- result_emission_ATB_FTIR %>%
        select(DATE.TIME, matches("^emission_NH3"), matches("^emission_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

result_emission_ATB_CRDS <- result_emission_ATB_CRDS %>%
        select(DATE.TIME, matches("^emission_NH3"), matches("^emission_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

result_emission_LUFA_CRDS <- result_emission_LUFA_CRDS %>%
        select(DATE.TIME, matches("^emission_NH3"), matches("^emission_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

result_emission_UB_CRDS <- result_emission_UB_CRDS %>%
        select(DATE.TIME, matches("^emission_NH3"), matches("^emission_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)


# Repeat for all dataframes
final_emission_combined <- full_join(result_emission_LUFA_FTIR, result_emission_ANECO_FTIR, by = "DATE.TIME") %>%
        full_join(result_emission_MBBM_FTIR, by = "DATE.TIME") %>%
        full_join(result_emission_ATB_FTIR, by = "DATE.TIME") %>%
        full_join(result_emission_ATB_CRDS, by = "DATE.TIME") %>%
        full_join(result_emission_LUFA_CRDS, by = "DATE.TIME") %>%
        full_join(result_emission_UB_CRDS, by = "DATE.TIME")

write.csv(final_emission_combined,    "20250408-15_final_ringversuch_emission_results_combined.csv",    row.names = FALSE)


########## Statistics and Data Visualization of emissions ############
final_emission_combined <- final_emission_combined %>%
        filter(DATE.TIME <= as.POSIXct("2025-04-15 00:00:00"))

colSums(!is.na(final_emission_combined)) #Total measurement time period must be 6.5 days or 157 hours

# Summary statistics
sumtable(final_emission_combined, 
         summ = c("mean", "sd", "median", "min", "max"),
         group = FALSE, # no grouping; use TRUE to split by a factor
         digits = 2)
