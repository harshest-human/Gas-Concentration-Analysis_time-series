getwd()
######## Load library #########
library(tidyverse)
library(reshape2)
library(lubridate)
library(psych)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(data.table)

######## Import Gas Data #########
#load gas concentration data
lab_combined <- read.csv("20250408-15_Ringversuche_combined_data.csv")

#create a new dataframe 
emission_data <- data.frame(
        DATE.TIME = as.POSIXct(character()),   # datetime column
        hour = integer(),                       # hour of the day
        n_animals = integer(),                 # number of animals
        m_weight_kg = numeric(),               # mean weight in kg
        p_pregnancy_day = numeric(),           # percentage in pregnancy days
        Y1_milk_prod_kg_day = numeric(),       # milk production (kg/day)
        CO2_out = numeric(),                   # CO2 concentration outside
        CO2_in = numeric(),                    # CO2 concentration inside
        NH3_out = numeric(),                   # NH3 concentration outside
        NH3_in = numeric(),                    # NH3 concentration inside
        CH4_out = numeric(),                   # CH4 concentration outside
        CH4_in = numeric(),                    # CH4 concentration inside
        Temperature = numeric()                # temperature in °C
        )

# Add constants
emission_data <- emission_data %>%
        mutate(
                P = 0.185,                      # CO2 per hpu (m³/h/hpu), from Pedersen 2008
                a = 0.22,                       # amplitude coefficient
                h_min = 2.9,                    # time of minimum activity
                CO2_Molmass = 44.01,            # CO2 molar mass (g/mol)
                NH3_Molmass = 17.031,           # NH3 molar mass (g/mol)
                CH4_Molmass = 16.04,            # CH4 molar mass (g/mol)
                R = 8.314472,                   # gas constant (J/mol·K)
                p = 1013                        # reference pressure (mbar)
        )


