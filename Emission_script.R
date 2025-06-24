getwd()
######## Load library #########
library(tidyverse)
library(reshape2)
library(lubridate)
library(psych)
library(ggplot2)
library(dplyr)
library(ggpubr)

######## Import Gas Data #########
#load gas concentration data
lab_combined <- read.csv("20250408-15_Ringversuche_lab_combined_data.csv")

#create a new dataframe 
emission_data <- data.frame(
        DATE.TIME = as.POSIXct(character()),   # datetime column
        hour = integer(),                       # hour of the day
        n_animals = integer(),                 # number of animals
        m_weight = numeric(),               # mean weight in kg
        p_pregnancy_day = numeric(),           # percentage in pregnancy days
        Y1_milk_prod = numeric(),       # milk production (kg/day)
        CO2_out = numeric(),                   # CO2 concentration outside
        CO2_in = numeric(),                    # CO2 concentration inside
        NH3_out = numeric(),                   # NH3 concentration outside
        NH3_in = numeric(),                    # NH3 concentration inside
        CH4_out = numeric(),                   # CH4 concentration outside
        CH4_in = numeric(),                    # CH4 concentration inside
        Temperature = numeric()                # temperature in °C
        )

# constants
emission_data <- emission_data %>%
        mutate(
                P_CO2_term = 0.185,                      # CO2 per hpu (m³/h/hpu), from Pedersen 2008
                a = 0.22,                       # amplitude coefficient
                h_min = 2.9,                    # time of minimum activity
                CO2_Molmass = 44.01,            # CO2 molar mass (g/mol)
                NH3_Molmass = 17.031,           # NH3 molar mass (g/mol)
                CH4_Molmass = 16.04,            # CH4 molar mass (g/mol)
                R_gas_constant = 8.314472,                   # gas constant (J/mol·K)
                p_pressure_ref = 1013                        # reference pressure (mbar)
        )


# Calculation of ventilation rate (Q) by indirect method
emission_data <- emission_data %>%
        mutate(
                # Animal activity corrected by time
                A_corr = 1 - a * 3 * sin((2 * pi / 24) * (hour + 6 - h_min)),
                
                # Heat production per animal (W)
                Phi_tot = 5.6 * m_weight^(0.75) + 22 * Y1_milk_prod + 1.6e-5 * p_pregnancy_day^3,
                
                # Heat production corrected for temperature (W/animal)
                Phi_T_corr = Phi_tot * (1 + 4e-5 * (20 - Temperature)^3),
                
                # hpu corrected for Temperature and Activity for all animals
                hpu_T_A_corr_all_animal = ((Phi_T_corr / 1000) * A_corr) * n_animals,
                
                # Number of Livestock Units
                n_LU = (n_animals * m_weight) / 500,
                
                # CO2 production (m³/h) corrected for T and A
                P_CO2_T_A_all_animal = hpu_T_A_corr_all_animal * P_CO2_term,
                
                # Total ventilation rate (m³/h)
                Q_Vent_rate = P_CO2_T_A_all_animal / ((CO2_in - CO2_out) * 1e-6),
                
                # NH3 concentration difference (mg/m³)
                delta_NH3 = (0.1 * NH3_Molmass * p_pressure_ref * (NH3_in - NH3_out)) / ((Temperature + 273.15) * R_gas_constant),
                
                # NH3 emission (g/h per barn)
                emission_NH3 = (delta_NH3 * Q_Vent_rate) / 1000,
                
                # NH3 emission (kg/year per barn)
                emission_NH3_per_year = emission_NH3 * 24 * 365 / 1000,
                
                # CH4 concentration difference (mg/m³)
                delta_CH4 = (0.1 * CH4_Molmass * p_pressure_ref * (CH4_in - CH4_out)) / ((Temperature + 273.15) * R_gas_constant),
                
                # CH4 emission (g/h per barn)
                emission_CH4 = (delta_CH4 * Q_Vent_rate) / 1000,
                
                # CH4 emission (kg/year per barn)
                emission_CH4_per_year = emission_CH4 * 24 * 365 / 1000
        )


