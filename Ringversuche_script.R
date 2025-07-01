getwd()
######## Load library #########
library(tidyverse)
library(reshape2)
library(lubridate)
library(psych)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(scales)
library(gridExtra)
library(magick)
library(summarytools)
source("function_indirect.CO2.balance.method.R")

######## Import Gas Data #########
# Load animal and temperature data
animal_temp <- read.csv("20250408-15_LVAT_Animal_Temp_data.csv")

# Read and convert DATE.TIME for FTIR data
ATB_FTIR <- read.csv("20250408-15_long_ATB_FTIR.1.csv")
LUFA_FTIR <- read.csv("20250408-15_long_LUFA_FTIR.2.csv")
ANECO_FTIR <- read.csv("20250408-15_long_ANECO_FTIR.3.csv")
MBBM_FTIR <- read.csv("20250408-15_long_MBBM_FTIR.4.csv")

# Read and convert DATE.TIME for CRDS data
ATB_CRDS <- read.csv("20250408-15_long_ATB_CRDS.1.csv")
UB_CRDS <- read.csv("20250408-15_long_UB_CRDS.2.csv")
LUFA_CRDS <- read.csv("20250408-15_long_LUFA_CRDS.3.csv")

# Combine all data set
gas_data <- bind_rows(LUFA_FTIR, ANECO_FTIR, MBBM_FTIR, ATB_FTIR, ATB_CRDS, LUFA_CRDS, UB_CRDS)
input_combined <-full_join(gas_data, animal_temp, by = "DATE.TIME")

# Write csv
input_combined <- input_combined %>% select(DATE.TIME, hour, everything())
write.csv(input_combined, "20250408-15_ringversuche_input_combined_data.csv", row.names = FALSE)

######## Computation of ventilation rates and emissions #########
# Convert DATE.TIME format
input_combined <- input_combined %>% filter(DATE.TIME >= "2025-04-08 12:00:00" & DATE.TIME <= "2025-04-14 23:00:00")

# Calculate emissions using the function
emission_combined  <- indirect.CO2.balance(input_combined)

# Write csv
write.csv(emission_combined, "20250408-15_ringversuche_emission_combined_data.csv", row.names = FALSE)


########## Data visualization ############
#Development of plot functions
emicon.plot <- function(df, x, y) {
        library(dplyr)
        library(ggpubr)
        library(ggplot2)
        library(scales)
        
        # Define y-axis labels with "computed by" for ventilation and background info for others
        ylabels <- list(
                # Concentrations ppm
                "delta_NH3_N_ppm" = "NH3 Concentration Δ (ppm) (computed by North background)",
                "delta_CH4_N_ppm" = "CH4 Concentration Δ (ppm) (computed by North background)",
                "delta_CO2_N_ppm" = "CO2 Concentration Δ (ppm) (computed by North background)",
                "delta_NH3_S_ppm" = "NH3 Concentration Δ (ppm) (computed by South background)",
                "delta_CH4_S_ppm" = "CH4 Concentration Δ (ppm) (computed by South background)",
                "delta_CO2_S_ppm" = "CO2 Concentration Δ (ppm) (computed by South background)",
                
                # Concentrations mg/m³
                "delta_NH3_N_mgm3" = "NH3 Concentration Δ (mg/m³) (computed by North background)",
                "delta_CH4_N_mgm3" = "CH4 Concentration Δ (mg/m³) (computed by North background)",
                "delta_CO2_N_mgm3" = "CO2 Concentration Δ (mg/m³) (computed by North background)",
                "delta_NH3_S_mgm3" = "NH3 Concentration Δ (mg/m³) (computed by South background)",
                "delta_CH4_S_mgm3" = "CH4 Concentration Δ (mg/m³) (computed by South background)",
                "delta_CO2_S_mgm3" = "CO2 Concentration Δ (mg/m³) (computed by South background)",
                
                # Ventilation rates
                "Q_Vent_rate_N" = "Ventilation Rate (m³/h) (computed by North background)",
                "Q_Vent_rate_S" = "Ventilation Rate (m³/h) (computed by South background)",
                
                # Emissions g/h
                "e_NH3_N" = "NH3 Emission (g/h) (North background)",
                "e_CH4_N" = "CH4 Emission (g/h) (North background)",
                "e_CO2_N" = "CO2 Emission (g/h) (North background)",
                "e_NH3_S" = "NH3 Emission (g/h) (South background)",
                "e_CH4_S" = "CH4 Emission (g/h) (South background)",
                "e_CO2_S" = "CO2 Emission (g/h) (South background)"
        )
        
        # Check if x exists in df
        if (!(x %in% names(df))) stop(paste("x variable", x, "not found in data frame"))
        if (!(y %in% names(df))) stop(paste("y variable", y, "not found in data frame"))
        
        # Set plot title y-label from list or default to y
        y_label <- ifelse(!is.null(ylabels[[y]]), ylabels[[y]], y)
        
        # Basic plot with ggline
        p <- ggline(df, x = x, y = y,
                    add = "mean_se",
                    color = "analyzer") +
                labs(title = paste(y_label, "Trends (mean ± SE)"),
                     x = x,
                     y = y_label,
                     color = "Analyzer") +
                scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
                theme_light() +
                theme(legend.position = "top")
        
        # If x is datetime, add datetime scale and rotate labels
        if (inherits(df[[x]], "POSIXct") || inherits(df[[x]], "POSIXt")) {
                p <- p +
                        scale_x_datetime(date_breaks = "6 hours", date_labels = "%d.%m %H:%M") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
        
        # If x = "hour" (numeric from 0 to 23), set breaks for every hour
        if (x == "hour") {
                p <- p + 
                        scale_x_continuous(breaks = 0:23) +
                        theme(axis.text.x = element_text(angle = 0))
        }
        
        return(p)
}


# DATE.TIME conversion
result <- emission_combined %>% mutate(DATE.TIME <- as.POSIXct(emission_combined$DATE.TIME, format = "%Y-%m-%d %H:%M:%S"))

# Concentration plots in ppm
emicon.plot(df = result, x = "hour", y = "delta_NH3_N_ppm")
emicon.plot(df = result, x = "hour", y = "delta_NH3_S_ppm")
emicon.plot(df = result, x = "hour", y = "delta_CH4_N_ppm")
emicon.plot(df = result, x = "hour", y = "delta_CH4_S_ppm")
emicon.plot(df = result, x = "hour", y = "delta_CO2_N_ppm")
emicon.plot(df = result, x = "hour", y = "delta_CO2_S_ppm")

# Concentration plots in mgm3
emicon.plot(df = result, x = "hour", y = "delta_NH3_N_mgm3")
emicon.plot(df = result, x = "hour", y = "delta_NH3_S_mgm3")
emicon.plot(df = result, x = "hour", y = "delta_CH4_N_mgm3")
emicon.plot(df = result, x = "hour", y = "delta_CH4_S_mgm3")
emicon.plot(df = result, x = "hour", y = "delta_CO2_N_mgm3")
emicon.plot(df = result, x = "hour", y = "delta_CO2_S_mgm3")


# Ventilation rate plots
emicon.plot(df = result, x = "hour", y = "Q_Vent_rate_N")
emicon.plot(df = result, x = "hour", y = "Q_Vent_rate_N")


#emission plots
emicon.plot(df = result, x = "hour", y = "e_CH4_N")
emicon.plot(df = result, x = "hour", y = "e_CH4_S")
emicon.plot(df = result, x = "hour", y = "e_NH3_N")
emicon.plot(df = result, x = "hour", y = "e_NH3_S")


# save plots
y_vars <- c("Q_Vent_rate_N", "Q_Vent_rate_S",
            "delta_NH3_N_ppm", "delta_NH3_S_ppm",
            "delta_CH4_N_ppm", "delta_CH4_S_ppm",
            "delta_CO2_N_ppm", "delta_CO2_S_ppm",
            "delta_NH3_N_mgm3", "delta_NH3_S_mgm3",
            "delta_CH4_N_mgm3", "delta_CH4_S_mgm3",
            "delta_CO2_N_mgm3", "delta_CO2_S_mgm3",
            "e_CH4_N", "e_CH4_S")

for (y_var in y_vars) {
        p <- emicon.plot(df = result, x = "hour", y = y_var)
        
        file_name <- paste0(y_var, ".png")
        
        ggsave(filename = file_name,
               plot = p,
               width = 8,
               height = 5,
               dpi = 300)
        
        cat("Plot saved as", file_name, "\n")
}

