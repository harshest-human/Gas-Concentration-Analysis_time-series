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

######## Import Gas Data #########
#Load processed datasets
LUFA_FTIR <- read.csv("20250408-15_hourly_LUFA_FTIR.csv")
ANECO_FTIR <- read.csv("20250408-15_hourly_ANECO_FTIR.csv")
MBBM_FTIR <- read.csv("20250408-15_hourly_MBBM_FTIR.csv")
ATB_FTIR <- read.csv("20250408-15_hourly_ATB_FTIR.1.csv")
ATB_CRDS <- read.csv("20250408-15_hourly_ATB_CRDS.P8.csv")
LUFA_CRDS <- read.csv("20250408-15_hourly_LUFA_CRDS.P8.csv")
UB_CRDS <- read.csv("20250408-15_hourly_UB_CRDS.P8.csv")

#combine all data set
lab_combined <- bind_rows(LUFA_FTIR, ANECO_FTIR, MBBM_FTIR, ATB_FTIR, ATB_CRDS, LUFA_CRDS, UB_CRDS)
lab_combined$lab.analyzer <- paste(lab_combined$lab, lab_combined$analyzer, sep = "_")
lab_combined$DATE.TIME <- as.POSIXct(lab_combined$DATE.TIME)
lab_combined <- lab_combined %>% select(DATE.TIME, location, lab.analyzer, CO2, CH4, NH3)
lab_combined <- lab_combined %>% filter(DATE.TIME >= "2025-04-08 12:00:00" & DATE.TIME <= "2025-04-15 12:00:00")
lab_combined <- lab_combined %>% mutate(hour = hour(DATE.TIME))
lab_combined$hour <- as.factor(lab_combined$hour)

#write csv
lab_combined <- lab_combined %>% select(DATE.TIME, hour, everything())
write.csv(lab_combined, "20250408-15_final_ringversuche_concentration_combined_data.csv", row.names = FALSE)

####### Data Visualization Weekly ##########
plot_concentration <- function(df, x, y) {
        #load libraries
        library(dplyr)
        library(ggpubr)
        library(ggplot2)
        library(scales)
        
        # Check if x and y exist in the dataframe
        if (!(x %in% colnames(df))) stop(paste("X-axis column", x, "not found in dataframe."))
        if (!(y %in% colnames(df))) stop(paste("Y-axis column", y, "not found in dataframe."))
        
        # Filter out NA or non-finite values from y
        df <- df %>% filter(!is.na(.data[[y]]), is.finite(.data[[y]]))
        
        # Start building plot
        p <- ggline(df,
                    x = x,
                    y = y,
                    add = "mean_se",
                    color = "lab.analyzer",
                    facet.by = "location") +
                labs(
                        title = paste(y, "Concentration Trends (Mean ± SE)"),
                        x = x,
                        y = paste(y, "(ppmv)"),
                        color = "Laboratory"
                ) +
                theme_light() +
                theme(legend.position = "top")
        
        # If x is datetime, add datetime scale and rotate x-axis labels
        if (inherits(df[[x]], "POSIXct")) {
                p <- p +
                        scale_x_datetime(date_breaks = "6 hours", date_labels = "%d.%m %H:%M") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
        
        return(p)
}


# --- Generate DATE.TIME plots ---
plot_CO2_conc_all <- plot_concentration(lab_combined, x = "DATE.TIME", y = "CO2")
plot_CH4_conc_all <- plot_concentration(lab_combined, x = "DATE.TIME", y = "CH4")
plot_NH3_conc_all <- plot_concentration(lab_combined, x = "DATE.TIME", y = "NH3")

# --- View DATE.TIME plots ---
plot_CO2_conc_all
plot_CH4_conc_all
plot_NH3_conc_all

# --- Save DATE.TIME plots ---
ggsave("plot_CO2_conc_all.pdf", plot_CO2_conc_all, width = 10, height = 6)
ggsave("plot_CH4_conc_all.pdf", plot_CH4_conc_all, width = 10, height = 6)
ggsave("plot_NH3_conc_all.pdf", plot_NH3_conc_all, width = 10, height = 6)

# --- Generate hour-based plots ---
plot_CO2_conc_hour <- plot_concentration(lab_combined, x = "hour", y = "CO2")
plot_CH4_conc_hour <- plot_concentration(lab_combined, x = "hour", y = "CH4")
plot_NH3_conc_hour <- plot_concentration(lab_combined, x = "hour", y = "NH3")

# --- View hour-based plots ---
plot_CO2_conc_hour
plot_CH4_conc_hour
plot_NH3_conc_hour

# --- Save hour-based plots ---
ggsave("plot_CO2_conc_hour.pdf", plot_CO2_conc_hour, width = 10, height = 6)
ggsave("plot_CH4_conc_hour.pdf", plot_CH4_conc_hour, width = 10, height = 6)
ggsave("plot_NH3_conc_hour.pdf", plot_NH3_conc_hour, width = 10, height = 6)
