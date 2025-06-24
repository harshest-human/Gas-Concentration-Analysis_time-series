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
lab_combined$lab.anaylzer <- paste(lab_combined$lab, lab_combined$analyzer, sep = "_")
lab_combined$DATE.TIME <- as.POSIXct(lab_combined$DATE.TIME)
lab_combined <- lab_combined %>% select(DATE.TIME, location, lab.anaylzer, CO2, CH4, NH3)
lab_combined <- lab_combined %>% filter(DATE.TIME >= "2025-04-08 12:00:00" & DATE.TIME <= "2025-04-15 12:00:00")
lab_combined <- lab_combined %>% mutate(hour = hour(DATE.TIME))
lab_combined$hour <- as.factor(lab_combined$hour)

#write csv
lab_combined <- lab_combined %>% select(DATE.TIME, hour, everything())
write.csv(lab_combined, "20250408-15_Ringversuche_lab_combined_data.csv", row.names = FALSE)

####### Data Visualization Weekly ##########
#CO2
ggline(lab_combined,
       x = "DATE.TIME",
       y = "CO2",
       add = "mean_se",
       color = "lab.anaylzer",
       facet.by = "location",
       xlab = "Date",
       ylab = "CO2 (ppmv)",
       title = "CO2 (ppmv) over Time (Mean ± SE)") +
        scale_x_datetime(date_breaks = "24 hours", date_labels = "%d.%m", expand = c(0.1, 0.1))  +
        scale_y_continuous(breaks = seq(0, 10000, by = 10)) +
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size =8))


#CH4
ggline(lab_combined,
       x = "DATE.TIME",
       y = "CH4",
       add = "mean_se",
       color = "lab.anaylzer",
       facet.by = "location",
       xlab = "Date",
       ylab = "CH4 (ppmv)",
       title = "CH4 (ppmv) over Time (Mean ± SE)") +
        scale_x_datetime(date_breaks = "24 hours", date_labels = "%d.%m", expand = c(0.1, 0.1))  +
        scale_y_continuous(breaks = seq(0, 1000, by = 1)) +
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size =8))

#NH3
ggline(lab_combined,
       x = "DATE.TIME",
       y = "NH3",
       add = "mean_se",
       color = "lab.anaylzer",
       facet.by = "location",
       xlab = "Date",
       ylab = "NH3 (ppmv)",
       title = "NH3 (ppmv) over Time (Mean ± SE)") +
        scale_x_datetime(date_breaks = "24 hours", date_labels = "%d.%m", expand = c(0.1, 0.1))  +
        scale_y_continuous(breaks = seq(0, 100, by = 0.1)) +
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size =8))


####### Data Visualization Hourly ##########
#CO2
ggline(lab_combined,
       x = "hour",
       y = "CO2",
       add = "mean_se",
       color = "lab.anaylzer",
       facet.by = "location",
       xlab = "Hour",
       ylab = "CO2 (ppmv)",
       title = "CO2 (ppmv) over Time (Mean ± SE)") +
        scale_y_continuous(breaks = seq(0, 10000, by = 10)) +
        theme_light()


#CH4
ggline(lab_combined,
       x = "hour",
       y = "CH4",
       add = "mean_se",
       color = "lab.anaylzer",
       facet.by = "location",
       xlab = "Hour",
       ylab = "CH4 (ppmv)",
       title = "CH4 (ppmv) over Time (Mean ± SE)") +
        scale_y_continuous(breaks = seq(0, 10000, by = 10)) +
        theme_light()


#NH3
ggline(lab_combined,
       x = "hour",
       y = "NH3",
       add = "mean_se",
       color = "lab.anaylzer",
       facet.by = "location",
       xlab = "Hour",
       ylab = "NH3 (ppmv)",
       title = "NH3 (ppmv) over Time (Mean ± SE)") +
        scale_y_continuous(breaks = seq(0, 10000, by = 0.1)) +
        theme_light()




