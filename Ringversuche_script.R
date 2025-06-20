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
ANECO_FTIR <- read.csv("20250408-15_hourly_ANECO_FTIR.csv")
MBBM_FTIR <- read.csv("20250408-15_hourly_MBBM_FTIR.csv")
ATB_FTIR <- read.csv("20250408-15_hourly_ATB_FTIR.1.csv")
ATB_CRDS <- read.csv("20250408-15_hourly_ATB_CRDS.P8.csv")
LUFA_CRDS <- read.csv("20250408-15_hourly_LUFA_CRDS.P8.csv")
UB_CRDS <- read.csv("20250408-15_hourly_UB_CRDS.P8.csv")
lab_combined <- bind_rows(ANECO_FTIR, MBBM_FTIR, ATB_FTIR, ATB_CRDS, LUFA_CRDS, UB_CRDS)

####### Data Visualization ##########
lab_combined.1 <- lab_combined %>% filter(DATE.TIME >= "2025-04-08 12:00:00" & DATE.TIME <= "2025-04-09 12:00:00")
lab_combined.1 <- lab_combined %>% filter(analyzer != "FTIR.1")

ggplot(lab_combined.1, aes(x = DATE.TIME, y = CO2, color = lab)) +
        geom_point() +
        facet_wrap(~ location) +
        labs(title = "CO2 (ppmv) by Lab over Time", x = "Hour", y = "CO2 (ppmv)") +
        theme_minimal()



