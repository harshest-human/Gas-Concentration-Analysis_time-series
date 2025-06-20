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
ANECO_FTIR <- read.csv("D:/Data Analysis/Gas-Concentration-Analysis_time-series/Ringversuche_clean/20250408-15_hourly_ANECO_FTIR.csv")
MBBM_FTIR <- read.csv("D:/Data Analysis/Gas-Concentration-Analysis_time-series/Ringversuche_clean/20250408-15_hourly_MBBM_FTIR.csv")
ATB_FTIR <- read.csv("D:/Data Analysis/Gas-Concentration-Analysis_time-series/Ringversuche_clean/20250408-15_hourly_ATB_FTIR1.csv")
ATB_CRDS <- read.csv("D:/Data Analysis/Gas-Concentration-Analysis_time-series/Ringversuche_clean/20250408-15_hourly_ATB_CRDS.P8.csv")
LUFA_CRDS <- read.csv("D:/Data Analysis/Gas-Concentration-Analysis_time-series/Ringversuche_clean/20250408-15_hourly_LUFA_CRDS.P8.csv")
UB_CRDS <- read.csv("D:/Data Analysis/Gas-Concentration-Analysis_time-series/Ringversuche_clean/20250408-15_hourly_UB_CRDS.P8.csv")