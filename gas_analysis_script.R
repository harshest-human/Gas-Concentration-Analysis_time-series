getwd()
######## Load library #########
library(tidyverse)
library(reshape2)
library(hablar)
library(lubridate)
library(psych)
library(ggplot2)
library(readxl)
library(dplyr)
library(ggpubr)
library(readr)
library(data.table)



######## Import Gas Data #########
FTIR.comb <- read.csv("D:/Data Analysis/GasmetCX4000_FTIR_Gas_Measurement/FTIR.comb.csv")
CRDS.comb <- CRDS.comb <- read.csv("D:/Data Analysis/Picarro-G2508_CRDS_gas_measurement/CRDS.comb.csv")


######## Data combining ##########
data.table::setDT(FTIR.comb)
data.table::setDT(CRDS.comb)

FTIR.comb$DATE.TIME = as.POSIXct(FTIR.comb$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")
CRDS.comb$DATE.TIME = as.POSIXct(CRDS.comb$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

####### Read and write ##########
GAS.comb <- FTIR.comb[CRDS.comb, on = .(DATE.TIME), roll = "nearest"]
write.csv(FTIR.comb, "GAS.comb.csv", row.names = FALSE)
GAS.comb <- read.csv("GAS.comb.csv")


####### Data Analysis ########