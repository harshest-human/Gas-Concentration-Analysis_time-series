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
FTIR.comb <- fread("D:/Data Analysis/GasmetCX4000_FTIR_Gas_Measurement/FTIR.comb.csv")
CRDS.comb <- fread("D:/Data Analysis/Picarro-G2508_CRDS_gas_measurement/CRDS.comb.csv")


######## Data combining ##########
# Format Date and time
FTIR.comb$DATE.TIME = as.POSIXct(FTIR.comb$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")
CRDS.comb$DATE.TIME = as.POSIXct(CRDS.comb$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

# convert into data.table
data.table::setDT(FTIR.comb)
data.table::setDT(CRDS.comb)

####### Read and write ##########
GAS.comb <- FTIR.comb[CRDS.comb, on = .(DATE.TIME), roll = "nearest"]
write.csv(FTIR.comb, "GAS.comb.csv", row.names = FALSE)
GAS.comb$DATE.TIME = as.POSIXct(GAS.comb$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")


####### Data Analysis ########
# Plotting using ggplot2 
ggplot(GAS.comb, aes(x = factor(Messstelle.F1), y = CO2.F1)) +
        geom_line(aes(color = "Messstelle.F1"), alpha = 0.5) +
        geom_line(aes(x = factor(Messstelle.F2), y = CO2.F2, color = "Messstelle.F2"), alpha = 0.5) +
        theme_minimal()


# Plotting using ggline
GAS.comb$Messstelle.F1 <- as.factor(GAS.comb$Messstelle.F1)
GAS.comb$Messstelle.F2 <- as.factor(GAS.comb$Messstelle.F2)


# Plotting using ggplot2 
ggplot(GAS.comb, aes(x = factor(MPVPosition.P8), y = CO2.P8)) +
        geom_boxplot(aes(color = "MPVPosition.P8"), alpha = 0.5) +
        geom_boxplot(aes(x = factor(MPVPosition.P9), y = CO2.P9, color = "MPVPosition.P9"), alpha = 0.5) +
        labs(x = "MPVPosition", y = "CO2 Mean") +
        scale_color_manual(values = c("MPVPosition.P8" = "blue", "MPVPosition.P9" = "red")) +
        theme_minimal()

# Plotting using ggplot2 
ggplot(GAS.comb, aes(x = factor(Messstelle.F1), y = CO2.F1)) +
        geom_boxplot(aes(color = "Messstelle.F1"), alpha = 0.5) +
        geom_boxplot(aes(x = factor(Messstelle.F2), y = CO2.F2, color = "Messstelle.F2"), alpha = 0.5) +
        scale_color_manual(values = c("Messstelle.F1" = "blue", "Messstelle.F2" = "red")) +
        theme_minimal()


# Plotting using ggline
ggline(GAS.comb, x = "MPVPosition.P8", y = "CO2.P8",
       add = "mean_se",
       linetype = "solid") 

# Plotting using ggline
ggline(GAS.comb, x = "MPVPosition.P9", y = "CO2.P9",
       add = "mean_se",
       linetype = "solid")

# Plotting using ggline
ggline(GAS.comb, x = "Messstelle.F2", y = "NH3.F2",
       add = "mean_se",
       linetype = "solid",
       legend = "right")

# Plotting using ggline
ggline(GAS.comb, x = "Messstelle.F1", y = "NH3.F1",
       add = "mean_se",
       linetype = "solid",
       legend = "right")

