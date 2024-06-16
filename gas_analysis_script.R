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
FTIR.comb <- fread("D:/Data Analysis/GasmetCX4000_FTIR_Gas_Measurement/FTIR.comb.csv")
CRDS.comb <- fread("D:/Data Analysis/Picarro-G2508_CRDS_gas_measurement/CRDS.comb.csv")


######## Data combining ##########
# Format Date and time
FTIR.comb$DATE.TIME = as.POSIXct(FTIR.comb$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")
CRDS.comb$DATE.TIME = as.POSIXct(CRDS.comb$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

# convert into data.table
data.table::setDT(FTIR.comb)
data.table::setDT(CRDS.comb)

# combine FTIR and CRDS
GAS.comb <- FTIR.comb[CRDS.comb, on = .(DATE.TIME), roll = "nearest"]

# write the combined dataframe
write.csv(GAS.comb, "GAS.comb.csv", row.names = FALSE)



######## Data reshaping ##########
# Import the final combined dataframe
GAS.comb <- read.csv("D:/Data Analysis/Gas-Concentration-Analysis_time-series/GAS.comb.csv")
GAS.comb$DATE.TIME = as.POSIXct(GAS.comb$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

# Convert GAS.comb to a data.table if it's not already
setDT(GAS.comb)

# Subset the columns of interest using data.table syntax
GAS.comb <- GAS.comb[, .(DATE.TIME,
                         MPVPosition.P9,
                         Messstelle.F2,
                         MPVPosition.P8,
                         Messstelle.F1,
                         CO2.P9, NH3.P9, CH4.P9, H2O.P9,
                         CO2.F2, NH3.F2, CH4.F2, H2O.F2,
                         CO2.P8, NH3.P8, CH4.P8, H2O.P8,
                         CO2.F1, NH3.F1, CH4.F1, H2O.F1)]

# Convert columns to factors with specified levels and labels
GAS.comb$MPVPosition.P9 <- factor(GAS.comb$MPVPosition.P9, levels = 1:16, labels = 1:16)
GAS.comb$Messstelle.F2 <- factor(GAS.comb$Messstelle.F2, levels = 1:10, labels = 17:26)
GAS.comb$MPVPosition.P8 <- factor(GAS.comb$MPVPosition.P8, levels = 1:16, labels = 27:42)
GAS.comb$Messstelle.F1 <- factor(GAS.comb$Messstelle.F1, levels = 1:10, labels = 43:52)


# Create a new dataframe with the desired structure
GAS.long <- data.table(
        DATE.TIME = GAS.comb$DATE.TIME,
        ID = rep(c("MPVPosition.P9", "Messstelle.F2", "MPVPosition.P8", "Messstelle.F1"), each = nrow(GAS.comb)),
        sampling.point = c(GAS.comb$MPVPosition.P9, rep(GAS.comb$Messstelle.F2, 1), GAS.comb$MPVPosition.P8, rep(GAS.comb$Messstelle.F1, 1)),
        CO2 = c(GAS.comb$CO2.P9, GAS.comb$CO2.F2, GAS.comb$CO2.P8, GAS.comb$CO2.F1),
        CH4 = c(GAS.comb$CH4.P9, GAS.comb$CH4.F2, GAS.comb$CH4.P8, GAS.comb$CH4.F1),
        NH3 = c(GAS.comb$NH3.P9, GAS.comb$NH3.F2, GAS.comb$NH3.P8, GAS.comb$NH3.F1),
        H2O = c(GAS.comb$H2O.P9, GAS.comb$H2O.F2, GAS.comb$H2O.P8, GAS.comb$H2O.F1)
)

# Convert 'GAS.long' to data.table
setDT(GAS.long)

# write after arranging columns
GAS.long <- GAS.long[, .(DATE.TIME, ID, sampling.point, CO2, CH4, NH3, H2O)]

write.csv(GAS.long, "GAS.long.csv", row.names = FALSE)


####### Data Analysis ########
GAS.long <- fread("GAS.long.csv")

GAS.long <- GAS.long[sampling.point != 52]
GAS.long$sampling.point <- as.factor(GAS.long$sampling.point)

# Calculate average values for each gas
avg_CO2 <- GAS.long[, mean(CO2, na.rm = TRUE)]
avg_CH4 <- GAS.long[, mean(CH4, na.rm = TRUE)]
avg_NH3 <- GAS.long[, mean(NH3, na.rm = TRUE)]
avg_H2O <- GAS.long[, mean(H2O, na.rm = TRUE)]

# Calculate relative errors for each gas
GAS.long[, Err_CO2 := ((CO2 - avg_CO2) / avg_CO2) * 100, by = sampling.point]
GAS.long[, Err_CH4 := ((CH4 - avg_CH4) / avg_CH4) * 100, by = sampling.point]
GAS.long[, Err_NH3 := ((NH3 - avg_NH3) / avg_NH3) * 100, by = sampling.point]
GAS.long[, Err_H2O := ((H2O - avg_H2O) / avg_H2O) * 100, by = sampling.point]


# Ensure GAS.long is a data.table
setDT(GAS.long)

# Define the vertical_groups list
vertical_groups <- list(
        top = c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49),
        mid = c(2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50),
        bottom = c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51)
)

# Assign groups based on sampling point
GAS.long[, vertical := ifelse(sampling.point %in% vertical_groups$top, "top",
                              ifelse(sampling.point %in% vertical_groups$mid, "mid",
                                     "bottom"))]

# Define colors for vertical groups
point_fill <- c("top" = "orange", "mid" = "green3", "bottom" = "steelblue1")


# Plot CO2 standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = Err_CO2, fill = vertical)) +
        geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
        geom_point(stat = "summary", fun = "mean", size = 3, shape = 21) +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        labs(x = "Sampling Point", y = "Relative Error (%)", title = "Relative Error of CO2 Concentration by Sampling Point") +
        scale_fill_manual(values = point_fill) +
        scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 20)) +
        theme_minimal() + guides(fill = FALSE)  


# Plot CH4 standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = Err_CH4, fill = vertical)) +
        geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
        geom_point(stat = "summary", fun = "mean", size = 3, shape = 21) +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        labs(x = "Sampling Point", y = "Relative Error (%)", title = "Relative Error of CH4 Concentration by Sampling Point") +
        scale_fill_manual(values = point_fill) +
        scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 20)) +
        theme_minimal() + guides(fill = FALSE)  

# Plot NH3 standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = Err_NH3, fill = vertical)) +
        geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
        geom_point(stat = "summary", fun = "mean", size = 3, shape = 21) +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        labs(x = "Sampling Point", y = "Relative Error (%)", title = "Relative Error of NH3 Concentration by Sampling Point") +
        scale_fill_manual(values = point_fill) +
        scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 20)) +
        theme_minimal() + guides(fill = FALSE)  

# Plot H2O standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = Err_H2O, fill = vertical)) +
        geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
        geom_point(stat = "summary", fun = "mean", size = 3, shape = 21) +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        labs(x = "Sampling Point", y = "Relative Error (%)", title = "Relative Error of H2O % by Sampling Point") +
        scale_fill_manual(values = point_fill) +
        scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 20)) +
        theme_minimal() + guides(fill = FALSE)  
