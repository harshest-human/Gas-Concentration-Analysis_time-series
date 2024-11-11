# Load necessary libraries
########### Load packages ############
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggpubr)

# Load the data
########### Import combined dataframe of CRDS gas analysers ############
CRDS.comb <- fread("2024_Nov_06_to_11_CRDS.comb.csv")

# Reshape the data for ggplot
CRDS.long <- data.table(
        DATE.TIME = CRDS.comb$DATE.TIME,
        ID = rep(c("MPVPosition.P8", "MPVPosition.P9"), each = nrow(CRDS.comb)),
        sampling.point = c(CRDS.comb$MPVPosition.P9, rep(CRDS.comb$Messstelle.F2, 1), CRDS.comb$MPVPosition.P8, rep(CRDS.comb$Messstelle.F1, 1)),
        CO2 = c(CRDS.comb$CO2.P8, CRDS.comb$CO2.P9),
        CH4 = c(CRDS.comb$CH4.P8, CRDS.comb$CH4.P9),
        NH3 = c(CRDS.comb$NH3.P8, CRDS.comb$NH3.P9),
        H2O = c(CRDS.comb$H2O.P8, CRDS.comb$H2O.P9))


# Convert 'CRDS.long' to data.table
setDT(CRDS.long)

# write after arranging columns
CRDS.long <- CRDS.long[, .(DATE.TIME, ID, sampling.point, CO2, CH4, NH3, H2O)]

write.csv(GAS.long, "2024_Nov_06_to_11_CRDS.long.csv", row.names = FALSE)

# Convert DATE.TIME to datetime format if necessary
CRDS.long$DATE.TIME <- as.POSIXct(CRDS.long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

########### Import reshaped dataframe ############
CRDS.long <- fread("2024_Nov_06_to_11_CRDS.long.csv")

CRDS.long$sampling.point = as.factor(CRDS.long$sampling.point)

# Filter rows by 'DATE.TIME' range before merging
start_date <- as.POSIXct("2024-11-07 13:00:00")  # Set start date
end_date <- as.POSIXct("2024-11-07 15:00:00")  # Set end date


CRDS.long <- CRDS.long[DATE.TIME >= start_date & DATE.TIME <= end_date]


# Create the ggplot
ggplot(CRDS.long, aes(x = DATE.TIME, y = CO2, color = ID, facet(sampling.point))) +
        geom_line() + 
        theme_minimal()


# CalculGeomLine# Calculate average values for each gas
#avg_CO2 <- CRDS.long[, mean(CO2, na.rm = TRUE)]
#avg_CH4 <- CRDS.long[, mean(CH4, na.rm = TRUE)]
#avg_NH3 <- CRDS.long[, mean(NH3, na.rm = TRUE)]
#avg_H2O <- CRDS.long[, mean(H2O, na.rm = TRUE)]


# Calculate relative errors for each gas
#CRDS.long[, Err_CO2 := ((CO2 - avg_CO2) / avg_CO2) * 100, by = sampling.point]
#CRDS.long[, Err_CH4 := ((CH4 - avg_CH4) / avg_CH4) * 100, by = sampling.point]
#CRDS.long[, Err_NH3 := ((NH3 - avg_NH3) / avg_NH3) * 100, by = sampling.point]
#CRDS.long[, Err_H2O := ((H2O - avg_H2O) / avg_H2O) * 100, by = sampling.point]


# Calculate ratio
#CRDS.long[, mean_NH3 := (mean(NH3)), by = sampling.point]
#CRDS.long[, mean_CO2 := (mean(CO2)), by = sampling.point]
#CRDS.long[, ratio_NH3_CO2 := (mean_NH3 / mean_CO2) * 10^3]

########### Data Visualization ggplot2 ############
# Plot CO2 standard error bar
#ggplot(CRDS.long, aes(x = DATE.TIME, y = Err_CO2, fill = sampling.point)) +
        #geom_point(stat = "summary", fun = mean, size = 2, shape = 21) +
        #geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2) +
        #theme_minimal() +
        #guides(fill = FALSE)

########### Data Visualization ggline::ggpubr #############
# Ensure sampling.point is treated as a factor for discrete coloring and shaping
#CRDS.long$sampling.point <- as.factor(CRDS.long$sampling.point)

# Plot CO2 standard error bar using ggline
ggline(CRDS.long, x = "DATE.TIME", y = "Err_CO2", 
       add = "mean_se", 
       color = "ID",    
       palette = "jco",
       point.size = 0.5) +
        theme_minimal() +
        guides(color = FALSE, shape = FALSE)

ggline(CRDS.long, x = "DATE.TIME", y = "CH4", 
       add = "mean_se", 
       color = "ID", 
       facet = "sampling.point",   
       palette = "jco",
       point.size = 0.5) +
        theme_minimal() +
        guides(color = FALSE, shape = FALSE)

ggline(CRDS.long, x = "DATE.TIME", y = "NH3", 
       add = "mean_se", 
       color = "ID", 
       facet = "sampling.point",   
       palette = "jco",
       point.size = 0.5) +
        theme_minimal() +
        guides(color = FALSE, shape = FALSE)


