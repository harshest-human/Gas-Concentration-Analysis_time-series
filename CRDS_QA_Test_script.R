# Load necessary libraries
########### Load packages ############
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggpubr)
library(lubridate)

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

# Filter rows by 'DATE.TIME' range before merging
CRDS.long$DATE.TIME <- as.POSIXct(CRDS.long$DATE.TIME, format="%Y-%m-%d %H:%M:%S")  
CRDS.long <- CRDS.long[DATE.TIME >= "2024-11-06 14:00:00" & DATE.TIME <= "2024-11-11 10:00:00"]


########### Data processing (hourly averages) ############
# Extract the hour from the DATE.TIME column
CRDS.long$Hour <- floor_date(CRDS.long$DATE.TIME, unit = "4 hours")

# Calculate the hourly average for each ID, sampling.point, and hour
CRDS.long <- CRDS.long %>% group_by(ID, sampling.point, Hour) %>%
        summarise(CO2 = mean(CO2, na.rm = TRUE),
                CH4 = mean(CH4, na.rm = TRUE),
                NH3 = mean(NH3, na.rm = TRUE),
                H2O = mean(H2O, na.rm = TRUE))


########### Data Visualization ggplot2 ############
CRDS.long$sampling.point = as.factor(CRDS.long$sampling.point)
CRDS.long$ID = as.factor(CRDS.long$ID)

ggplot(CRDS.long, aes(x = Hour, y = CO2, color = ID)) + 
        geom_line() +  
        scale_x_datetime(date_breaks = "4 hour", date_labels = "%Y-%m-%d %H:%M") +  
        scale_color_manual(values = c("red", "blue")) +
        theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))

# Plotting CH4 trends with DATE.TIME on x-axis, faceted by sampling.point and colored by ID
ggplot(CRDS.long, aes(x = Hour, y = CH4, color = ID)) + 
        geom_line() +  
        scale_x_datetime(date_breaks = "2 hour", date_labels = "%Y-%m-%d %H:%M") +  
        scale_color_manual(values = c("red", "blue")) +
        theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))

# Plotting NH3 trends with DATE.TIME on x-axis, faceted by sampling.point and colored by ID
ggplot(CRDS.long, aes(x = Hour, y = NH3, color = ID)) + 
        geom_line() +  
        scale_x_datetime(date_breaks = "2 hour", date_labels = "%Y-%m-%d %H:%M") +  
        scale_color_manual(values = c("red", "blue")) +
        theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))









