# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

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

# Create the ggplot
ggplot(CRDS.long, aes(x = DATE.TIME, y = CO2, color = sampling.point)) +
        geom_point() + 
        theme_minimal()

# Plot CO2 standard error bar
ggplot(GAS.long, aes(x = DATE.TIME, y = CO2, fill = sampling.point)) +
        geom_line(stat = "summary", fun = mean, aes(group = 1)) +
        geom_point(stat = "summary", fun = mean, size = 3, shape = 21) +
        geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2) +
        labs(x = "Sampling Point", y = "CO2") +
        theme_minimal() +
        guides(fill = FALSE) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red")




# Load ggpubr if not already loaded
library(ggpubr)

# Ensure sampling.point is treated as a factor for discrete coloring and shaping
CRDS.long$sampling.point <- as.factor(CRDS.long$sampling.point)

# Plot CO2 standard error bar using ggline
ggline(CRDS.long, x = "DATE.TIME", y = "CO2", 
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


