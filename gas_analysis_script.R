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
library(viridis)
library(gplots)
library(RColorBrewer)



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


######## 0. Plotting errors using ggplot2 #########

# Plot CO2 standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = Err_CO2)) +
        geom_line(stat = "summary", fun.y = "mean", aes(group = 1), size = 1) +
        geom_point(stat = "summary", fun.y = "mean", size = 3, shape = 21, fill = "yellow2") +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        labs(x = "Sampling Point", y = "Relative Error (%)", title = "Relative Error of CO2 Concentration by Sampling Point") +
        theme_minimal()

# Plot CH4 standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = Err_CH4)) +
        geom_line(stat = "summary", fun.y = "mean", aes(group = 1), size = 1) +
        geom_point(stat = "summary", fun.y = "mean", size = 3, shape = 21, fill = "green2") +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        labs(x = "Sampling Point", y = "Relative Error (%)", title = "Relative Error of CH4 Concentration by Sampling Point") +
        theme_minimal()

# Plot NH3 standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = Err_NH3)) +
        geom_line(stat = "summary", fun.y = "mean", aes(group = 1), size = 1) +
        geom_point(stat = "summary", fun.y = "mean", size = 3, shape = 21, fill = "orange2") +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        labs(x = "Sampling Point", y = "Relative Error (%)", title = "Relative Error of NH3 Concentration by Sampling Point") +
        theme_minimal()

# Plot H2O standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = Err_H2O)) +
        geom_line(stat = "summary", fun.y = "mean", aes(group = 1), size = 1) +
        geom_point(stat = "summary", fun.y = "mean", size = 3, shape = 21, fill = "blue2") +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        labs(x = "Sampling Point", y = "Relative Error (%)", title = "Relative Error of H2O Concentration by Sampling Point") +
        theme_minimal()

####### 1. Plotting using ggplot2 ########
ggplot(GAS.long, aes(x = sampling.point, y = CO2, group = 1)) +
        geom_line(color = "blue") +
        geom_point(color = "blue", size = 2) +
        labs(
                title = "Variations in CO2 at Various Sampling Points",
                x = "Sampling Point",
                y = "CO2 Concentration",
                caption = "Data source: Your Source" ) +
        theme_minimal() +
        theme(  plot.title = element_text(size = 16, face = "bold"),
                axis.title = element_text(size = 14),
                axis.text.x = element_text(angle = 45, hjust = 1))

####### 2. Plotting using ggplot2 hour of the day########
# Example data setup (replace with your actual GAS.long dataframe)
# You may need to adjust column names or data types based on your actual dataframe structure.

# Example: Creating a sequence of sampling.location based on every 3 sampling.point
GAS.long$sampling.location <- as.factor((as.numeric(as.factor(GAS.long$sampling.point)) - 1) %/% 3 + 1)

# Convert DATE.TIME to POSIXct format if it's not already
GAS.long$DATE.TIME <- as.POSIXct(GAS.long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

# Extract hour of the day from DATE.TIME
GAS.long$hour <- hour(GAS.long$DATE.TIME)

# Aggregate data by sampling.location and hour, calculating mean CO2 concentration
agg_data <- GAS.long %>%
        group_by(sampling.location, hour) %>%
        summarise(mean_CO2 = mean(CO2))

# Plotting the graph
ggplot(agg_data, aes(x = hour, y = mean_CO2, group = sampling.location, color = sampling.location)) +
        geom_line(size = 1) +
        labs(
                title = "Variation of CO2 Concentration Across 17 Sampling Locations Every Hour",
                x = "Hour of the Day",
                y = "Mean CO2 Concentration",
                color = "Sampling Location"
        ) +
        scale_color_discrete(name = "Sampling Location") +
        theme_minimal() +
        theme(
                plot.title = element_text(size = 16, face = "bold"),
                axis.title = element_text(size = 14),
                legend.position = "bottom"
        )


####### 3. Plotting using ggplot2 hour of the day########
# Convert DATE.TIME to POSIXct format for extracting hours
GAS.long$DATE.TIME <- as.POSIXct(GAS.long$DATE.TIME)

# Extract hour of the day and convert to factor
GAS.long$hour <- factor(hour(GAS.long$DATE.TIME), levels = 0:23)  # Assuming hours range from 0 to 23

# Aggregate CO2 concentrations by hour and sampling point
agg_data <- GAS.long %>%
        group_by(sampling.point, hour) %>%
        summarise(mean_CO2 = mean(CO2),
                  se_CO2 = sd(CO2) / sqrt(n()))  # Calculate standard error

# Plotting with mean line and standard error
ggplot(agg_data, aes(x = hour, y = mean_CO2, color = sampling.point, group = sampling.point)) +
        geom_line() +
        geom_errorbar(aes(ymin = mean_CO2 - se_CO2, ymax = mean_CO2 + se_CO2), width = 0.1) +  # Error bars for standard error
        labs(
                title = "Variation of CO2 Concentration by Hour and Sampling Point",
                x = "Hour of the Day",
                y = "Mean CO2 Concentration",
                color = "Sampling Point"
        ) +
        theme_minimal() +
        theme(
                plot.title = element_text(size = 16, face = "bold"),
                axis.title = element_text(size = 14)
        )


####### 4. Plotting using gplots and RColorBrewer ########
library(gplots)
library(RColorBrewer)

# Ensure DATE.TIME and sampling.point are factors and in correct order
GAS.long$DATE.TIME <- as.factor(GAS.long$DATE.TIME)
GAS.long$sampling.point <- factor(GAS.long$sampling.point, levels = unique(GAS.long$sampling.point))

# Create a matrix for heatmap with sampling.point as rows and DATE.TIME as columns
heatmap_data <- acast(GAS.long, sampling.point ~ DATE.TIME, value.var = "CO2", fun.aggregate = mean)

# Set color palette for heatmap
my_palette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))(100)

# Plot heatmap
heatmap.2(as.matrix(heatmap_data),
          Rowv = NULL, Colv = NULL, dendrogram = "none",
          col = my_palette, trace = "none", margins = c(5, 10),
          main = "Heatmap of CO2 Concentration by Sampling Point",
          xlab = "DATE.TIME", ylab = "Sampling Point",
          key.title = NA, keysize = 1.5, cexRow = 0.7, cexCol = 0.7)


####### 5. Plotting using ggplot2 and virdis ########
library(ggplot2)
library(viridis)

# Assuming GAS.long is your dataframe with the required structure

# Categorize sampling points into three groups: top, mid, bottom
# Define the grouping based on your description
vertical <- cut(as.numeric(GAS.long$sampling.point), breaks = c(0, 16, 32, 51), labels = c("top", "mid", "bottom"))

# Add the vertical grouping to GAS.long
GAS.long <- cbind(GAS.long, vertical)

# Convert sampling.point and vertical to factors for correct ordering in heatmap
GAS.long$sampling.point <- factor(GAS.long$sampling.point, levels = unique(GAS.long$sampling.point))
GAS.long$vertical <- factor(GAS.long$vertical, levels = c("top", "mid", "bottom"))

# Create heatmap using ggplot2
ggplot(GAS.long, aes(x = vertical, y = sampling.point, fill = CO2)) +
        geom_tile(color = "white") +
        scale_fill_viridis(name = "CO2", option = "plasma", limits = range(GAS.long$CO2, na.rm = TRUE)) +
        labs(
                title = "Heatmap of CO2 Variation by Sampling Point Groups",
                x = "Vertical Groups",
                y = "Sampling Point",
                fill = "CO2 Concentration"
        ) +
        theme_minimal() +
        theme(
                plot.title = element_text(size = 16, face = "bold"),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12)
        )


####### 6. Plotting using ggplot2 ########
# Convert sampling.point to factor to ensure correct ordering in plots
GAS.long$sampling.point <- factor(GAS.long$sampling.point, levels = 1:51)

# Define groups for vertical categorization
vertical_groups <- list(
        top = c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49),
        mid = c(2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50),
        bottom = c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51)
)

# Convert to long format and calculate average CO2 by vertical group
average_CO2 <- GAS.long %>%
        mutate(vertical_group = case_when(
                sampling.point %in% vertical_groups$top ~ "top",
                sampling.point %in% vertical_groups$mid ~ "mid",
                sampling.point %in% vertical_groups$bottom ~ "bottom"
        )) %>%
        group_by(vertical_group, DATE.TIME) %>%
        summarise(avg_CO2 = mean(CO2, na.rm = TRUE))

# Plotting average CO2 by vertical group
ggplot(average_CO2, aes(x = DATE.TIME, y = avg_CO2, color = vertical_group, group = vertical_group)) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        labs(
                title = "Average CO2 Concentration by Vertical Group",
                x = "Date/Time",
                y = "Average CO2 Concentration",
                color = "Vertical Group"
        ) +
        theme_minimal() +
        theme(
                plot.title = element_text(size = 16, face = "bold"),
                axis.title = element_text(size = 14),
                legend.position = "bottom"
        )













                


