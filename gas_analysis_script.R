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
library(gganimate)


######## Import Gas Data #########
FTIR.comb <- fread("2024_June_FTIR.comb.csv")
CRDS.comb <- fread("2024_June_CRDS.comb.csv")


# convert into data.table
data.table::setDT(FTIR.comb)
data.table::setDT(CRDS.comb)

# combine FTIR and CRDS
GAS.comb <- FTIR.comb[CRDS.comb, on = .(DATE.TIME), roll = "nearest"]

# write the combined dataframe
write.csv(GAS.comb, "2024_June_GAS.comb.csv", row.names = FALSE)

######## Data reshaping ##########
# Import the final combined dataframe
GAS.comb <- fread("2024_June_GAS.comb.csv")
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

# Count observations for each unique position 
count.P8 <- GAS.comb %>%
        group_by(MPVPosition.P8) %>%
        summarise(count = n())

count.P9 <- GAS.comb %>%
        group_by(MPVPosition.P9) %>%
        summarise(count = n())

count.F1 <- GAS.comb %>%
        group_by(Messstelle.F1) %>%
        summarise(count = n())

count.F2 <- GAS.comb %>%
        group_by(Messstelle.F2) %>%
        summarise(count = n())

# Convert columns to factors with specified levels and labels
GAS.comb$MPVPosition.P9 <- factor(GAS.comb$MPVPosition.P9, levels = 1:16, labels = 1:16)
GAS.comb$Messstelle.F2 <- factor(GAS.comb$Messstelle.F2, levels = 1:10, labels = 17:26)
GAS.comb$MPVPosition.P8 <- factor(GAS.comb$MPVPosition.P8, levels = 1:16, labels = 27:42)
GAS.comb$Messstelle.F1 <- factor(GAS.comb$Messstelle.F1, levels = 1:10, labels = 43:52)


# Create a new dataframe with the desired structure
GAS.long <- data.table(
        DATE.TIME = GAS.comb$DATE.TIME,
        ID = rep(c("MPVPosition.P9", "Messstelle.F2", "MPVPosition.P8", "Messstelle.F1"), each = nrow(GAS.comb)),
        sampling.point = c(GAS.comb$MPVPosition.P9, GAS.comb$Messstelle.F2, GAS.comb$MPVPosition.P8, GAS.comb$Messstelle.F1),
        CO2 = c(GAS.comb$CO2.P9, GAS.comb$CO2.F2, GAS.comb$CO2.P8, GAS.comb$CO2.F1),
        CH4 = c(GAS.comb$CH4.P9, GAS.comb$CH4.F2, GAS.comb$CH4.P8, GAS.comb$CH4.F1),
        NH3 = c(GAS.comb$NH3.P9, GAS.comb$NH3.F2, GAS.comb$NH3.P8, GAS.comb$NH3.F1),
        H2O = c(GAS.comb$H2O.P9, GAS.comb$H2O.F2, GAS.comb$H2O.P8, GAS.comb$H2O.F1))

# Convert 'GAS.long' to data.table
setDT(GAS.long)

# write after arranging columns
GAS.long <- GAS.long[, .(DATE.TIME, ID, sampling.point, CO2, CH4, NH3, H2O)]
setorder(GAS.long, DATE.TIME)

write.csv(GAS.long, "2024_June_GAS.long.csv", row.names = FALSE)


####### Data Analysis ########
GAS.long <- fread("2024_June_GAS.long.csv")

#count.52 <- GAS.long %>% group_by(sampling.point) %>% summarise(count = n())

setDT(GAS.long)
GAS.long <- GAS.long[sampling.point != 52]
GAS.long$sampling.point <- as.factor(GAS.long$sampling.point)


# Calculate coefficient of variation (CV) at each sampling point
GAS.long[, CV_CO2 := (sd(CO2, na.rm = TRUE) / mean(CO2, na.rm = TRUE)) / 51, by = sampling.point]
GAS.long[, CV_CH4 := (sd(CH4, na.rm = TRUE) / mean(CH4, na.rm = TRUE)) / 51, by = sampling.point]
GAS.long[, CV_NH3 := (sd(NH3, na.rm = TRUE) / mean(NH3, na.rm = TRUE)) / 51, by = sampling.point]
GAS.long[, CV_H2O := (sd(H2O, na.rm = TRUE) / mean(H2O, na.rm = TRUE)) / 51, by = sampling.point]

# Calculate ratio
GAS.long[, ratio_NH3_CO2 := (NH3 / CO2) * 1000]


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


######## Plotting relative error between sampling points #############
# Plot CO2 standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = Err_CO2, fill = vertical)) +
        geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
        geom_point(stat = "summary", fun = "mean", size = 3, shape = 21) +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        labs(x = "Sampling Point", y = " CO2 Relative Error (%)") +
        scale_fill_manual(values = point_fill) +
        scale_y_continuous(limits = c(-60, 60), breaks = seq(-60, 60, by = 10)) +
        theme_minimal() + guides(fill = FALSE) + geom_hline(yintercept = 0, linetype = "dashed", color = "red") 


# Plot CH4 standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = Err_CH4, fill = vertical)) +
        geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
        geom_point(stat = "summary", fun = "mean", size = 3, shape = 21) +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        labs(x = "Sampling Point", y = " CH4 Relative Error (%)") +
        scale_fill_manual(values = point_fill) +
        scale_y_continuous(limits = c(-60, 60), breaks = seq(-60, 60, by = 10)) +
        theme_minimal() + guides(fill = FALSE) + geom_hline(yintercept = 0, linetype = "dashed", color = "red")     

# Plot NH3 standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = Err_NH3, fill = vertical)) +
        geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
        geom_point(stat = "summary", fun = "mean", size = 3, shape = 21) +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        labs(x = "Sampling Point", y = "NH3 Relative Error (%)") +
        scale_fill_manual(values = point_fill) +
        scale_y_continuous(limits = c(-60, 60), breaks = seq(-60, 60, by = 10)) +
        theme_minimal() + guides(fill = FALSE) + geom_hline(yintercept = 0, linetype = "dashed", color = "red")    


######## Plotting coefficient of variation (CV) at each sampling point #############
# Plot CO2 standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = CV_CO2, fill = vertical)) +
        geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
        geom_point(stat = "summary", fun = "mean", size = 3, shape = 21) +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        scale_fill_manual(values = point_fill) +
        theme_minimal() + guides(fill = FALSE) 

# Plot CH4 standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = CV_CH4, fill = vertical)) +
        geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
        geom_point(stat = "summary", fun = "mean", size = 3, shape = 21) +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        scale_fill_manual(values = point_fill) +
        theme_minimal() + guides(fill = FALSE) 

# Plot NH3 standard error bar
ggplot(GAS.long, aes(x = sampling.point, y = CV_NH3, fill = vertical)) +
        geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
        geom_point(stat = "summary", fun = "mean", size = 3, shape = 21) +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        scale_fill_manual(values = point_fill) +
        theme_minimal() + guides(fill = FALSE) 


######## Plotting diel variations #############
# Calculate Hour of Day
GAS.long$Hour <- hour(GAS.long$DATE.TIME)

# Aggregate data by hour for CO2 (similarly for CH4, NH3)
hourly_summary <- GAS.long %>%
        group_by(sampling.point, Hour) %>%
        summarise(mean_CO2 = mean(CO2),
                  mean_CH4 = mean(CH4),
                  mean_NH3 = mean(NH3))

# Reshape data into long format
hourly_summary_long <- hourly_summary %>%
        pivot_longer(cols = starts_with("mean_"), names_to = "Gas", values_to = "Mean_Concentration")

# Plot diel variation for all gases with faceting
ggplot(hourly_summary_long, aes(x = Hour, y = Mean_Concentration, group = sampling.point, color = sampling.point)) +
        geom_line() +
        geom_point() +
        labs(x = "Hour of Day", y = "Mean Concentration", 
             title = "Diel Variation in Gas Concentrations by Sampling Point") +
        scale_x_continuous(breaks = seq(0, 23, by = 1), labels = sprintf("%02d:00", seq(0, 23, by = 1))) +
        facet_wrap(~Gas, scales = "free_y", ncol = 1, labeller = as_labeller(c(mean_CO2 = "CO2", mean_CH4 = "CH4", mean_NH3 = "NH3"))) +
        theme_minimal() +
        theme(legend.position = "none")

        

########## Statistical tests ###########
# Normailty 
hist(GAS.long$CO2, main="Histogram of CO2")
hist(GAS.long$CH4, main="Histogram of CH4")
hist(GAS.long$NH3, main="Histogram of NH3")


# Perform ANOVA 
summary(aov(CV_CO2 ~ sampling.point, data = GAS.long))
summary(aov(CV_CH4 ~ sampling.point, data = GAS.long))
summary(aov(CV_NH3 ~ sampling.point, data = GAS.long))

# Perform Linear Regression
summary(lm(CV_CO2 ~ sampling.point, data = GAS.long))
summary(lm(CV_CH4 ~ sampling.point, data = GAS.long))
summary(lm(CV_NH3 ~ sampling.point, data = GAS.long))


######## CRDS vs FTIR TEST #########
# Import Gas Data
FTIR.test <- fread("D:/Data Analysis/GasmetCX4000_FTIR_Gas_Measurement/FTIR2_test.csv")
CRDS.test <- fread("D:/Data Analysis/Picarro-G2508_CRDS_gas_measurement/CRDS.test.csv")

# Rename levels
FTIR.test$Messstelle.F2 <- factor(FTIR.test$Messstelle.F2, labels = "FTIR")
CRDS.test$MPVPosition.P8 <- factor(CRDS.test$MPVPosition.P8, labels = "CRDS")

# convert into data.table
data.table::setDT(FTIR.test)
data.table::setDT(CRDS.test)

# combine FTIR and CRDS
GAS.test <- FTIR.test[CRDS.test, on = .(DATE.TIME), roll = "nearest"]

# Create a new dataframe with the desired structure
GAS.test.long <- data.table(
        DATE.TIME = GAS.test$DATE.TIME,
        sampling.point = c(GAS.test$Messstelle.F2, rep(GAS.test$MPVPosition.P8)),
        CO2 = c(GAS.test$CO2.F2, GAS.test$CO2.P8),
        CH4 = c(GAS.test$CH4.F2, GAS.test$CH4.P8),
        NH3 = c(GAS.test$NH3.F2, GAS.test$NH3.P8),
        H2O = c(GAS.test$H2O.F2, GAS.test$H2O.P8))

GAS.test.long <- GAS.test.long %>%
        filter(DATE.TIME >= "2024-05-15 12:36:00" & DATE.TIME <= "2024-05-16 14:36:00") %>%
        select(DATE.TIME, sampling.point, CO2, CH4, NH3, H2O)


# write the combined dataframe
write.csv(GAS.test.long, "GAS.test.long.csv", row.names = FALSE)

# Plot for CO2 with blue and red colors
ggplot(GAS.test.long, aes(x = DATE.TIME, y = CO2, colour = as.factor(sampling.point))) +
        geom_line(size = 1) +
        scale_colour_manual(values = c("blue3", "green3")) +  # Setting colors to blue and red
        labs(colour = NULL) + guides(colour = FALSE) +
        theme_minimal()

# Plot for CH4 with blue and red colors
ggplot(GAS.test.long, aes(x = DATE.TIME, y = CH4, colour = as.factor(sampling.point))) +
        geom_line(size = 1) +
        scale_colour_manual(values = c("blue3", "green3")) +  # Setting colors to blue and red
        labs(colour = NULL)  + guides(colour = FALSE) +
        theme_minimal()

# Plot for NH3 with blue and red colors
ggplot(GAS.test.long, aes(x = DATE.TIME, y = NH3, colour = as.factor(sampling.point))) +
        geom_line(size = 1) +
        scale_colour_manual(values = c("blue3", "green3")) +  # Setting colors to blue and red
        labs(colour = NULL)  + guides(colour = FALSE) +
        theme_minimal()

# Compute relative error
GAS.error <- GAS.test %>%
        mutate(Err_CO2 = (CO2.P8 - CO2.F2) / CO2.F2,
               Err_CH4 = (CH4.P8 - CH4.F2) / CH4.F2,
               Err_NH3 = (NH3.P8 - NH3.F2) / NH3.F2)

mean(GAS.error$Err_CO2)
mean(GAS.error$Err_CH4)
mean(GAS.error$Err_NH3)

ggplot(GAS.error, aes(x = DATE.TIME)) +
        geom_line(aes(y = Err_CO2, color = "CO2")) +
        geom_line(aes(y = Err_CH4, color = "CH4")) +
        geom_line(aes(y = Err_NH3, color = "NH3")) +
        scale_color_manual(name = "Errors", values = c("CO2" = "blue", "CH4" = "red", "NH3" = "green")) +
        labs(x = "Time", y = "Error") +
        theme_minimal()




