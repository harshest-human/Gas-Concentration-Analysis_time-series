getwd()
######## Load library #########
library(tidyverse)
library(reshape2)
library(lubridate)
library(psych)
library(ggplot2)
library(dplyr)
library(data.table)

########### Import combined dataframe of all four gas analysers ############
GAS.long <- fread("2024_Sep_Oct_GAS.long.csv")

######### Calculate hourly averages for Emission values ##########
GAS.long.DT <- GAS.long %>% arrange(DATE.TIME)

count.52 <- GAS.long %>%
        group_by(sampling.point) %>%
        summarise(count = n())

GAS.long.DT <- GAS.long %>% filter(sampling.point != 27) #removed 27th because the CRDS.P8 had error and hence no observation in Sep and October

# Hourly average for sampling points inside
GAS.long_in <- GAS.long.DT %>% 
        filter(sampling.point != 52)

GAS.long_in <- GAS.long_in %>%
        mutate(hour = floor_date(DATE.TIME, unit = "hour")) %>%
        group_by(hour) %>%
        summarise(
                CO2.in = mean(CO2, na.rm = TRUE),
                CH4.in = mean(CH4, na.rm = TRUE),
                NH3.in = mean(NH3, na.rm = TRUE),
                H2O.in = mean(H2O, na.rm = TRUE),
                .groups = "drop")


# Hourly average for sampling points outside
GAS.long_out <- GAS.long.DT %>% 
        filter(sampling.point == 52)

GAS.long_out <- GAS.long_out %>%
        mutate(hour = floor_date(DATE.TIME, unit = "hour")) %>%
        group_by(hour) %>%
        summarise(
                CO2.out = mean(CO2, na.rm = TRUE),
                CH4.out = mean(CH4, na.rm = TRUE),
                NH3.out = mean(NH3, na.rm = TRUE),
                H2O.out = mean(H2O, na.rm = TRUE),
                .groups = "drop")

# convert into data.table
data.table::setDT(GAS.long_in)
data.table::setDT(GAS.long_out)

# combine FTIR and CRDS
GAS.in_out <- GAS.long_in[GAS.long_out, on = .(hour), roll = "nearest"]

# write the combined dataframe
write.csv(GAS.in_out, "2024_Sep_Oct_GAS.in_out.csv", row.names = FALSE)



