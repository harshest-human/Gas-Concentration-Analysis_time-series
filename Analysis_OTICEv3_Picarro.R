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

######## Import Gas Data #####
CRDS.data <- fread("D:/Data Analysis/Gas_data/Clean_data/CRDS_clean/20240729_Picarro_G2509.csv")

CRDS.data <- CRDS.data %>%
        mutate(DATE.TIME = as.POSIXct(paste(DATE, TIME), format = "%Y-%m-%d %H:%M:%S")) %>%
        select(DATE.TIME, MPVPosition, N2O, CO2, CH4, H2O, NH3)

CRDS.data <- CRDS.data %>%
        filter(DATE.TIME >= "2024-07-30 13:11:00" & DATE.TIME <= "2024-07-30 15:01:00") %>%
        mutate(DATE.TIME = floor_date(DATE.TIME, unit = "2 minutes"))%>%
        group_by(DATE.TIME) %>%
        summarise(
                MPVPosition = mean(MPVPosition, na.rm = TRUE),
                N2O = mean(N2O, na.rm = TRUE),
                CO2 = mean(CO2, na.rm = TRUE),
                CH4 = mean(CH4, na.rm = TRUE),
                H2O = mean(H2O, na.rm = TRUE),
                NH3 = mean(NH3, na.rm = TRUE))


OTICEv3.data <- fread("D:/Data Analysis/Gas_data/Clean_data/OTICE_clean/20240730_OTICEv3_data.csv")
OTICEv3.data$DATE.TIME <- floor_date(OTICEv3.data$DATE.TIME, unit = "2 minutes")
OTICEv3.data$Node <- as.factor(OTICEv3.data$Node)


####### Plotting ######
# Create individual plots
plot_OTICE <- ggline(OTICEv3.data, x = "DATE.TIME", y = "NH3",
                     add = "mean_se",
                     title = "NH3 Levels by Node (OTICEv3)",
                     xlab = "Time",
                     ylab = "NH3 Concentration",
                     color = "Node")

plot_CRDS <- ggline(CRDS.data, x = "DATE.TIME", y = "NH3",
                    add = "mean_se",
                    title = "NH3 Levels by Node (CRDS)",
                    xlab = "Time",
                    ylab = "NH3 Concentration",
                    color = "black")

# Combine the plots
ggarrange(plot_OTICE, plot_CRDS, ncol = 1, nrow = 2)
