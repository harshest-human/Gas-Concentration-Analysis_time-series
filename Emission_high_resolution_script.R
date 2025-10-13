getwd()
######## Load library #########
library(tidyverse)
library(reshape2)
library(lubridate)
library(psych)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(scales)
library(gridExtra)
library(magick)
library(summarytools)
library(rlang)
library(DescTools)
library(rstatix)
library(ggcorrplot)
library(patchwork)
library(networkD3)
library(ggridges)
library(rstatix)
library(multcompView)
library(viridis)
library(lme4)
library(patchwork)
library(kableExtra)
library(knitr)

###### Import Data ########
animal_temp_data <- read.csv("20250928-20250930_LVAT_Animal_Temp_data.csv")
gas_data <- read.csv("20251010_high_resolution_gas_concentration_data.csv")

# Define colors for vertical groups
point_fill <- c("top" = "orange", "mid" = "green3", "bottom" = "steelblue1")


############# DATA ANALYSIS ##########

############# DATA Visualization ##########
ggplot(gas_data  %>% filter(DATE.TIME >= "2024-07-31 11:05:58"), aes(x = location, y = NH3, fill = vgroup)) +
        geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
        geom_point(stat = "summary", fun = "mean", size = 3, shape = 21) +
        geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
        scale_fill_manual(values = point_fill) +
        theme_minimal() + guides(fill = FALSE)
