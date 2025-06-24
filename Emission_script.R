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
source("function_indirect.CO2.balance.method.R")

######## Import Data #########
#Load processed gas datasets
LUFA_FTIR_long <- read.csv("20250408-15_long_LUFA_FTIR.csv")
LUFA_FTIR_long$DATE.TIME <- as.POSIXct(LUFA_FTIR_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

ANECO_FTIR_long <- read.csv("20250408-15_long_ANECO_FTIR.csv")
ANECO_FTIR_long$DATE.TIME <- as.POSIXct(ANECO_FTIR_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

MBBM_FTIR_long <- read.csv("20250408-15_long_MBBM_FTIR.csv")
MBBM_FTIR_long$DATE.TIME <- as.POSIXct(MBBM_FTIR_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

ATB_FTIR_long <- read.csv("20250408-15_long_ATB_FTIR.1.csv")
ATB_FTIR_long$DATE.TIME <- as.POSIXct(ATB_FTIR_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

ATB_CRDS_long <- read.csv("20250408-15_long_ATB_CRDS.P8.csv")
ATB_CRDS_long$DATE.TIME <- as.POSIXct(ATB_CRDS_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

LUFA_CRDS_long <- read.csv("20250408-15_long_LUFA_CRDS.P8.csv")
LUFA_CRDS_long$DATE.TIME <- as.POSIXct(LUFA_CRDS_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

UB_CRDS_long <- read.csv("20250408-15_long_UB_CRDS.P8.csv")
UB_CRDS_long$DATE.TIME <- as.POSIXct(UB_CRDS_long$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")

#load animal and temperature data
animal_temp <- read.csv("20250408-15_LVAT_Animal_Temp_data.csv")
animal_temp$DATE.TIME <- as.POSIXct(animal_temp$DATE.TIME, format = "%Y-%m-%d %H:%M:%S")


# Join animal & temp data with gas datasets
emission_LUFA_FTIR  <- left_join(animal_temp, LUFA_FTIR_long,  by = "DATE.TIME")
emission_ANECO_FTIR <- left_join(animal_temp, ANECO_FTIR_long, by = "DATE.TIME")
emission_MBBM_FTIR  <- left_join(animal_temp, MBBM_FTIR_long,  by = "DATE.TIME")
emission_ATB_FTIR   <- left_join(animal_temp, ATB_FTIR_long,   by = "DATE.TIME")
emission_ATB_CRDS   <- left_join(animal_temp, ATB_CRDS_long,   by = "DATE.TIME")
emission_LUFA_CRDS  <- left_join(animal_temp, LUFA_CRDS_long,  by = "DATE.TIME")
emission_UB_CRDS    <- left_join(animal_temp, UB_CRDS_long,    by = "DATE.TIME")


########### Computation of ventilation rates and emissions #################
# Calculate emissions using the function
result_emission_LUFA_FTIR  <- indirect.CO2.balance.method(emission_LUFA_FTIR)
result_emission_ANECO_FTIR <- indirect.CO2.balance.method(emission_ANECO_FTIR)
result_emission_MBBM_FTIR  <- indirect.CO2.balance.method(emission_MBBM_FTIR)
result_emission_ATB_FTIR   <- indirect.CO2.balance.method(emission_ATB_FTIR)
result_emission_ATB_CRDS   <- indirect.CO2.balance.method(emission_ATB_CRDS)
result_emission_LUFA_CRDS  <- indirect.CO2.balance.method(emission_LUFA_CRDS)
result_emission_UB_CRDS    <- indirect.CO2.balance.method(emission_UB_CRDS)

# Write emissions results to CSV
write.csv(result_emission_LUFA_FTIR,  "result_emission_LUFA_FTIR.csv",  row.names = FALSE)
write.csv(result_emission_ANECO_FTIR, "result_emission_ANECO_FTIR.csv", row.names = FALSE)
write.csv(result_emission_MBBM_FTIR,  "result_emission_MBBM_FTIR.csv",  row.names = FALSE)
write.csv(result_emission_ATB_FTIR,   "result_emission_ATB_FTIR.csv",   row.names = FALSE)
write.csv(result_emission_ATB_CRDS,   "result_emission_ATB_CRDS.csv",   row.names = FALSE)
write.csv(result_emission_LUFA_CRDS,  "result_emission_LUFA_CRDS.csv",  row.names = FALSE)
write.csv(result_emission_UB_CRDS,    "result_emission_UB_CRDS.csv",    row.names = FALSE)

# Extract NH3 and CH4 emission columns from each result
result_emission_LUFA_FTIR <- result_emission_LUFA_FTIR %>%
        select(DATE.TIME, matches("^e_NH3"), matches("^e_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

result_emission_ANECO_FTIR <- result_emission_ANECO_FTIR %>%
        select(DATE.TIME, matches("^e_NH3"), matches("^e_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

result_emission_MBBM_FTIR <- result_emission_MBBM_FTIR %>%
        select(DATE.TIME, matches("^e_NH3"), matches("^e_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

result_emission_ATB_FTIR <- result_emission_ATB_FTIR %>%
        select(DATE.TIME, matches("^e_NH3"), matches("^e_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

result_emission_ATB_CRDS <- result_emission_ATB_CRDS %>%
        select(DATE.TIME, matches("^e_NH3"), matches("^e_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

result_emission_LUFA_CRDS <- result_emission_LUFA_CRDS %>%
        select(DATE.TIME, matches("^e_NH3"), matches("^e_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

result_emission_UB_CRDS <- result_emission_UB_CRDS %>%
        select(DATE.TIME, matches("^e_NH3"), matches("^e_CH4")) %>%
        distinct(DATE.TIME, .keep_all = TRUE)

########### Cleaned hourly emissions #################
# Combine for all result dataframes
e_CH4 <- final_emission_combined %>%
        select(DATE.TIME, contains("e_CH4_")) %>%
        select(-contains("_per_year")) %>%
        pivot_longer(
                cols = -DATE.TIME,
                names_to = "lab",
                values_to = "e_CH4"
        ) %>%
        mutate(bg_direction = case_when(
                        str_detect(lab, "_N_") ~ "North",
                        str_detect(lab, "_S_") ~ "South",
                        TRUE ~ "Unknown"
                ),
                lab.analyzer = case_when(
                        str_detect(lab, "LUFA_FTIR") ~ "LUFA_FTIR",
                        str_detect(lab, "ANECO_FTIR") ~ "ANECO_FTIR",
                        str_detect(lab, "ATB_FTIR.1") ~ "ATB_FTIR.1",
                        str_detect(lab, "MBM_FTIR") ~ "MBBM_FTIR",
                        str_detect(lab, "ATB_CRDS.P8") ~ "ATB_CRDS.P8",
                        str_detect(lab, "UB_CRDS.P8") ~ "UB_CRDS.P8",
                        str_detect(lab, "LUFA_CRDS.P8") ~ "LUFA_CRDS.P8",
                        TRUE ~ "Unknown"
                ),
                hour = as.factor(hour(DATE.TIME))
        )

e_CH4 <- e_CH4 %>% select(-lab)
        
e_NH3 <- final_emission_combined %>%
        select(DATE.TIME, contains("e_NH3_")) %>%
        select(-contains("_per_year")) %>%
        pivot_longer(
                cols = -DATE.TIME,
                names_to = "lab",
                values_to = "e_NH3"
        ) %>%
        mutate(bg_direction = case_when(
                str_detect(lab, "_N_") ~ "North",
                str_detect(lab, "_S_") ~ "South",
                TRUE ~ "Unknown"
        ),
        lab.analyzer = case_when(
                str_detect(lab, "LUFA_FTIR") ~ "LUFA_FTIR",
                str_detect(lab, "ANECO_FTIR") ~ "ANECO_FTIR",
                str_detect(lab, "ATB_FTIR.1") ~ "ATB_FTIR.1",
                str_detect(lab, "MBM_FTIR") ~ "MBBM_FTIR",
                str_detect(lab, "ATB_CRDS.P8") ~ "ATB_CRDS.P8",
                str_detect(lab, "UB_CRDS.P8") ~ "UB_CRDS.P8",
                str_detect(lab, "LUFA_CRDS.P8") ~ "LUFA_CRDS.P8",
                TRUE ~ "Unknown"
        ),
        hour = as.factor(hour(DATE.TIME))
        )

e_NH3 <- e_NH3 %>% select(-lab)

# --- Plotting function with 'background' ---
plot_emission <- function(df, x, y, bg_direction) {
        library(dplyr)
        library(ggpubr)
        library(ggplot2)
        library(scales)
        
        df_filtered <- df %>%
                filter(bg_direction == bg_direction)
        
        p <- ggline(df_filtered, x = x, y = y,
                    add = "mean_se",
                    color = "lab.analyzer") +
                labs(title = paste(y, "Emission Trends (mean ± SE) – Background:", bg_direction),
                     x = x,
                     y = paste(y, "Emission (g/h)"),
                     color = "Laboratory") +
                scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
                theme_light() +
                theme(legend.position = "top")
        
        # If x is datetime, add scale_x_datetime and rotate labels
        if (inherits(df[[x]], "POSIXct") || inherits(df[[x]], "POSIXt")) {
                p <- p +
                        scale_x_datetime(date_breaks = "6 hours", date_labels = "%d.%m %H:%M") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
        
        return(p)
}



# --- Generate DATE.TIME plots ---
plot_NH3_north <- plot_emission(e_NH3, x = "DATE.TIME", y = "e_NH3", bg_direction = "North")
plot_NH3_south <- plot_emission(e_NH3, x = "DATE.TIME", y = "e_NH3", bg_direction = "South")
plot_CH4_north <- plot_emission(e_CH4, x = "DATE.TIME", y = "e_CH4", bg_direction = "North")
plot_CH4_south <- plot_emission(e_CH4, x = "DATE.TIME", y = "e_CH4", bg_direction = "South")

# --- View DATE.TIME plots ---
plot_NH3_north
plot_NH3_south
plot_CH4_north
plot_CH4_south


# --- Generate hour-based plots ---
plot_NH3_north_hour <- plot_emission(e_NH3, x = "hour", y = "e_NH3", bg_direction = "North")
plot_NH3_south_hour <- plot_emission(e_NH3, x = "hour", y = "e_NH3", bg_direction = "South")
plot_CH4_north_hour <- plot_emission(e_CH4, x = "hour", y = "e_CH4", bg_direction = "North")
plot_CH4_south_hour <- plot_emission(e_CH4, x = "hour", y = "e_CH4", bg_direction = "South")

# --- View hour-based plots ---
plot_NH3_north_hour
plot_NH3_south_hour
plot_CH4_north_hour
plot_CH4_south_hour

# Define named list of emission plots
emission_plots <- list(
        plot_NH3_north = plot_NH3_north,
        plot_NH3_south = plot_NH3_south,
        plot_CH4_north_hour = plot_CH4_north_hour,
        plot_CH4_south_hour = plot_CH4_south_hour
)

# Save each emission plot as high-res PNG
for (name in names(emission_plots)) {
        ggsave(
                filename = paste0(name, ".png"),
                plot = emission_plots[[name]],
                width = 14,
                height = 8,
                dpi = 600
        )
}

# Read the emission PNGs and combine into a single PDF
png_files <- paste0(names(emission_plots), ".png")
img_list <- magick::image_read(png_files)
magick::image_write(image = img_list, path = "Ringversuche_emission_plots.pdf", format = "pdf")

