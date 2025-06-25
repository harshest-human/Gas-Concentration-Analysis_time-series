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

######## Import Gas Data #########
#Load processed datasets
LUFA_FTIR <- read.csv("20250408-15_hourly_LUFA_FTIR.csv")
ANECO_FTIR <- read.csv("20250408-15_hourly_ANECO_FTIR.csv")
MBBM_FTIR <- read.csv("20250408-15_hourly_MBBM_FTIR.csv")
ATB_FTIR <- read.csv("20250408-15_hourly_ATB_FTIR.1.csv")
ATB_CRDS <- read.csv("20250408-15_hourly_ATB_CRDS.P8.csv")
LUFA_CRDS <- read.csv("20250408-15_hourly_LUFA_CRDS.P8.csv")
UB_CRDS <- read.csv("20250408-15_hourly_UB_CRDS.P8.csv")

#combine all data set
lab_combined <- bind_rows(LUFA_FTIR, ANECO_FTIR, MBBM_FTIR, ATB_FTIR, ATB_CRDS, LUFA_CRDS, UB_CRDS)
lab_combined$lab.analyzer <- paste(lab_combined$lab, lab_combined$analyzer, sep = "_")
lab_combined$DATE.TIME <- as.POSIXct(lab_combined$DATE.TIME)
lab_combined <- lab_combined %>% select(DATE.TIME, location, lab.analyzer, CO2, CH4, NH3)
lab_combined <- lab_combined %>% filter(DATE.TIME >= "2025-04-08 12:00:00" & DATE.TIME <= "2025-04-15 12:00:00")
lab_combined <- lab_combined %>% mutate(hour = hour(DATE.TIME))
lab_combined$hour <- as.factor(lab_combined$hour)

#write csv
lab_combined <- lab_combined %>% select(DATE.TIME, hour, everything())
write.csv(lab_combined, "20250408-15_final_ringversuche_concentration_combined_data.csv", row.names = FALSE)

####### Data Visualization Weekly ##########
plot_concentration <- function(df, x, y) {
        #load libraries
        library(dplyr)
        library(ggpubr)
        library(ggplot2)
        library(scales)
        
        # Check if x and y exist in the dataframe
        if (!(x %in% colnames(df))) stop(paste("X-axis column", x, "not found in dataframe."))
        if (!(y %in% colnames(df))) stop(paste("Y-axis column", y, "not found in dataframe."))
        
        # Filter out NA or non-finite values from y
        df <- df %>% filter(!is.na(.data[[y]]), is.finite(.data[[y]]))
        
        # Start building plot
        p <- ggline(df,
                    x = x,
                    y = y,
                    add = "mean_se",
                    color = "lab.analyzer",
                    facet.by = "location") +
                labs(
                        title = paste(y, "Concentration Trends (Mean ± SE)"),
                        x = x,
                        y = paste(y, "(ppmv)"),
                        color = "Laboratory"
                ) +
                theme_light() +
                theme(legend.position = "top")
        
        # If x is datetime, add datetime scale and rotate x-axis labels
        if (inherits(df[[x]], "POSIXct")) {
                p <- p +
                        scale_x_datetime(date_breaks = "6 hours", date_labels = "%d.%m %H:%M") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
        
        return(p)
}


# --- Generate DATE.TIME plots ---
plot_CO2_conc_all <- plot_concentration(lab_combined, x = "DATE.TIME", y = "CO2")
plot_CH4_conc_all <- plot_concentration(lab_combined, x = "DATE.TIME", y = "CH4")
plot_NH3_conc_all <- plot_concentration(lab_combined, x = "DATE.TIME", y = "NH3")

# --- View DATE.TIME plots ---
plot_CO2_conc_all
plot_CH4_conc_all
plot_NH3_conc_all


# --- Generate hour-based plots ---
plot_CO2_conc_hour <- plot_concentration(lab_combined, x = "hour", y = "CO2")
plot_CH4_conc_hour <- plot_concentration(lab_combined, x = "hour", y = "CH4")
plot_NH3_conc_hour <- plot_concentration(lab_combined, x = "hour", y = "NH3")

# --- View hour-based plots ---
plot_CO2_conc_hour
plot_CH4_conc_hour
plot_NH3_conc_hour


# Define a named list of concentration plots
plots <- list(
        plot_CO2_conc_all  = plot_CO2_conc_all,
        plot_CH4_conc_all  = plot_CH4_conc_all,
        plot_NH3_conc_all  = plot_NH3_conc_all,
        plot_CO2_conc_hour = plot_CO2_conc_hour,
        plot_CH4_conc_hour = plot_CH4_conc_hour,
        plot_NH3_conc_hour = plot_NH3_conc_hour
)


# Save each plot as high-res PNG
for (name in names(plots)) {
        ggsave(
                filename = paste0(name, ".png"),
                plot = plots[[name]],
                width = 14,
                height = 8,
                dpi = 600
        )
}

# Combine PNGs into PDF using magick
png_files <- paste0(names(plots), ".png")
img_list <- magick::image_read(png_files)
magick::image_write(image = img_list, path = "Ringversuche_concentration_plots.pdf", format = "pdf")


######### Statistical analysis ########
# Read each dataset, convert DATE.TIME and rename columns (except DATE.TIME) with suffix:
LUFA_FTIR_long <- read.csv("20250408-15_long_LUFA_FTIR.csv") %>%
        mutate(DATE.TIME = as.POSIXct(DATE.TIME, format = "%Y-%m-%d %H:%M:%S")) %>%
        rename_with(~paste0(., "_FTIR"), -DATE.TIME)

ANECO_FTIR_long <- read.csv("20250408-15_long_ANECO_FTIR.csv") %>%
        mutate(DATE.TIME = as.POSIXct(DATE.TIME, format = "%Y-%m-%d %H:%M:%S")) %>%
        rename_with(~paste0(., "_FTIR"), -DATE.TIME)

MBBM_FTIR_long <- read.csv("20250408-15_long_MBBM_FTIR.csv") %>%
        mutate(DATE.TIME = as.POSIXct(DATE.TIME, format = "%Y-%m-%d %H:%M:%S")) %>%
        rename_with(~paste0(., "_FTIR"), -DATE.TIME)

ATB_FTIR_long <- read.csv("20250408-15_long_ATB_FTIR.1.csv") %>%
        mutate(DATE.TIME = as.POSIXct(DATE.TIME, format = "%Y-%m-%d %H:%M:%S")) %>%
        rename_with(~paste0(., "_FTIR.1"), -DATE.TIME)

ATB_CRDS_long <- read.csv("20250408-15_long_ATB_CRDS.P8.csv") %>%
        mutate(DATE.TIME = as.POSIXct(DATE.TIME, format = "%Y-%m-%d %H:%M:%S")) %>%
        rename_with(~paste0(., "_CRDS.P8"), -DATE.TIME)

LUFA_CRDS_long <- read.csv("20250408-15_long_LUFA_CRDS.P8.csv") %>%
        mutate(DATE.TIME = as.POSIXct(DATE.TIME, format = "%Y-%m-%d %H:%M:%S")) %>%
        rename_with(~paste0(., "_CRDS.P8"), -DATE.TIME)

UB_CRDS_long <- read.csv("20250408-15_long_UB_CRDS.P8.csv") %>%
        mutate(DATE.TIME = as.POSIXct(DATE.TIME, format = "%Y-%m-%d %H:%M:%S")) %>%
        rename_with(~paste0(., "_CRDS.P8"), -DATE.TIME)


final_long_combined <- full_join(LUFA_FTIR_long, ANECO_FTIR_long, by = "DATE.TIME") %>%
        full_join(MBBM_FTIR_long, by = "DATE.TIME") %>%
        full_join(ATB_FTIR_long, by = "DATE.TIME") %>%
        full_join(ATB_CRDS_long, by = "DATE.TIME") %>%
        full_join(LUFA_CRDS_long, by = "DATE.TIME") %>%
        full_join(UB_CRDS_long, by = "DATE.TIME") %>%
        select(-contains("analyzer"))

write.csv(final_long_combined, "20250408-15_final_concentration_long_combined.csv",
          row.names = FALSE)


# Calculate CO2 emissions percentage errors relative to ATB_CRDS.P8
CO2_err <- final_long_combined %>%
        mutate(
                CO2_in_LUFA_FTIR_err    = 100 * (CO2_in_LUFA_FTIR     - CO2_in_ATB_CRDS.P8) / CO2_in_ATB_CRDS.P8,
                CO2_in_ANECO_FTIR_err   = 100 * (CO2_in_ANECO_FTIR    - CO2_in_ATB_CRDS.P8) / CO2_in_ATB_CRDS.P8,
                CO2_in_MBBM_FTIR_err    = 100 * (CO2_in_MBBM_FTIR     - CO2_in_ATB_CRDS.P8) / CO2_in_ATB_CRDS.P8,
                CO2_in_ATB_FTIR.1_err   = 100 * (CO2_in_ATB_FTIR.1    - CO2_in_ATB_CRDS.P8) / CO2_in_ATB_CRDS.P8,
                CO2_in_UB_CRDS.P8_err   = 100 * (CO2_in_UB_CRDS.P8    - CO2_in_ATB_CRDS.P8) / CO2_in_ATB_CRDS.P8,
                CO2_in_LUFA_CRDS.P8_err = 100 * (CO2_in_LUFA_CRDS.P8  - CO2_in_ATB_CRDS.P8) / CO2_in_ATB_CRDS.P8,
                CO2_N_LUFA_FTIR_err    = 100 * (CO2_N_LUFA_FTIR     - CO2_N_ATB_CRDS.P8) / CO2_N_ATB_CRDS.P8,
                CO2_N_ANECO_FTIR_err   = 100 * (CO2_N_ANECO_FTIR    - CO2_N_ATB_CRDS.P8) / CO2_N_ATB_CRDS.P8,
                CO2_N_MBBM_FTIR_err    = 100 * (CO2_N_MBBM_FTIR     - CO2_N_ATB_CRDS.P8) / CO2_N_ATB_CRDS.P8,
                CO2_N_ATB_FTIR.1_err   = 100 * (CO2_N_ATB_FTIR.1    - CO2_N_ATB_CRDS.P8) / CO2_N_ATB_CRDS.P8,
                CO2_N_UB_CRDS.P8_err   = 100 * (CO2_N_UB_CRDS.P8    - CO2_N_ATB_CRDS.P8) / CO2_N_ATB_CRDS.P8,
                CO2_N_LUFA_CRDS.P8_err = 100 * (CO2_N_LUFA_CRDS.P8  - CO2_N_ATB_CRDS.P8) / CO2_N_ATB_CRDS.P8,
                CO2_S_LUFA_FTIR_err    = 100 * (CO2_S_LUFA_FTIR     - CO2_S_ATB_CRDS.P8) / CO2_S_ATB_CRDS.P8,
                CO2_S_ANECO_FTIR_err   = 100 * (CO2_S_ANECO_FTIR    - CO2_S_ATB_CRDS.P8) / CO2_S_ATB_CRDS.P8,
                CO2_S_MBBM_FTIR_err    = 100 * (CO2_S_MBBM_FTIR     - CO2_S_ATB_CRDS.P8) / CO2_S_ATB_CRDS.P8,
                CO2_S_ATB_FTIR.1_err   = 100 * (CO2_S_ATB_FTIR.1    - CO2_S_ATB_CRDS.P8) / CO2_S_ATB_CRDS.P8,
                CO2_S_UB_CRDS.P8_err   = 100 * (CO2_S_UB_CRDS.P8    - CO2_S_ATB_CRDS.P8) / CO2_S_ATB_CRDS.P8,
                CO2_S_LUFA_CRDS.P8_err = 100 * (CO2_S_LUFA_CRDS.P8  - CO2_S_ATB_CRDS.P8) / CO2_S_ATB_CRDS.P8
        ) %>%
        select(
                DATE.TIME,
                CO2_in_LUFA_FTIR_err, 
                CO2_in_ANECO_FTIR_err, 
                CO2_in_MBBM_FTIR_err, 
                CO2_in_ATB_FTIR.1_err,
                CO2_in_UB_CRDS.P8_err,
                CO2_in_LUFA_CRDS.P8_err,
                CO2_N_LUFA_FTIR_err, 
                CO2_N_ANECO_FTIR_err, 
                CO2_N_MBBM_FTIR_err, 
                CO2_N_ATB_FTIR.1_err,
                CO2_N_UB_CRDS.P8_err,
                CO2_N_LUFA_CRDS.P8_err,
                CO2_S_LUFA_FTIR_err,
                CO2_S_ANECO_FTIR_err,
                CO2_S_MBBM_FTIR_err,
                CO2_S_ATB_FTIR.1_err,
                CO2_S_UB_CRDS.P8_err,
                CO2_S_LUFA_CRDS.P8_err
        )

# Calculate CH4 emissions percentage errors relative to ATB_CRDS.P8
CH4_err <- final_long_combined %>%
        mutate(
                CH4_in_LUFA_FTIR_err    = 100 * (CH4_in_LUFA_FTIR     - CH4_in_ATB_CRDS.P8) / CH4_in_ATB_CRDS.P8,
                CH4_in_ANECO_FTIR_err   = 100 * (CH4_in_ANECO_FTIR    - CH4_in_ATB_CRDS.P8) / CH4_in_ATB_CRDS.P8,
                CH4_in_MBBM_FTIR_err    = 100 * (CH4_in_MBBM_FTIR     - CH4_in_ATB_CRDS.P8) / CH4_in_ATB_CRDS.P8,
                CH4_in_ATB_FTIR.1_err   = 100 * (CH4_in_ATB_FTIR.1    - CH4_in_ATB_CRDS.P8) / CH4_in_ATB_CRDS.P8,
                CH4_in_UB_CRDS.P8_err   = 100 * (CH4_in_UB_CRDS.P8    - CH4_in_ATB_CRDS.P8) / CH4_in_ATB_CRDS.P8,
                CH4_in_LUFA_CRDS.P8_err = 100 * (CH4_in_LUFA_CRDS.P8  - CH4_in_ATB_CRDS.P8) / CH4_in_ATB_CRDS.P8,
                CH4_N_LUFA_FTIR_err    = 100 * (CH4_N_LUFA_FTIR     - CH4_N_ATB_CRDS.P8) / CH4_N_ATB_CRDS.P8,
                CH4_N_ANECO_FTIR_err   = 100 * (CH4_N_ANECO_FTIR    - CH4_N_ATB_CRDS.P8) / CH4_N_ATB_CRDS.P8,
                CH4_N_MBBM_FTIR_err    = 100 * (CH4_N_MBBM_FTIR     - CH4_N_ATB_CRDS.P8) / CH4_N_ATB_CRDS.P8,
                CH4_N_ATB_FTIR.1_err   = 100 * (CH4_N_ATB_FTIR.1    - CH4_N_ATB_CRDS.P8) / CH4_N_ATB_CRDS.P8,
                CH4_N_UB_CRDS.P8_err   = 100 * (CH4_N_UB_CRDS.P8    - CH4_N_ATB_CRDS.P8) / CH4_N_ATB_CRDS.P8,
                CH4_N_LUFA_CRDS.P8_err = 100 * (CH4_N_LUFA_CRDS.P8  - CH4_N_ATB_CRDS.P8) / CH4_N_ATB_CRDS.P8,
                CH4_S_LUFA_FTIR_err    = 100 * (CH4_S_LUFA_FTIR     - CH4_S_ATB_CRDS.P8) / CH4_S_ATB_CRDS.P8,
                CH4_S_ANECO_FTIR_err   = 100 * (CH4_S_ANECO_FTIR    - CH4_S_ATB_CRDS.P8) / CH4_S_ATB_CRDS.P8,
                CH4_S_MBBM_FTIR_err    = 100 * (CH4_S_MBBM_FTIR     - CH4_S_ATB_CRDS.P8) / CH4_S_ATB_CRDS.P8,
                CH4_S_ATB_FTIR.1_err   = 100 * (CH4_S_ATB_FTIR.1    - CH4_S_ATB_CRDS.P8) / CH4_S_ATB_CRDS.P8,
                CH4_S_UB_CRDS.P8_err   = 100 * (CH4_S_UB_CRDS.P8    - CH4_S_ATB_CRDS.P8) / CH4_S_ATB_CRDS.P8,
                CH4_S_LUFA_CRDS.P8_err = 100 * (CH4_S_LUFA_CRDS.P8  - CH4_S_ATB_CRDS.P8) / CH4_S_ATB_CRDS.P8
        ) %>%
        select(
                DATE.TIME,
                CH4_in_LUFA_FTIR_err, 
                CH4_in_ANECO_FTIR_err, 
                CH4_in_MBBM_FTIR_err, 
                CH4_in_ATB_FTIR.1_err,
                CH4_in_UB_CRDS.P8_err,
                CH4_in_LUFA_CRDS.P8_err,
                CH4_N_LUFA_FTIR_err, 
                CH4_N_ANECO_FTIR_err, 
                CH4_N_MBBM_FTIR_err, 
                CH4_N_ATB_FTIR.1_err,
                CH4_N_UB_CRDS.P8_err,
                CH4_N_LUFA_CRDS.P8_err,
                CH4_S_LUFA_FTIR_err,
                CH4_S_ANECO_FTIR_err,
                CH4_S_MBBM_FTIR_err,
                CH4_S_ATB_FTIR.1_err,
                CH4_S_UB_CRDS.P8_err,
                CH4_S_LUFA_CRDS.P8_err
        )



# Calculate NH3 emissions percentage errors relative to ATB_CRDS.P8
NH3_err <- final_long_combined %>%
        mutate(
                NH3_in_LUFA_FTIR_err    = 100 * (NH3_in_LUFA_FTIR     - NH3_in_ATB_CRDS.P8) / NH3_in_ATB_CRDS.P8,
                NH3_in_ANECO_FTIR_err   = 100 * (NH3_in_ANECO_FTIR    - NH3_in_ATB_CRDS.P8) / NH3_in_ATB_CRDS.P8,
                NH3_in_MBBM_FTIR_err    = 100 * (NH3_in_MBBM_FTIR     - NH3_in_ATB_CRDS.P8) / NH3_in_ATB_CRDS.P8,
                NH3_in_ATB_FTIR.1_err   = 100 * (NH3_in_ATB_FTIR.1    - NH3_in_ATB_CRDS.P8) / NH3_in_ATB_CRDS.P8,
                NH3_in_UB_CRDS.P8_err   = 100 * (NH3_in_UB_CRDS.P8    - NH3_in_ATB_CRDS.P8) / NH3_in_ATB_CRDS.P8,
                NH3_in_LUFA_CRDS.P8_err = 100 * (NH3_in_LUFA_CRDS.P8  - NH3_in_ATB_CRDS.P8) / NH3_in_ATB_CRDS.P8,
                NH3_N_LUFA_FTIR_err    = 100 * (NH3_N_LUFA_FTIR     - NH3_N_ATB_CRDS.P8) / NH3_N_ATB_CRDS.P8,
                NH3_N_ANECO_FTIR_err   = 100 * (NH3_N_ANECO_FTIR    - NH3_N_ATB_CRDS.P8) / NH3_N_ATB_CRDS.P8,
                NH3_N_MBBM_FTIR_err    = 100 * (NH3_N_MBBM_FTIR     - NH3_N_ATB_CRDS.P8) / NH3_N_ATB_CRDS.P8,
                NH3_N_ATB_FTIR.1_err   = 100 * (NH3_N_ATB_FTIR.1    - NH3_N_ATB_CRDS.P8) / NH3_N_ATB_CRDS.P8,
                NH3_N_UB_CRDS.P8_err   = 100 * (NH3_N_UB_CRDS.P8    - NH3_N_ATB_CRDS.P8) / NH3_N_ATB_CRDS.P8,
                NH3_N_LUFA_CRDS.P8_err = 100 * (NH3_N_LUFA_CRDS.P8  - NH3_N_ATB_CRDS.P8) / NH3_N_ATB_CRDS.P8,
                NH3_S_LUFA_FTIR_err    = 100 * (NH3_S_LUFA_FTIR     - NH3_S_ATB_CRDS.P8) / NH3_S_ATB_CRDS.P8,
                NH3_S_ANECO_FTIR_err   = 100 * (NH3_S_ANECO_FTIR    - NH3_S_ATB_CRDS.P8) / NH3_S_ATB_CRDS.P8,
                NH3_S_MBBM_FTIR_err    = 100 * (NH3_S_MBBM_FTIR     - NH3_S_ATB_CRDS.P8) / NH3_S_ATB_CRDS.P8,
                NH3_S_ATB_FTIR.1_err   = 100 * (NH3_S_ATB_FTIR.1    - NH3_S_ATB_CRDS.P8) / NH3_S_ATB_CRDS.P8,
                NH3_S_UB_CRDS.P8_err   = 100 * (NH3_S_UB_CRDS.P8    - NH3_S_ATB_CRDS.P8) / NH3_S_ATB_CRDS.P8,
                NH3_S_LUFA_CRDS.P8_err = 100 * (NH3_S_LUFA_CRDS.P8  - NH3_S_ATB_CRDS.P8) / NH3_S_ATB_CRDS.P8
        ) %>%
        select(
                DATE.TIME,
                NH3_in_LUFA_FTIR_err, 
                NH3_in_ANECO_FTIR_err, 
                NH3_in_MBBM_FTIR_err, 
                NH3_in_ATB_FTIR.1_err,
                NH3_in_UB_CRDS.P8_err,
                NH3_in_LUFA_CRDS.P8_err,
                NH3_N_LUFA_FTIR_err, 
                NH3_N_ANECO_FTIR_err, 
                NH3_N_MBBM_FTIR_err, 
                NH3_N_ATB_FTIR.1_err,
                NH3_N_UB_CRDS.P8_err,
                NH3_N_LUFA_CRDS.P8_err,
                NH3_S_LUFA_FTIR_err,
                NH3_S_ANECO_FTIR_err,
                NH3_S_MBBM_FTIR_err,
                NH3_S_ATB_FTIR.1_err,
                NH3_S_UB_CRDS.P8_err,
                NH3_S_LUFA_CRDS.P8_err
        )


# Pivot CO2 emission errors longer# Pivot CO2 emission errors longer
CO2_err <- CO2_err %>%
        pivot_longer(
                cols = -DATE.TIME,
                names_to = "lab",
                values_to = "CO2_err"
        ) %>%
        mutate(
                bg_direction = case_when(
                        str_detect(lab, "_in_") ~ "Inside",
                        str_detect(lab, "_N_") ~ "North",
                        str_detect(lab, "_S_") ~ "South",
                        TRUE ~ "Unknown"
                ),
                lab.analyzer = case_when(
                        str_detect(lab, "LUFA_FTIR") ~ "LUFA_FTIR",
                        str_detect(lab, "ANECO_FTIR") ~ "ANECO_FTIR",
                        str_detect(lab, "ATB_FTIR.1") ~ "ATB_FTIR.1",
                        str_detect(lab, "MBBM_FTIR") ~ "MBBM_FTIR",
                        str_detect(lab, "ATB_CRDS.P8") ~ "ATB_CRDS.P8",
                        str_detect(lab, "UB_CRDS.P8") ~ "UB_CRDS.P8",
                        str_detect(lab, "LUFA_CRDS.P8") ~ "LUFA_CRDS.P8",
                        TRUE ~ "Unknown"
                ),
                hour = as.factor(hour(DATE.TIME))
        ) %>%
        select(-lab)


# Pivot CH4 emission errors longer# Pivot CO2 emission errors longer
CH4_err <- CH4_err %>%
        pivot_longer(
                cols = -DATE.TIME,
                names_to = "lab",
                values_to = "CH4_err"
        ) %>%
        mutate(
                bg_direction = case_when(
                        str_detect(lab, "_in_") ~ "Inside",
                        str_detect(lab, "_N_") ~ "North",
                        str_detect(lab, "_S_") ~ "South",
                        TRUE ~ "Unknown"
                ),
                lab.analyzer = case_when(
                        str_detect(lab, "LUFA_FTIR") ~ "LUFA_FTIR",
                        str_detect(lab, "ANECO_FTIR") ~ "ANECO_FTIR",
                        str_detect(lab, "ATB_FTIR.1") ~ "ATB_FTIR.1",
                        str_detect(lab, "MBBM_FTIR") ~ "MBBM_FTIR",
                        str_detect(lab, "ATB_CRDS.P8") ~ "ATB_CRDS.P8",
                        str_detect(lab, "UB_CRDS.P8") ~ "UB_CRDS.P8",
                        str_detect(lab, "LUFA_CRDS.P8") ~ "LUFA_CRDS.P8",
                        TRUE ~ "Unknown"
                ),
                hour = as.factor(hour(DATE.TIME))
        ) %>%
        select(-lab)


# Pivot NH3 emission errors longer
NH3_err <- NH3_err %>%
        pivot_longer(
                cols = -DATE.TIME,
                names_to = "lab",
                values_to = "NH3_err"
        ) %>%
        mutate(
                bg_direction = case_when(
                        str_detect(lab, "_in_") ~ "Inside",
                        str_detect(lab, "_N_") ~ "North",
                        str_detect(lab, "_S_") ~ "South",
                        TRUE ~ "Unknown"
                ),
                lab.analyzer = case_when(
                        str_detect(lab, "LUFA_FTIR") ~ "LUFA_FTIR",
                        str_detect(lab, "ANECO_FTIR") ~ "ANECO_FTIR",
                        str_detect(lab, "ATB_FTIR.1") ~ "ATB_FTIR.1",
                        str_detect(lab, "MBBM_FTIR") ~ "MBBM_FTIR",
                        str_detect(lab, "ATB_CRDS.P8") ~ "ATB_CRDS.P8",
                        str_detect(lab, "UB_CRDS.P8") ~ "UB_CRDS.P8",
                        str_detect(lab, "LUFA_CRDS.P8") ~ "LUFA_CRDS.P8",
                        TRUE ~ "Unknown"
                ),
                hour = as.factor(hour(DATE.TIME))
        ) %>%
        select(-lab)


# calculate mean relative error
avg_CO2_err <- CO2_err %>%
        group_by(bg_direction, lab.analyzer) %>%
        summarise(mean_CO2_err = mean(CO2_err, na.rm = TRUE), .groups = "drop")%>%
        mutate(bar_color = ifelse(mean_CO2_err >= 0, "Positive", "Negative"))

avg_CH4_err <- CH4_err %>%
        group_by(bg_direction, lab.analyzer) %>%
        summarise(mean_CH4_err = mean(CH4_err, na.rm = TRUE), .groups = "drop")%>%
        mutate(bar_color = ifelse(mean_CH4_err >= 0, "Positive", "Negative"))

avg_NH3_err <- NH3_err %>%
        group_by(bg_direction, lab.analyzer) %>%
        summarise(mean_NH3_err = mean(NH3_err, na.rm = TRUE), .groups = "drop")%>%
        mutate(bar_color = ifelse(mean_NH3_err >= 0, "Positive", "Negative"))


# Create a column for bar color
plot_CO2_err <- ggplot(avg_CO2_err, aes(x = lab.analyzer, y = mean_CO2_err, fill = bar_color)) +
        geom_bar(stat = "identity", width = 0.7) +
        geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.7) +
        facet_wrap(~ bg_direction) +
        scale_fill_manual(values = c("Positive" = "lightgreen", "Negative" = "red4")) +
        labs(
                title = "Mean Relative Errors of CO2 Concentrations      (Ref. = ATB_CRDS.P8)",
                x = "Laboratory",
                y = "CO2 Relative Error (%)",
                fill = "Error Sign"
        ) +
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_CH4_err <- ggplot(avg_CH4_err, aes(x = lab.analyzer, y = mean_CH4_err, fill = bar_color)) +
        geom_bar(stat = "identity", width = 0.7) +
        geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.7) +
        facet_wrap(~ bg_direction) +
        scale_fill_manual(values = c("Positive" = "lightgreen", "Negative" = "red4")) +
        labs(
                title = "Mean Relative Errors of CH4 Concentrations      (Ref. = ATB_CRDS.P8)",
                x = "Laboratory",
                y = "CH4 Relative Error (%)",
                fill = "Error Sign"
        ) +
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))


plot_NH3_err <- ggplot(avg_NH3_err, aes(x = lab.analyzer, y = mean_NH3_err, fill = bar_color)) +
        geom_bar(stat = "identity", width = 0.7) +
        geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.7) +
        facet_wrap(~ bg_direction) +
        scale_fill_manual(values = c("Positive" = "lightgreen", "Negative" = "red4")) +
        labs(
                title = "Mean Relative Errors of NH3 Concentrations      (Ref. = ATB_CRDS.P8)",
                x = "Laboratory",
                y = "NH3 Relative Error (%)",
                fill = "Error Sign"
        ) +
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Name the plot list
err_c_plots <- list(CO2 = plot_CO2_err, CH4 = plot_CH4_err, NH3 = plot_NH3_err)

# Save each plot as high-res PNG
for (name in names(err_c_plots)) {
        ggsave(
                filename = paste0(name, "_c_err.png"),
                plot = err_c_plots[[name]],
                width = 14,
                height = 8,
                dpi = 600
        )
}

# Combine saved PNGs into a single PDF
err_c_png_files <- paste0(names(err_c_plots), "_c_err.png")
err_c_img_list <- magick::image_read(err_c_png_files)
magick::image_write(image = err_c_img_list, path = "Ringversuche_concentration_err_plots.pdf", format = "pdf")

