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
source("D:/Data Analysis/Picarro-G2508_CRDS_gas_measurement/remove_outliers_function.R")

######## Functions #########
indirect.CO2.balance <- function(
                df,
                CO2_ppm_in = NULL, CO2_ppm_out = NULL,
                NH3_ppm_in = NULL, NH3_ppm_out = NULL,
                CH4_ppm_in = NULL, CH4_ppm_out = NULL
) {
        library(dplyr)
        
        # ppm → mg/m³ at 0°C (273.15K) and 1 atm
        ppm_to_mgm3 <- function(ppm, molar_mass) {
                T_K <- 273.15      # Kelvin
                P   <- 101325      # Pa
                R   <- 8.314472    # J/mol/K
                (ppm * 1e-6) * molar_mass * 1e3 * P / (R * T_K)
        }
        
        # ---- SAFETY CHECKS ----
        if (is.null(CO2_ppm_in) | is.null(CO2_ppm_out))
                stop("Please provide both CO2_ppm_in and CO2_ppm_out column names.")
        
        if (is.null(NH3_ppm_in) | is.null(NH3_ppm_out))
                stop("Please provide both NH3_ppm_in and NH3_ppm_out column names.")
        
        if (is.null(CH4_ppm_in) | is.null(CH4_ppm_out))
                stop("Please provide both CH4_ppm_in and CH4_ppm_out column names.")
        
        
        df <- df %>%
                mutate(
                        # Hour of day
                        hour = as.numeric(format(DATE.TIME, "%H")),
                        
                        # Animal activity constants
                        a = 0.22,
                        h_min = 2.9,
                        
                        # Heat production and correction factors
                        phi = 5.6 * m_weight^0.75 + 22 * Y1_milk_prod + 1.6e-5 * p_pregnancy_day^3,
                        t_factor = 1 + 4e-5 * (20 - temp_in)^3,
                        phi_T_cor = phi * t_factor,
                        A_cor = 1 - a * (sin((2*pi/24) * (hour + 6 - h_min))),
                        hpu_T_A_cor_per_cow = phi_T_cor * A_cor,
                        
                        # CO2 production (mg/s → mg/h)
                        PCO2 = (0.185 * hpu_T_A_cor_per_cow) * 1000,
                        
                        # Convert ppm to mg/m³
                        CO2_mgm3_in  = ppm_to_mgm3(.data[[CO2_ppm_in]], 44.01),
                        CO2_mgm3_out = ppm_to_mgm3(.data[[CO2_ppm_out]], 44.01),
                        
                        NH3_mgm3_in  = ppm_to_mgm3(.data[[NH3_ppm_in]], 17.031),
                        NH3_mgm3_out = ppm_to_mgm3(.data[[NH3_ppm_out]], 17.031),
                        
                        CH4_mgm3_in  = ppm_to_mgm3(.data[[CH4_ppm_in]], 16.04),
                        CH4_mgm3_out = ppm_to_mgm3(.data[[CH4_ppm_out]], 16.04),
                        
                        # Deltas
                        delta_CO2 = CO2_mgm3_in - CO2_mgm3_out,
                        delta_NH3 = NH3_mgm3_in - NH3_mgm3_out,
                        delta_CH4 = CH4_mgm3_in - CH4_mgm3_out,
                        
                        # Ventilation rate
                        Q_vent = ifelse(delta_CO2 != 0, PCO2 / delta_CO2, NA_real_),
                        
                        # Emissions (g/h)
                        e_NH3_gh = (delta_NH3 * Q_vent / 1000) * n_dairycows_in,
                        e_CH4_gh = (delta_CH4 * Q_vent / 1000) * n_dairycows_in,
                        
                        # Emissions per livestock unit (kg/year)
                        e_NH3_ghLU = (e_NH3_gh * 500) / (n_dairycows_in * m_weight),
                        e_CH4_ghLU = (e_CH4_gh * 500) / (n_dairycows_in * m_weight)
                )
        
        return(df)
}

###### Import Data ########
# Animal Parameters Data
animal_temp_data <- read.csv("20250928-20250930_LVAT_Animal_Temp_data.csv") %>%
        mutate(DATE.TIME = dmy_hm(DATE.TIME),
               DATE.TIME = floor_date(DATE.TIME, unit = "hour"))

animal_daily_avg <- animal_temp_data %>%
        group_by(DATE.TIME) %>%
        summarise(temp_in = mean(temp_inside, na.rm = TRUE),
                  n_dairycows_in = mean(n_dairycows, na.rm = TRUE),
                  m_weight = mean(m_weight, na.rm = TRUE),
                  p_pregnancy_day = mean(p_pregnancy_day, na.rm = TRUE),
                  Y1_milk_prod = mean(Y1_milk_prod, na.rm = TRUE),
                  .groups = "drop") %>%
        arrange(DATE.TIME)

# Gas Concentration Data
gas_data <- read.csv("20251010_high_resolution_gas_concentration_data.csv") %>%
        mutate(DATE.TIME = ymd_hms(DATE.TIME)) 

gas_avg <- gas_data %>%
        filter(DATE.TIME >= ymd_hms("2025-08-28 13:00:00"))%>%
        mutate(DATE.TIME = floor_date(DATE.TIME, unit = "hour")) %>%
        group_by(DATE.TIME, location, analyzer) %>%
        summarise(CO2_ppm = mean(CO2, na.rm = TRUE),
                  CH4_ppm = mean(CH4, na.rm = TRUE),
                  NH3_ppm = mean(NH3, na.rm = TRUE),
                  .groups = "drop") %>%
        pivot_wider(names_from = c(location, analyzer),
                    values_from = c(CO2_ppm, CH4_ppm, NH3_ppm),
                    names_sep = "_") %>%
        arrange(DATE.TIME) %>%
        remove_outliers(exclude_cols = c("DATE.TIME"))

# Combined data
input_combined <- gas_avg %>%
        left_join(animal_daily_avg, by = "DATE.TIME") %>%
        arrange(DATE.TIME) %>%
        distinct()

# Define in/out locations for each analyzer
in_locs_CRDS9 <- c(1,3,4,6,7,9,10,12,13,15,16,18,22,24,53)
in_locs_CRDS8 <- c(28,30,31,33,34,36,37,39,40,42,43,45,46,48,49,51)
out_loc_CRDS9 <- 52
# Fixed out columns from CRDS9
out_cols_CRDS9 <- list(
        CO2 = "CO2_ppm_52_CRDS9",
        NH3 = "NH3_ppm_52_CRDS9",
        CH4 = "CH4_ppm_52_CRDS9"
)

# Define a helper function to run emission for a given analyzer
run_emission <- function(input_df, analyzer, in_locs, out_cols_fixed = NULL) {
        results_list <- list()
        
        for (loc in in_locs) {
                # Build in columns for this analyzer
                CO2_in_col <- paste0("CO2_ppm_", loc, "_", analyzer)
                NH3_in_col <- paste0("NH3_ppm_", loc, "_", analyzer)
                CH4_in_col <- paste0("CH4_ppm_", loc, "_", analyzer)
                
                # Use fixed out columns if provided
                if (!is.null(out_cols_fixed)) {
                        CO2_out_col <- out_cols_fixed$CO2
                        NH3_out_col <- out_cols_fixed$NH3
                        CH4_out_col <- out_cols_fixed$CH4
                } else {
                        stop("No out column specified!")
                }
                
                # Run the indirect.CO2.balance function
                results_list[[paste0(analyzer, "_in_", loc)]] <- indirect.CO2.balance(
                        input_df,
                        CO2_ppm_in  = CO2_in_col,  CO2_ppm_out  = CO2_out_col,
                        NH3_ppm_in  = NH3_in_col,  NH3_ppm_out  = NH3_out_col,
                        CH4_ppm_in  = CH4_in_col,  CH4_ppm_out  = CH4_out_col
                )
        }
        
        # Combine results into a single dataframe
        combined <- bind_rows(results_list, .id = "in_location")
        combined$analyzer <- analyzer
        return(combined)
}

# Run emissions for each analyzer
emission_CRDS9 <- run_emission(input_combined, "CRDS9",
                               in_locs_CRDS9, out_cols_fixed = out_cols_CRDS9) %>%
        remove_outliers(exclude_cols = c("DATE.TIME"))

emission_CRDS8 <- run_emission(input_combined, "CRDS8",
                               in_locs_CRDS8, out_cols_fixed = out_cols_CRDS9) %>%
        remove_outliers(exclude_cols = c("DATE.TIME"))

# Combine both results safely
emission_result <- bind_rows(emission_CRDS9, emission_CRDS8) %>%
        # Extract numeric location from in_location (e.g. "CRDS9_in_18" → 18)
        mutate(location = as.numeric(gsub(".*_in_(\\d+)", "\\1", in_location))) %>%
        # Arrange columns in desired order
        select(
                DATE.TIME, location,
                delta_CO2, delta_NH3, delta_CH4,
                Q_vent,
                e_NH3_gh, e_NH3_ghLU,
                e_CH4_gh, e_CH4_ghLU
        ) %>%
        arrange(DATE.TIME, location)


############# DATA Analysis ##########
gas_avg_loc <- gas_data %>%
        filter(DATE.TIME >= ymd_hms("2025-08-28 13:00:00"),
               !location %in% c("52")) %>%
        group_by(location, vgroup) %>%
        summarise(CO2_ppm = mean(CO2, na.rm = TRUE),
                  CH4_ppm = mean(CH4, na.rm = TRUE),
                  NH3_ppm = mean(NH3, na.rm = TRUE),
                  NHCO    = mean(NHCO, na.rm = TRUE),
                  CHCO    = mean(CHCO, na.rm = TRUE),
                  .groups = "drop") %>%       
        mutate(location = as.character(location)) %>%
        bind_rows(tibble(location = "baseline",
                        CO2_ppm = mean(.$CO2_ppm, na.rm = TRUE),
                        CH4_ppm = mean(.$CH4_ppm, na.rm = TRUE),
                        NH3_ppm = mean(.$NH3_ppm, na.rm = TRUE),
                        NHCO    = mean(.$NHCO, na.rm = TRUE),
                        CHCO    = mean(.$CHCO, na.rm = TRUE),
                        vgroup = "baseline"))

# Define numeric columns for PRE
numeric_cols <- c("CO2_ppm", "CH4_ppm", "NH3_ppm", "NHCO", "CHCO")

# Separate baseline and other locations
baseline_row <- gas_avg_loc %>% 
        filter(location == "baseline") %>%
        mutate(across(all_of(numeric_cols), ~ 0, .names = "{.col}_PRE"))

other_rows <- gas_avg_loc %>% filter(location != "baseline")

# Calculate PRE vs baseline
other_rows <- other_rows %>%
        mutate(across(
                all_of(numeric_cols),
                ~ ((.x - baseline_row[[cur_column()]]) / baseline_row[[cur_column()]]) * 100,
                .names = "{.col}_PRE"
        ))

# Combine baseline and other rows
gas_avg_loc <- bind_rows(other_rows, baseline_row) %>%
        mutate(location_num = suppressWarnings(as.numeric(location))) %>%  # convert numeric locations
        arrange(location_num) %>%                                         # sort numeric locations first
        select(-location_num) %>%
        mutate(location = as.character(location)) %>%
        arrange(as.numeric(location)) %>% 
        mutate(location = factor(location, levels = unique(location)))                                             # remove helper column


# -------------------------------
# 1. Compute emission averages per location
# -------------------------------
emission_avg_loc <- emission_result %>%
        filter(!location %in% c("10", "48")) %>%
        group_by(location) %>%
        summarise(
                delta_CO2  = mean(delta_CO2, na.rm = TRUE),
                
                delta_NH3  = mean(delta_NH3, na.rm = TRUE),
                
                delta_CH4  = mean(delta_CH4, na.rm = TRUE),
                
                Q_vent     = mean(Q_vent, na.rm = TRUE),
                
                e_CH4_ghLU = mean(e_CH4_ghLU, na.rm = TRUE),
                
                e_NH3_ghLU = mean(e_NH3_ghLU, na.rm = TRUE),
                
                .groups = "drop")


# -------------------------------
# 2. Assign vertical groups
# -------------------------------
vertical_groups <- list(
        top       = c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49),
        mid       = c(2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50),
        bottom    = c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51),
        ref.ring  = 53,
        out       = 52
)

emission_avg_loc <- emission_avg_loc %>%
        mutate(location = as.character(location),
               vgroup = case_when(
                       location %in% vertical_groups$top    ~ "top",
                       location %in% vertical_groups$mid    ~ "mid",
                       location %in% vertical_groups$bottom ~ "bottom",
                       location == as.character(vertical_groups$`ref.ring`) ~ "ref.ring",
                       location == as.character(vertical_groups$out)        ~ "out"
               ))

# -------------------------------
# 3. Compute baseline (mean across locations)
# -------------------------------
baseline_row <- emission_avg_loc %>%
        summarise(
                location     = "baseline",
                delta_CO2    = mean(delta_CO2, na.rm = TRUE),
                delta_NH3    = mean(delta_NH3, na.rm = TRUE),
                delta_CH4    = mean(delta_CH4, na.rm = TRUE),
                Q_vent       = mean(Q_vent, na.rm = TRUE),
                e_CH4_ghLU   = mean(e_CH4_ghLU, na.rm = TRUE),
                e_NH3_ghLU   = mean(e_NH3_ghLU, na.rm = TRUE),
                vgroup       = "baseline"
        )

# -------------------------------
# 4. Separate baseline and other rows
# -------------------------------
numeric_cols <- c("delta_CO2", "delta_NH3", "delta_CH4", "Q_vent", "e_CH4_ghLU", "e_NH3_ghLU")

other_rows <- emission_avg_loc %>% filter(location != "baseline")

# Calculate PRE relative to baseline
other_rows <- other_rows %>%
        mutate(across(
                all_of(numeric_cols),
                ~ ((.x - baseline_row[[cur_column()]]) / baseline_row[[cur_column()]]) * 100,
                .names = "{.col}_PRE"
        ))

# Baseline PRE columns are zero
baseline_row <- baseline_row %>%
        mutate(across(all_of(numeric_cols), ~ 0, .names = "{.col}_PRE"))

# -------------------------------
# 5. Combine and arrange sequentially
# -------------------------------
emission_avg_loc <- bind_rows(other_rows, baseline_row) %>%
        mutate(location_num = suppressWarnings(as.numeric(location))) %>%
        # put baseline at the end
        mutate(location_num = ifelse(is.na(location_num), max(location_num, na.rm = TRUE) + 1, location_num)) %>%
        arrange(location_num) %>%
        select(-location_num) %>%
        mutate(location = factor(location, levels = c(sort(as.numeric(other_rows$location)), "baseline")))

############# DATA Visualization ##########
# Define fill colors
point_fill <- c("top" = "orange", "mid" = "green3", "bottom" = "steelblue1", "ref.ring" = "purple")

# Define y-axis limits and breaks
# y-axis limits and breaks
y_limits <- c(-60, 60)
y_breaks <- seq(-60, 60, by = 10)

# Directory to save plots
save_dir <- "plots"
if(!dir.exists(save_dir)) dir.create(save_dir)

#-----------------------------
# Data sources for each variable
#-----------------------------
CO2_ppm_PRE_data    <- gas_avg_loc     %>% filter(!location %in% c("10","48","baseline")) %>% arrange(location) %>% mutate(location = factor(location, levels = location))
NH3_ppm_PRE_data    <- gas_avg_loc     %>% filter(!location %in% c("10","48","baseline")) %>% arrange(location) %>% mutate(location = factor(location, levels = location))
CH4_ppm_PRE_data    <- gas_avg_loc     %>% filter(!location %in% c("10","48","baseline")) %>% arrange(location) %>% mutate(location = factor(location, levels = location))
Q_vent_PRE_data     <- emission_avg_loc %>% filter(!location %in% c("10","48","baseline")) %>% arrange(location) %>% mutate(location = factor(location, levels = location))
e_CH4_ghLU_PRE_data <- emission_avg_loc %>% filter(!location %in% c("10","48","baseline")) %>% arrange(location) %>% mutate(location = factor(location, levels = location))
e_NH3_ghLU_PRE_data <- emission_avg_loc %>% filter(!location %in% c("10","48","baseline")) %>% arrange(location) %>% mutate(location = factor(location, levels = location))

#-----------------------------
# Helper: Lollipop plot generator
#-----------------------------
make_lollipop <- function(data, var, label) {
        ggplot(data, aes(x = location, y = .data[[var]], color = vgroup)) +
                geom_segment(aes(xend = location, y = 0, yend = .data[[var]]), linewidth = 0.8) +
                geom_point(size = 4) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                scale_color_manual(values = point_fill) +
                scale_y_continuous(limits = y_limits, breaks = y_breaks) +
                theme_minimal(base_size = 13) +
                theme(
                        axis.text.x = element_text(angle = 0, vjust = 0.5),
                        legend.position = "none"
                ) +
                labs(x = "Sampling point", y = label)
}

#-----------------------------
# Create individual lollipop plots
#-----------------------------
p_CO2_ppm_PRE    <- make_lollipop(CO2_ppm_PRE_data,    "CO2_ppm_PRE",    "CO2 Relative Error (%)")
p_NH3_ppm_PRE    <- make_lollipop(NH3_ppm_PRE_data,    "NH3_ppm_PRE",    "NH3 Relative Error (%)")
p_CH4_ppm_PRE    <- make_lollipop(CH4_ppm_PRE_data,    "CH4_ppm_PRE",    "CH4 Relative Error (%)")
p_Q_vent_PRE     <- make_lollipop(Q_vent_PRE_data,     "Q_vent_PRE",     "Q_vent Relative Error (%)")
p_e_CH4_ghLU_PRE <- make_lollipop(e_CH4_ghLU_PRE_data, "e_CH4_ghLU_PRE", "CH4 Emission Relative Error (%)")
p_e_NH3_ghLU_PRE <- make_lollipop(e_NH3_ghLU_PRE_data, "e_NH3_ghLU_PRE", "NH3 Emission Relative Error (%)")

#-----------------------------
# Save all plots individually
#-----------------------------
#ggsave(paste0(save_dir, "/CO2_ppm_PRE2024.png"),    p_CO2_ppm_PRE,    width = 12, height = 4)
#ggsave(paste0(save_dir, "/NH3_ppm_PRE2024.png"),    p_NH3_ppm_PRE,    width = 12, height = 4)
#ggsave(paste0(save_dir, "/CH4_ppm_PRE2024.png"),    p_CH4_ppm_PRE,    width = 12, height = 4)
ggsave(paste0(save_dir, "/CO2_ppm_PRE.png"),    p_CO2_ppm_PRE,    width = 12, height = 4)
ggsave(paste0(save_dir, "/NH3_ppm_PRE.png"),    p_NH3_ppm_PRE,    width = 12, height = 4)
ggsave(paste0(save_dir, "/CH4_ppm_PRE.png"),    p_CH4_ppm_PRE,    width = 12, height = 4)
ggsave(paste0(save_dir, "/Q_vent_PRE.png"),     p_Q_vent_PRE,     width = 12, height = 4)
ggsave(paste0(save_dir, "/e_CH4_ghLU_PRE.png"), p_e_CH4_ghLU_PRE, width = 12, height = 4)
ggsave(paste0(save_dir, "/e_NH3_ghLU_PRE.png"), p_e_NH3_ghLU_PRE, width = 12, height = 4)

############ CV Plots ##########
emission_cv <- emission_result %>%
        filter(!location %in% c("10", "48")) %>%
        group_by(location) %>%
        summarise(
                delta_CO2_cv    = (sd(delta_CO2, na.rm = TRUE) / mean(delta_CO2, na.rm = TRUE)) * 100,
                delta_NH3_cv    = (sd(delta_NH3, na.rm = TRUE) / mean(delta_NH3, na.rm = TRUE)) * 100,
                delta_CH4_cv    = (sd(delta_CH4, na.rm = TRUE) / mean(delta_CH4, na.rm = TRUE)) * 100,
                Q_vent_cv       = (sd(Q_vent, na.rm = TRUE) / mean(Q_vent, na.rm = TRUE)) * 100,
                e_CH4_ghLU_cv   = (sd(e_CH4_ghLU, na.rm = TRUE) / mean(e_CH4_ghLU, na.rm = TRUE)) * 100,
                e_NH3_ghLU_cv   = (sd(e_NH3_ghLU, na.rm = TRUE) / mean(e_NH3_ghLU, na.rm = TRUE)) * 100,
                .groups = "drop"
        ) %>%
        mutate(location = as.numeric(location)) %>%
        arrange(location)


ggplot(emission_cv, aes(x = factor(location), y = delta_CO2_cv)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = round(delta_CO2_cv, 1)), vjust = -0.5, size = 3) +
        theme_minimal(base_size = 12) +
        labs(x = "Sampling point", y = "CV (%)") +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5))
ggsave("cv_delta_CO2.png", width = 6, height = 4, dpi = 300)

ggplot(emission_cv, aes(x = factor(location), y = delta_CH4_cv)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = round(delta_CH4_cv, 1)), vjust = -0.5, size = 3) +
        theme_minimal(base_size = 12) +
        labs(x = "Sampling point", y = "CV (%)") +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5))
ggsave("cv_delta_CH4.png", width = 6, height = 4, dpi = 300)

ggplot(emission_cv, aes(x = factor(location), y = delta_NH3_cv)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = round(delta_NH3_cv, 1)), vjust = -0.5, size = 3) +
        theme_minimal(base_size = 12) +
        labs(x = "Sampling point", y = "CV (%)") +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5))
ggsave("cv_delta_NH3.png", width = 6, height = 4, dpi = 300)

ggplot(emission_cv, aes(x = factor(location), y = Q_vent_cv)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = round(Q_vent_cv, 1)), vjust = -0.5, size = 3) +
        theme_minimal(base_size = 12) +
        labs(x = "Sampling point", y = "CV (%)") +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5))
ggsave("cv_Q_vent.png", width = 6, height = 4, dpi = 300)

ggplot(emission_cv, aes(x = factor(location), y = e_CH4_ghLU_cv)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = round(e_CH4_ghLU_cv, 1)), vjust = -0.5, size = 3) +
        theme_minimal(base_size = 12) +
        labs(x = "Sampling point", y = "CV (%)") +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5))
ggsave("cv_e_CH4_ghLU.png", width = 6, height = 4, dpi = 300)

ggplot(emission_cv, aes(x = factor(location), y = e_NH3_ghLU_cv)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = round(e_NH3_ghLU_cv, 1)), vjust = -0.5, size = 3) +
        theme_minimal(base_size = 12) +
        labs(x = "Sampling point", y = "CV (%)") +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5))
ggsave("cv_e_NH3_ghLU.png", width = 6, height = 4, dpi = 300)



 