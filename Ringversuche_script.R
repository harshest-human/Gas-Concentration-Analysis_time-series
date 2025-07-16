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


######## Development of functions #######


# Development of indirect.CO2.balance function
indirect.CO2.balance <- function(df) {
        df %>%
                mutate(
                        P_CO2_term = 0.185,
                        a = 0.22,
                        h_min = 2.9,
                        CO2_Molmass = 44.01,
                        NH3_Molmass = 17.031,
                        CH4_Molmass = 16.04,
                        R_gas_constant = 8.314472,
                        p_pressure_ref = 1013,
                        hour = hour(DATE.TIME),
                        
                        A_corr = 1 - a * 3 * sin((2 * pi / 24) * (hour + 6 - h_min)),
                        Phi_tot = 5.6 * (m_weight ^ 0.75) + 22 * Y1_milk_prod + 1.6e-5 * (p_pregnancy_day ^ 3),
                        Phi_T_corr = Phi_tot * (1 + 4e-5 * (20 - Temperature)^3),
                        hpu_T_A_corr_all_animal = ((Phi_T_corr / 1000) * A_corr) * n_animals,
                        n_LU = (n_animals * m_weight) / 500,
                        P_CO2_T_A_all_animal = hpu_T_A_corr_all_animal * P_CO2_term,
                        
                        # Ventilation rate North
                        Q_Vent_rate_N = P_CO2_T_A_all_animal / ((CO2_in - CO2_N) * 1e-6),
                        # Ventilation rate South
                        Q_Vent_rate_S = P_CO2_T_A_all_animal / ((CO2_in - CO2_S) * 1e-6),
                        
                        # Delta in ppm (original concentration differences)
                        delta_NH3_N_ppm = NH3_in - NH3_N,
                        delta_CH4_N_ppm = CH4_in - CH4_N,
                        delta_CO2_N_ppm = CO2_in - CO2_N,
                        
                        delta_NH3_S_ppm = NH3_in - NH3_S,
                        delta_CH4_S_ppm = CH4_in - CH4_S,
                        delta_CO2_S_ppm = CO2_in - CO2_S,
                        
                        # Delta converted to mg/m3
                        delta_NH3_N_mgm3 = (0.1 * NH3_Molmass * p_pressure_ref * delta_NH3_N_ppm) / ((Temperature + 273.15) * R_gas_constant),
                        delta_CH4_N_mgm3 = (0.1 * CH4_Molmass * p_pressure_ref * delta_CH4_N_ppm) / ((Temperature + 273.15) * R_gas_constant),
                        delta_CO2_N_mgm3 = (0.1 * CO2_Molmass * p_pressure_ref * delta_CO2_N_ppm) / ((Temperature + 273.15) * R_gas_constant),
                        
                        delta_NH3_S_mgm3 = (0.1 * NH3_Molmass * p_pressure_ref * delta_NH3_S_ppm) / ((Temperature + 273.15) * R_gas_constant),
                        delta_CH4_S_mgm3 = (0.1 * CH4_Molmass * p_pressure_ref * delta_CH4_S_ppm) / ((Temperature + 273.15) * R_gas_constant),
                        delta_CO2_S_mgm3 = (0.1 * CO2_Molmass * p_pressure_ref * delta_CO2_S_ppm) / ((Temperature + 273.15) * R_gas_constant),
                        
                        # Emissions North
                        e_NH3_N = (delta_NH3_N_mgm3 * Q_Vent_rate_N) / 1000,
                        e_CH4_N = (delta_CH4_N_mgm3 * Q_Vent_rate_N) / 1000,
                        e_CO2_N = (delta_CO2_N_mgm3 * Q_Vent_rate_N) / 1000,
                        
                        # Emissions South
                        e_NH3_S = (delta_NH3_S_mgm3 * Q_Vent_rate_S) / 1000,
                        e_CH4_S = (delta_CH4_S_mgm3 * Q_Vent_rate_S) / 1000,
                        e_CO2_S = (delta_CO2_S_mgm3 * Q_Vent_rate_S) / 1000,
                        )
}

# Development of function rm_outliers_IQR
rm_outliers_IQR <- function(df, cols) {
        for (col in cols) {
                Q1 <- quantile(df[[col]], 0.25, na.rm = TRUE)
                Q3 <- quantile(df[[col]], 0.75, na.rm = TRUE)
                IQR_val <- Q3 - Q1
                lower <- Q1 - 1.5 * IQR_val
                upper <- Q3 + 1.5 * IQR_val
                
                # Replace outliers with NA (or any flag value)
                df[[col]] <- ifelse(df[[col]] < lower | df[[col]] > upper, NA, df[[col]])
        }
        return(df)
}


# Development of function stat_table
stat_table <- function(data, response_vars, group_var) {
        # Load required libraries
        require(dplyr)
        require(DescTools)
        require(tidyselect)
        
        data %>%
                group_by(.data[[group_var]]) %>%
                summarise(
                        n = n(),
                        across(
                                all_of(response_vars),
                                list(
                                        mean = ~mean(., na.rm = TRUE),
                                        sd   = ~sd(., na.rm = TRUE),
                                        cv   = ~DescTools::CoefVar(., na.rm = TRUE) * 100
                                ),
                                .names = "{.fn}_{.col}"
                        ),
                        .groups = "drop"
                ) %>%
                rename_with(~paste0(.x, " (ppm)"), .cols = starts_with(c("mean_", "sd_"))) %>%
                rename_with(~paste0(.x, " (%)"),       .cols = starts_with("cv_"))
}


# Development of function HSD_table
HSD_table <- function(data, response_vars, group_var) {
        
        # Helper function to get labeled Tukey results for one variable
        get_tukey_labels <- function(var) {
                formula <- as.formula(paste(var, "~", group_var))
                tukey_res <- rstatix::tukey_hsd(data, formula)
                
                # Add label column with significance categories
                tukey_res <- tukey_res %>%
                        mutate(label = case_when(
                                is.na(p.adj) ~ "-",
                                p.adj <= 0.001 ~ "≤ 0.001",
                                p.adj <= 0.01 ~ "≤ 0.01",
                                p.adj <= 0.05 ~ "≤ 0.05",
                                TRUE ~ "> 0.05"
                        )) %>%
                        select(group1, group2, label) %>%
                        rename(!!var := label)
                
                return(tukey_res)
        }
        
        # Get results for all variables as a list of data frames
        results_list <- lapply(response_vars, get_tukey_labels)
        
        # Join all by group1 and group2 (pair of analyzers)
        combined <- Reduce(function(x, y) full_join(x, y, by = c("group1", "group2")), results_list)
        
        # Arrange and rename grouping columns
        combined <- combined %>%
                arrange(group1, group2) %>%
                rename(analyzer.1 = group1, analyzer.2 = group2)
        
        return(combined)
}


# Development of function relerror
relerror <- function(data, gases, reference) {
        library(dplyr)
        library(tidyr)
        library(stringr)
        
        vars_base <- c(gases, "Q_Vent_rate")
        
        build_patterns <- function(g) {
                c(
                        paste0("^", g, "(_in|_N|_S)?$"),
                        paste0("^delta_", g, "(_N|_S)(_ppm|_mgm3)?$"),
                        paste0("^e_", g, "(_N|_S)?$")
                )
        }
        
        all_patterns <- unlist(lapply(vars_base, build_patterns))
        
        selected_vars <- unique(unlist(
                lapply(all_patterns, function(p) grep(p, names(data), value = TRUE))
        ))
        
        needed_cols <- c("DATE.TIME", "day", "hour", "analyzer", selected_vars)
        
        # Average values by grouping metadata and analyzer (still keeping analyzer in grouping)
        err <- data %>%
                select(all_of(needed_cols)) %>%
                group_by(DATE.TIME, day, hour, analyzer) %>%
                summarise(across(all_of(selected_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
        
        # Now compute errors without pivoting (keep analyzer column separate)
        
        # Get reference rows (for chosen analyzer)
        ref_data <- err %>% filter(analyzer == reference) %>% select(-analyzer)
        
        # Join err data with ref_data by DATE.TIME, day, hour
        err_joined <- err %>%
                filter(analyzer != reference) %>%
                left_join(ref_data, by = c("DATE.TIME", "day", "hour"), suffix = c("", "_ref"))
        
        # For each selected_var compute relative error compared to *_ref
        for (var in selected_vars) {
                err_joined[[paste0(var, "_err")]] <- 100 * (err_joined[[var]] - err_joined[[paste0(var, "_ref")]]) / err_joined[[paste0(var, "_ref")]]
        }
        
        # Keep original columns plus error columns
        error_cols <- paste0(selected_vars, "_err")
        
        result <- err_joined %>%
                select(DATE.TIME, day, hour, analyzer, all_of(selected_vars), all_of(error_cols))
        
        return(result)
}


######## Import Gas Data #########
# Load animal and temperature data
animal_temp <- read.csv("20250408-15_LVAT_Animal_Temp_data.csv")

# Read and convert DATE.TIME for FTIR data
ATB_FTIR <- read.csv("20250408-15_long_ATB_FTIR.1.csv")
LUFA_FTIR <- read.csv("20250408-15_long_LUFA_FTIR.2.csv")
ANECO_FTIR <- read.csv("20250408-15_long_ANECO_FTIR.3.csv")
MBBM_FTIR <- read.csv("20250408-15_long_MBBM_FTIR.4.csv")

# Read and convert DATE.TIME for CRDS data
ATB_CRDS <- read.csv("20250408-15_long_ATB_CRDS.1.csv")
UB_CRDS <- read.csv("20250408-15_long_UB_CRDS.2.csv")
LUFA_CRDS <- read.csv("20250408-15_long_LUFA_CRDS.3.csv")

# Combine all data set
gas_data <- bind_rows(LUFA_FTIR, ANECO_FTIR, MBBM_FTIR, ATB_FTIR, ATB_CRDS, LUFA_CRDS, UB_CRDS)
input_combined <-full_join(gas_data, animal_temp, by = "DATE.TIME") %>%
        mutate(day = as.Date(DATE.TIME)) %>% 
        select(DATE.TIME, day, hour, everything())


# Write csv
input_combined <- input_combined %>% select(DATE.TIME, hour, everything())
write.csv(input_combined, "20250408-15_ringversuche_input_combined_data.csv", row.names = FALSE)

######## Computation of ventilation rates and emissions #########
# Convert DATE.TIME format
input_combined <- input_combined %>% filter(DATE.TIME >= "2025-04-08 12:00:00" & DATE.TIME <= "2025-04-14 10:10:00")

# Calculate emissions using the function
emission_combined  <- indirect.CO2.balance(input_combined)

# Write csv
write.csv(emission_combined, "20250408-15_ringversuche_emission_combined_data.csv", row.names = FALSE)

# Reload the final result
# DATE.TIME conversion
result <- emission_combined %>%
        mutate(
                DATE.TIME = as.POSIXct(DATE.TIME, format = "%Y-%m-%d %H:%M:%S"),
                analyzer = as.factor(analyzer),
                day = as.factor(day)
        )


######## Statistical Data Analysis #######
# Define Variables
vars <- c(
        "CO2_in", "CH4_in", "NH3_in",
        "CO2_N", "CH4_N", "NH3_N",
        "CO2_S", "CH4_S", "NH3_S",
        "Q_Vent_rate_N", "Q_Vent_rate_S",
        "e_CO2_N", "e_CH4_N", "e_NH3_N",
        "e_CO2_S", "e_CH4_S", "e_NH3_S",
        "delta_CO2_N_mgm3", "delta_CH4_N_mgm3", "delta_NH3_N_mgm3",
        "delta_CO2_S_mgm3", "delta_CH4_S_mgm3", "delta_NH3_S_mgm3"
)

vars_base <- c(
        "CO2", "CH4", "NH3",
        "Q_Vent_rate",
        "e_CO2", "e_CH4", "e_NH3",
        "delta_CO2", "delta_CH4", "delta_NH3"
)

# Descriptive statistics
result_stat_summary <- stat_table(data = result, response_vars = vars, group_var = "analyzer")
write_excel_csv(result_stat_summary, "20250408_20250414_stat_table.csv")

# Tukey HSD
result_HSD_summary <- HSD_table(data = result, response_vars = vars, group_var = "analyzer")
write_excel_csv(result_HSD_summary, "20250408_20250414_HSD_table.csv")

# Relative error
result_err_summary <- relerror(data =result, gases = c("CO2", "CH4", "NH3"), reference = "CRDS.1")
write_excel_csv(result_err_summary, "20250408_20250414_relerr_table.csv")


######## Data Visualization ########
# Plot and save all variables from `vars` (from `result`)
for (y_var in vars) {
        y_sym <- sym(y_var)
        plots <- emicon.plot(df = result, !!y_sym)
        
        for (day_name in names(plots)) {
                plot_obj <- plots[[day_name]]
                
                # Show the plot in the viewer
                print(plot_obj)
                
                # Save the plot
                ggsave(
                        filename = paste0(y_var, "_", day_name, ".png"),
                        plot = plot_obj,
                        width = 8, height = 5, dpi = 300
                )
                cat("Saved result plot:", paste0(y_var, "_", day_name, ".png"), "\n")
        }
}

# Plot and save error variables from `result_err_summary` / `err_long`
error_vars <- c(
        paste0("e_", c("CO2", "CH4", "NH3"), "_err"),
        paste0("delta_", c("CO2", "CH4", "NH3"), "_err")
)

for (y_var in error_vars) {
        y_sym <- sym(y_var)
        plots <- emicon.plot(df = err_long, !!y_sym)
        
        for (day_name in names(plots)) {
                plot_obj <- plots[[day_name]]
                
                # Show the plot
                print(plot_obj)
                
                # Save the plot
                ggsave(
                        filename = paste0(y_var, "_", day_name, ".png"),
                        plot = plot_obj,
                        width = 8, height = 5, dpi = 300
                )
                cat("Saved error plot:", paste0(y_var, "_", day_name, ".png"), "\n")
        }
}
