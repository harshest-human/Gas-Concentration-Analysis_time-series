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

######## Development of functions #######
# Development of indirect.CO2.balance function
indirect.CO2.balance <- function(df) {
        library(dplyr)
        
        # ppm → mg/m³ at 0°C (273.15K) and 1 atm
        ppm_to_mgm3 <- function(ppm, molar_mass) {
                T_K <- 273.15      # Kelvin
                P   <- 101325      # Pa
                R   <- 8.314472    # J/mol/K
                (ppm * 1e-6) * molar_mass * 1e3 * P / (R * T_K)  # mg/m³
        }
        
        df %>%
                mutate(
                        # Hour of day
                        hour = as.numeric(format(DATE.TIME, "%H")),
                        
                        # Animal activity constants
                        a = 0.22,
                        h_min = 2.9,  # hour of minimal activity
                        
                        # Baseline heat production per cow (W)
                        phi = 5.6 * m_weight^0.75 + 22 * Y1_milk_prod + 1.6e-5 * p_pregnancy_day^3,
                        
                        # Temperature correction factor
                        t_factor = 1 + 4e-5 * (20 - temp_inside)^3,
                        phi_T_cor = phi * t_factor,
                        
                        # Relative animal activity correction
                        A_cor = 1 - a * 3 * sin((2*pi/24) * (hour + 6 - h_min)),
                        
                        # Heat production per cow corrected for T and activity
                        hpu_T_A_cor_per_cow = phi_T_cor * A_cor,
                        
                        # Total heat production for all animals
                        hpu_T_A_cor_all = hpu_T_A_cor_per_cow * n_dairycows,
                        
                        # Convert CO2, NH3, CH4 from ppm → mg/m³ (0°C)
                        CO2_mgm3_in = ppm_to_mgm3(CO2_ppm_in, 44.01),
                        CO2_mgm3_N  = ppm_to_mgm3(CO2_ppm_N,  44.01),
                        CO2_mgm3_S  = ppm_to_mgm3(CO2_ppm_S,  44.01),
                        
                        NH3_mgm3_in = ppm_to_mgm3(NH3_ppm_in, 17.031),
                        NH3_mgm3_N  = ppm_to_mgm3(NH3_ppm_N,  17.031),
                        NH3_mgm3_S  = ppm_to_mgm3(NH3_ppm_S,  17.031),
                        
                        CH4_mgm3_in = ppm_to_mgm3(CH4_ppm_in, 16.04),
                        CH4_mgm3_N  = ppm_to_mgm3(CH4_ppm_N,  16.04),
                        CH4_mgm3_S  = ppm_to_mgm3(CH4_ppm_S,  16.04),
                        
                        # Delta concentrations (inside - outside)
                        delta_CO2_N = CO2_mgm3_in - CO2_mgm3_N,
                        delta_CO2_S = CO2_mgm3_in - CO2_mgm3_S,
                        
                        delta_NH3_N = NH3_mgm3_in - NH3_mgm3_N,
                        delta_NH3_S = NH3_mgm3_in - NH3_mgm3_S,
                        
                        delta_CH4_N = CH4_mgm3_in - CH4_mgm3_N,
                        delta_CH4_S = CH4_mgm3_in - CH4_mgm3_S,
                        
                        # Ventilation rate (m³/h)
                        Q_vent_N = ifelse(delta_CO2_N != 0, hpu_T_A_cor_all / delta_CO2_N, NA_real_),
                        Q_vent_S = ifelse(delta_CO2_S != 0, hpu_T_A_cor_all / delta_CO2_S, NA_real_),
                        
                        # Instantaneous emissions (g/h) divided by 1000 to convert mg to g
                        e_NH3_g_h_N = delta_NH3_N * Q_vent_N / 1000, 
                        e_CH4_g_h_N = delta_CH4_N * Q_vent_N / 1000,
                        e_NH3_g_h_S = delta_NH3_S * Q_vent_S / 1000,
                        e_CH4_g_h_S = delta_CH4_S * Q_vent_S / 1000,
                        
                        # Annual emissions (kg/year) divided by 1000 to convert g to kg
                        e_NH3_kg_yr_N = e_NH3_g_h_N * 24 * 365 / 1000,
                        e_CH4_kg_yr_N = e_CH4_g_h_N * 24 * 365 / 1000,
                        e_NH3_kg_yr_S = e_NH3_g_h_S * 24 * 365 / 1000,
                        e_CH4_kg_yr_S = e_CH4_g_h_S * 24 * 365 / 1000
                )
}
        
# Development of ratio calculator function
ratcal <- function(df) {
        library(dplyr)
        
        gases_to_compare <- c("CH4", "NH3")
        locations <- c("in", "N", "S")
        
        for (gas in gases_to_compare) {
                for (loc in locations) {
                        col_gas <- paste0(gas, "_mgm3_", loc)
                        col_CO2 <- paste0("CO2_mgm3_", loc)
                        ratio_col <- paste0("r_", gas, "/CO2_", loc)  # only one location suffix
                        
                        df <- df %>%
                                mutate(
                                        !!ratio_col := ifelse(.data[[col_CO2]] != 0,
                                                              .data[[col_gas]] / .data[[col_CO2]] * 100,
                                                              NA_real_)
                                )
                }
        }
        
        return(df)
}

# Development of relative percentage error calculator function
relerrcal <- function(df, cols_to_calc, group_cols = c("DATE.TIME", "analyzer")) {
        library(dplyr)
        library(tidyr)
        
        # Pivot longer for calculations
        df_long <- df %>%
                select(all_of(c(group_cols, cols_to_calc))) %>%
                pivot_longer(
                        cols = all_of(cols_to_calc),
                        names_to = "variable",
                        values_to = "value"
                )
        
        # Calculate daily (or timestamp) mean per variable
        df_summary <- df_long %>%
                group_by(across(all_of(group_cols)), variable) %>%
                summarise(daily_mean = mean(value, na.rm = TRUE), .groups = "drop")
        
        # Baseline per group (excluding analyzer)
        df_baseline <- df_summary %>%
                group_by(across(setdiff(group_cols, "analyzer")), variable) %>%
                summarise(baseline = mean(daily_mean, na.rm = TRUE), .groups = "drop")
        
        # Calculate percent error
        df_errors <- df_summary %>%
                left_join(df_baseline, by = c(setdiff(group_cols, "analyzer"), "variable")) %>%
                mutate(error_pct = ((daily_mean - baseline) / baseline) * 100)
        
        # Pivot back to wide: daily_mean
        df_values_wide <- df_errors %>%
                select(all_of(group_cols), variable, daily_mean) %>%
                pivot_wider(
                        id_cols = all_of(group_cols),
                        names_from = variable,
                        values_from = daily_mean
                ) %>%
                mutate(across(where(is.numeric), ~ round(.x, 2)))
        
        # Pivot back to wide: error_pct
        df_errors_wide <- df_errors %>%
                select(all_of(group_cols), variable, error_pct) %>%
                mutate(variable = paste0("err_", variable)) %>%
                pivot_wider(
                        id_cols = all_of(group_cols),
                        names_from = variable,
                        values_from = error_pct
                ) %>%
                mutate(across(where(is.numeric), ~ round(.x, 2)))
        
        # Merge both wide tables
        df_final <- left_join(df_values_wide, df_errors_wide, by = group_cols)
        
        return(df_final)
}

# Development of function reshape parameters
reparam <- function(data) {
        library(dplyr)
        library(tidyr)
        library(stringr)
        
        meta_cols <- c("DATE.TIME", "analyzer")
        
        # Match all non-meta columns
        cols_to_use <- setdiff(names(data), meta_cols)
        
        df_long <- data %>%
                pivot_longer(
                        cols = all_of(cols_to_use),
                        names_to = "variable",
                        values_to = "value"
                ) %>%
                mutate(
                        value = round(value, 2),
                        # Extract location suffix
                        location_code = str_extract(variable, "(in|N|S)$"),
                        location = case_when(
                                location_code == "in" ~ "Ringline inside",
                                location_code == "N"  ~ "North background",
                                location_code == "S"  ~ "South background",
                                TRUE ~ NA_character_
                        ),
                        # Remove only the location suffix
                        var = str_remove(variable, "_(in|N|S)$"),
                        
                        # Assign var_type
                        var_type = case_when(
                                # Concentrations
                                str_detect(var, "^(CO2|CH4|NH3)_mgm3$") ~ "concentration measured",
                                # Deltas
                                str_detect(var, "^delta_(CO2|CH4|NH3)$") ~ "concentration delta",
                                # Ratios
                                str_detect(var, "^r_.+/.*") ~ "concentration ratio percentage",
                                # Relative errors
                                str_detect(var, "^err_(CO2|CH4|NH3)_mgm3$") ~ "concentration relative error",
                                str_detect(var, "^err_delta_(CO2|CH4|NH3)$") ~ "concentration delta relative error",
                                str_detect(var, "^err_r_.+/.*") ~ "concentration ratio percentage relative error",
                                # Ventilation
                                str_detect(var, "^Q_vent") ~ "ventilation rate",
                                str_detect(var, "^err_Q_vent") ~ "ventilation rate error",
                                # Emissions per hour
                                str_detect(var, "^e_(CO2|CH4|NH3)_mgh$") ~ "emission per hour",
                                str_detect(var, "^e_(CH4|NH3)_g_h$") ~ "emission gram per hour",
                                str_detect(var, "^err_e_(CO2|CH4|NH3)_mgh$") ~ "emission per hour relative error",
                                # Emissions per year
                                str_detect(var, "^e_(CO2|CH4|NH3)_kg_yr$") ~ "emission per year",
                                str_detect(var, "^e_(CH4|NH3)_kg_h$") ~ "emission kilogram per hour",
                                str_detect(var, "^err_e_(CO2|CH4|NH3)_kg_yr$") ~ "emission per year relative error",
                                str_detect(var, "^err_e_(CH4|NH3)_g_h$") ~ "emission gram per hour relative error",
                                str_detect(var, "^err_e_(CH4|NH3)_kg_h$") ~ "emission kilogram per hour relative error",
                                TRUE ~ "other"
                        )
                ) %>%
                select(all_of(meta_cols), location, var, var_type, value)
        
        return(df_long)
}

# Development of function stat_table
stat_table <- function(data, response_vars, group_vars = NULL, var_type_filter = NULL, time_group = c("hour", "day")) {
        require(dplyr)
        require(lubridate)
        require(DescTools)
        
        time_group <- match.arg(time_group)
        
        df <- data
        
        # Filter by var_type if provided
        if (!is.null(var_type_filter)) {
                df <- df %>% filter(var_type %in% var_type_filter)
        }
        
        # Create time grouping column
        df <- df %>%
                mutate(
                        time_grp = if (time_group == "hour") {
                                lubridate::hour(DATE.TIME)   # numeric hour
                        } else {
                                as.Date(DATE.TIME)           # date for day
                        }
                )
        
        # Ensure response_vars exist in 'var' column
        df <- df %>% filter(var %in% response_vars)
        
        # Group and summarise
        df %>%
                group_by(across(all_of(c("time_grp", group_vars))), var) %>%
                summarise(
                        n    = n(),
                        mean = round(mean(value, na.rm = TRUE), 2),
                        sd   = round(sd(value, na.rm = TRUE), 2),
                        cv   = round(DescTools::CoefVar(value, na.rm = TRUE) * 100, 2),
                        .groups = "drop"
                ) %>%
                rename(time = time_grp)
}

# Development of function HSD_table
HSD_table <- function(data, response_vars, group_var) {
        
        # Helper function to get labeled Tukey results for one variable
        get_tukey_labels <- function(var) {
                formula <- as.formula(paste(var, "~", group_var))
                tukey_res <- rstatix::tukey_hsd(data, formula)
                
                # Add label column with significance categories
                tukey_res <- tukey_res %>%
                        mutate(
                                label = case_when(
                                        is.na(p.adj)   ~ "-",
                                        p.adj <= 0.001 ~ "***",
                                        p.adj <= 0.01  ~ "**",
                                        p.adj <= 0.05  ~ "*",
                                        p.adj > 0.05   ~ "ns"
                                ),
                                color = case_when(
                                        label == "ns"  ~ "red")
                        ) %>%
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

# Development of function emission and concentration plot
emiconplot <- function(data, y = NULL, var_type_filter = NULL, x = "day") {
        library(dplyr)
        library(ggplot2)
        library(scales)
        library(stringr)
        
        # Standardize column names
        if ("var" %in% names(data) && !"variable" %in% names(data)) {
                data <- data %>% rename(variable = var)
        }
        
        # Add day and hour columns if DATE.TIME exists
        if ("DATE.TIME" %in% names(data)) {
                data <- data %>%
                        mutate(
                                day = as.factor(as.Date(DATE.TIME)),
                                hour = as.factor(format(DATE.TIME, "%H"))
                        )
        }
        
        required_cols <- c("location", "analyzer", "var_type", "variable", "value", x)
        missing_cols <- setdiff(required_cols, names(data))
        if (length(missing_cols) > 0) {
                stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))
        }
        
        # Filter if needed
        if (!is.null(var_type_filter)) {
                data <- data %>% filter(var_type %in% var_type_filter)
        }
        if (!is.null(y)) {
                data <- data %>% filter(variable %in% y)
        }
        
        # Summary stats
        summary_data <- data %>%
                group_by(across(all_of(c(x, "location", "analyzer", "var_type", "variable")))) %>%
                summarise(
                        mean_val = mean(value, na.rm = TRUE),
                        se_val = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
                        .groups = "drop"
                )
        
        # Clean facet labels
        summary_data <- summary_data %>%
                mutate(
                        variable_clean = variable %>%
                                str_remove("^e_") %>%           # remove e_ prefix
                                str_remove("^delta_") %>%       # remove delta_ prefix
                                str_remove("_mgm3$") %>%        # remove mgm3 suffix
                                str_remove("_g_h$") %>%         # remove g_h suffix
                                str_remove("_kg_yr$") %>%       # remove kg_yr suffix
                                str_remove("^r_") %>%           # remove r_ prefix
                                str_replace("^Q_vent.*", "Q")   # rename ventilation
                )
        
        # Y-axis labels by var_type
        ylab_map <- list(
                "concentration measured" = expression(paste("Mean ± SD (mg ", m^{-3}, ")")),
                "concentration delta"    = expression(paste(Delta, " concentration (mg ", m^{-3}, ")")),
                "concentration ratio percentage" = "Mean ± SD (%)",
                "ventilation rate"       = expression(paste("Ventilation rate (m"^{3}, " h"^{-1}, ")")),
                "emission gram per hour" = expression(paste("Emission (g h"^{-1}, ")")),
                "emission per year"      = expression(paste("Emission (kg yr"^{-1}, ")")),
                "concentration relative error" = "Relative error (%)"
        )
        
        categories_present <- unique(summary_data$var_type)
        if (length(categories_present) == 1 && categories_present %in% names(ylab_map)) {
                ylab_to_use <- ylab_map[[categories_present]]
        } else {
                ylab_to_use <- "Mean ± SD"
        }
        
        # X-axis label
        xlab_to_use <- switch(
                x,
                "hour" = "Hour (aggregated by days)",
                "day" = "Day (aggregated by hours)",
                x
        )
        
        # Colors & shapes
        analyzer_colors <- c(
                "FTIR.1" = "#1b9e77", "FTIR.2" = "#d95f02", "FTIR.3" = "#7570b3",
                "FTIR.4" = "#e7298a", "CRDS.1" = "#66a61e", "CRDS.2" = "#e6ab02",
                "CRDS.3" = "#a6761d"
        )
        analyzer_shapes <- c(
                "FTIR.1" = 0, "FTIR.2" = 1, "FTIR.3" = 2,
                "FTIR.4" = 5, "CRDS.1" = 15, "CRDS.2" = 19, "CRDS.3" = 17
        )
        
        # Force all factor levels on x-axis
        summary_data[[x]] <- factor(summary_data[[x]], levels = unique(data[[x]]))
        
        p <- ggplot(summary_data, aes_string(x = x, y = "mean_val", color = "analyzer", shape = "analyzer", group = "analyzer")) +
                geom_line() +
                geom_point(size = 2) +
                geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), width = 0.2) +
                scale_color_manual(values = analyzer_colors) +
                scale_shape_manual(values = analyzer_shapes) +
                facet_grid(variable_clean ~ location, scales = "free_y", switch = "y") +
                scale_y_continuous(breaks = pretty_breaks(n = 8)) +
                labs(x = xlab_to_use, y = ylab_to_use) +
                guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
                theme_classic() +
                theme(
                        text = element_text(size = 12),
                        axis.text = element_text(size = 12),
                        axis.title = element_text(size = 12),
                        strip.text = element_text(size = 12),
                        panel.border = element_rect(colour = "black", fill = NA),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        legend.position = "bottom",
                        legend.direction = "horizontal",
                        legend.box = "horizontal",
                        legend.box.just = "left",
                        legend.title = element_blank(),
                        legend.text = element_text(size = 12),
                        legend.key.width = unit(1, "lines")
                )
        
        print(p)
        return(p)
}

# Development of function HSD_boxplot
emiboxplot <- function(data, response_vars, group_var = "analyzer", var_type_filter = NULL) {
        library(dplyr)
        library(ggpubr)
        library(rstatix)
        library(ggplot2)
        library(scales)
        
        # Fixed facet variables
        facet_x <- "location"
        facet_y <- "variable"
        
        df <- data
        
        # Apply var_type filter if specified
        if (!is.null(var_type_filter)) {
                df <- df %>% filter(var_type %in% var_type_filter)
        }
        
        # Filter data to only requested variables
        data_sub <- df %>%
                filter(.data[[facet_y]] %in% response_vars)
        
        # Remove outliers per group_var and variable (omit outlier values)
        data_no_outliers <- data_sub %>%
                group_by(!!sym(group_var), .data[[facet_y]]) %>%
                filter(!is_outlier(value)) %>%
                ungroup()
        
        # y-axis labels by var_type from emiconplot style
        ylab_map <- list(
                concentration = expression(paste("Mean ± SD (mg ", m^{-3}, ")")),
                ratio = "Mean ± SD (%)",
                ventilation = expression(paste("Ventilation rate (m"^{-3} * " h"^{-1}, ")")),
                emission = expression(paste("Emission (g h"^{-1}, ")"))
        )
        
        # Try to infer var_type for ylab; if none, default
        categories_present <- unique(data_no_outliers$var_type)
        if (length(categories_present) == 1 && categories_present %in% names(ylab_map)) {
                ylab_to_use <- ylab_map[[categories_present]]
        } else {
                ylab_to_use <- "Value"
        }
        
        # Colors and shapes for analyzers from emiconplot
        analyzer_colors <- c(
                "FTIR.1" = "#1b9e77", "FTIR.2" = "#d95f02", "FTIR.3" = "#7570b3",
                "FTIR.4" = "#e7298a", "CRDS.1" = "#66a61e", "CRDS.2" = "#e6ab02",
                "CRDS.3" = "#a6761d"
        )
        analyzer_shapes <- c(
                "FTIR.1" = 0, "FTIR.2" = 1, "FTIR.3" = 2,
                "FTIR.4" = 5, "CRDS.1" = 15, "CRDS.2" = 19, "CRDS.3" = 17
        )
        
        # Ensure group_var is a factor sorted alphabetically for ascending order on x axis
        data_no_outliers[[group_var]] <- factor(data_no_outliers[[group_var]], levels = sort(unique(data_no_outliers[[group_var]])))
        
        p <- ggplot(data_no_outliers, aes_string(x = group_var, y = "value", fill = group_var)) +
                geom_boxplot(alpha = 0.8, fatten = 1.5, color = "black") + 
                stat_summary(fun = mean, geom = "point", aes(shape = !!sym(group_var)), 
                             size = 1.5, color = "white", fill = "white") +  
                facet_grid(reformulate(facet_x, facet_y), scales = "free_y", switch = "y") +
                scale_fill_manual(values = analyzer_colors) +
                scale_shape_manual(values = analyzer_shapes) +
                guides(fill = guide_legend(nrow = 1),
                       shape = guide_legend(nrow = 1)) +
                theme_classic() +
                theme(
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                        axis.text.y = element_text(size = 12),
                        axis.title = element_text(size = 12),
                        strip.text = element_text(size = 12),
                        legend.position = "none",  # hide legend completely
                        panel.border = element_rect(colour = "black", fill = NA),
                        axis.ticks.length = unit(0.2, "cm")
                ) +
                scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
                labs(y = ylab_to_use, x = group_var)
        
        print(p)
        return(p)
}

# Development of function emiheatmap
emiheatmap <- function(data, response_vars, group_var = "analyzer", var_type_filter = NULL) {
        library(dplyr)
        library(ggplot2)
        library(viridis)
        
        facet_x <- "location"
        facet_y <- "variable"
        
        df <- data
        
        # Apply var_type filter if specified
        if (!is.null(var_type_filter)) {
                df <- df %>% filter(var_type %in% var_type_filter)
        }
        
        # Filter and prepare data
        data_sub <- df %>%
                filter(.data[[facet_y]] %in% response_vars,
                       .data[[group_var]] != "FTIR.4")
        
        data_sub[[group_var]] <- factor(data_sub[[group_var]], levels = sort(unique(data_sub[[group_var]])))
        
        # Plot
        p <- ggplot(data_sub, aes(x = factor(hour), y = !!sym(group_var), fill = cv)) +
                geom_tile(color = "black") +
                facet_grid(reformulate(facet_x, facet_y), scales = "free_y", switch = "y") +
                scale_fill_viridis_c(
                        option = "plasma",
                        name = "CV (%)",
                        limits = c(0, 200),
                        breaks = seq(0, 200, 50),
                        oob = scales::squish,
                        guide = guide_colorbar(
                                barwidth = 10,
                                barheight = 0.6,
                                title.position = "top",
                                title.hjust = 0.5
                        )
                ) +
                labs(
                        x = "Hour of Day",
                        y = group_var
                ) +
                theme_minimal(base_size = 12) +
                theme(
                        panel.grid = element_blank(),
                        axis.text.x = element_text(size = 12),
                        axis.text.y = element_text(size = 12),
                        panel.border = element_rect(color = "black", fill = NA),
                        panel.spacing = unit(0.1, "lines"),
                        strip.background = element_rect(color = "black", fill = NA),
                        strip.text = element_text(size = 12),
                        legend.position = "bottom",
                        legend.title = element_text(size = 12),
                        legend.text = element_text(size = 12)
                )
        
        print(p)
        return(p)
}

######## Import Data #########
# Read all gas data
ATB_FTIR <- read.csv("20250408-15_ATB_wide_FTIR.1.csv")
LUFA_FTIR <- read.csv("20250408-15_LUFA_wide_FTIR.2.csv")
MBBM_FTIR <- read.csv("20250408-15_MBBM_wide_FTIR.4.csv")
ANECO_FTIR <- read.csv("20250408-15_ANECO_wide_FTIR.4.csv")
ATB_CRDS <- read.csv("20250408-15_ATB_wide_CRDS.1.csv")
UB_CRDS <- read.csv("20250408-15_UB_wide_CRDS.2.csv")
LUFA_CRDS <- read.csv("20250408-15_LUFA_wide_CRDS.3.csv")

# Time period
start_time <- "2025-04-08 12:00:00"
end_time   <- "2025-04-14 12:00:00"

# Combine all data set
gas_data <- bind_rows(LUFA_FTIR, ANECO_FTIR, MBBM_FTIR, ATB_FTIR, ATB_CRDS, LUFA_CRDS, UB_CRDS) %>%
        mutate(DATE.TIME = ymd_hms(DATE.TIME, tz = "UTC")) %>%  # specify correct tz here
        filter(DATE.TIME >= ymd_hms(start_time, tz = "UTC") &
                       DATE.TIME <= ymd_hms(end_time, tz = "UTC")) %>%
        select(DATE.TIME, analyzer, everything(), -lab) %>%
        arrange(DATE.TIME) %>% distinct() 

gas_data %>% count(DATE.TIME, analyzer, name = "n_obs") #to check how many observations


# Read animal data and fix time zone
animal_temp <- read.csv("20250408-15_LVAT_Animal_Temp_data.csv") %>%
        mutate(DATE.TIME = dmy_hm(DATE.TIME, tz = "UTC")) %>%
        filter(DATE.TIME >= ymd_hms(start_time, tz = "UTC") &
                       DATE.TIME <= ymd_hms(end_time, tz = "UTC")) %>%
        group_by(DATE.TIME) %>%
        summarise(across(everything(), ~ first(.x)), .groups = "drop")


# Join with gas data
input_combined <- left_join(gas_data, animal_temp, by = "DATE.TIME", relationship = "many-to-many") %>%
        arrange(DATE.TIME) %>% distinct() 


# Write csv
input_combined <- input_combined %>% mutate(across(where(is.numeric), ~ round(.x, 2))) 

write_excel_csv(input_combined, "20250408-15_ringversuche_input_combined_data.csv")

######## Computation of emissions, ventilation rates and ratios #########
# Calculate emissions and ratios using the function
emission_result <- indirect.CO2.balance(input_combined)

# Remove constants and forumal inputs
emission_result <- emission_result %>% select(-n_dairycows, -m_weight, -p_pregnancy_day, 
                                              -Y1_milk_prod, -temp_inside, -a, -h_min,
                                              -phi, -t_factor, -phi_T_cor, -A_cor,
                                              -hpu_T_A_cor_per_cow, -hpu_T_A_cor_all)

# Calculate gas ratios using function
emission_and_ratio <- ratcal(emission_result) 

# Write csv
write_excel_csv(emission_and_ratio, "20250408-15_ringversuche_emission_and_ratio.csv")

######## Computation of relative errors #########
# Store all relevant columns in a variable
vars <- c(
        # mg/m3 concentrations
        "CO2_mgm3_in", "CO2_mgm3_N", "CO2_mgm3_S",
        "NH3_mgm3_in", "NH3_mgm3_N", "NH3_mgm3_S",
        "CH4_mgm3_in", "CH4_mgm3_N", "CH4_mgm3_S",
        
        # Ratios
        "r_CH4/CO2_in", "r_CH4/CO2_N", "r_CH4/CO2_S",
        "r_NH3/CO2_in", "r_NH3/CO2_N", "r_NH3/CO2_S",
        
        # Deltas
        "delta_CO2_N", "delta_CO2_S",
        "delta_NH3_N", "delta_NH3_S",
        "delta_CH4_N", "delta_CH4_S",
        
        # Ventilation rates
        "Q_vent_N", "Q_vent_S",
        
        # Instantaneous emissions (g/h)
        "e_NH3_g_h_N", "e_CH4_g_h_N",
        "e_NH3_g_h_S", "e_CH4_g_h_S",
        
        # Annual emissions (kg/yr)
        "e_NH3_kg_yr_N", "e_CH4_kg_yr_N", 
        "e_NH3_kg_yr_S","e_CH4_kg_yr_S")

# Calculate relative errors using function
emission_ratio_error <- relerrcal(emission_and_ratio, vars)

# Write csv
write_excel_csv(emission_ratio_error, "20250408-15_ringversuche_emission_result_ratio_error_data.csv")

######## Reshape the data #########
emission_reshaped <-  reparam(emission_ratio_error)

write_excel_csv(emission_reshaped, "20250408-15_ringversuche_emission_reshaped.csv")

########### Calculate Statistical Summary (mean, SD and CV) #########
# Concentration and relative error
c_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3",
                          "err_CO2_mgm3", "err_CH4_mgm3", "err_NH3_mgm3"),
        var_type_filter = c("concentration measured", "concentration relative error"),
        group_vars = c("analyzer", "location"),
        time_group = "day")

# Delta concentrations and relative error
d_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("delta_CO2", "delta_CH4", "delta_NH3",
                          "err_delta_CO2", "err_delta_CH4", "err_delta_NH3"),
        var_type_filter = c("concentration delta", "concentration delta relative error"),
        group_vars = c("analyzer", "location"),
        time_group = "day")

# Ratio concentrations and relative error
r_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("r_CH4_N/CO2", "r_CH4_S/CO2", "r_CH4_in/CH4", "r_CH4_in/CO2",
                          "r_CO2_in/CO2", "r_NH3_N/CO2", "r_NH3_S/CO2", "r_NH3_in/CO2", "r_NH3_in/NH3",
                          "err_r_CH4_N/CO2", "err_r_CH4_S/CO2", "err_r_CH4_in/CH4", "err_r_CH4_in/CO2",
                          "err_r_CO2_in/CO2", "err_r_NH3_N/CO2", "err_r_NH3_S/CO2", "err_r_NH3_in/CO2", 
                          "err_r_NH3_in/NH3"),
        var_type_filter = c("concentration ratio percentage", "concentration ratio percentage relative error"),
        group_vars = c("analyzer", "location"),
        time_group = "day")

# Ventilation rates and error
q_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("Q_vent", "err_Q_vent"),
        var_type_filter = c("ventilation rate", "ventilation rate error"),
        group_vars = c("analyzer", "location"),
        time_group = "day")

# Emissions per hour / year and relative error
e_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("e_CH4_g_h", "e_CH4_kg_yr", "e_NH3_g_h", "e_NH3_kg_yr",
                          "err_e_CH4_g_h", "err_e_CH4_kg_yr", "err_e_NH3_g_h", "err_e_NH3_kg_yr"),
        var_type_filter = c("emission gram per hour", "emission per year",
                            "emission per hour relative error", "emission per year relative error"),
        group_vars = c("analyzer", "location"),
        time_group = "day")


# write csv
stats_list <- list(
        concentration = c_stat_sum,
        delta_concentration = d_stat_sum,
        ratio_concentration = r_stat_sum,
        ventilation = q_stat_sum,
        emission = e_stat_sum
)

for(name in names(stats_list)){
        write_excel_csv(stats_list[[name]], paste0("statistical_summary_", name, "_day.csv"))
}

######## ANOVA and HSD Summary ########
vars <- c(
        # mg/m3 concentrations
        "CO2_mgm3_in", "CO2_mgm3_N", "CO2_mgm3_S",
        "NH3_mgm3_in", "NH3_mgm3_N", "NH3_mgm3_S",
        "CH4_mgm3_in", "CH4_mgm3_N", "CH4_mgm3_S",
        
        # Deltas
        "delta_CO2_N", "delta_CO2_S",
        "delta_NH3_N", "delta_NH3_S",
        "delta_CH4_N", "delta_CH4_S",
        
        # Ventilation rates
        "Q_vent_N", "Q_vent_S",
        
        # Instantaneous emissions (g/h)
        "e_NH3_g_h_N", "e_CH4_g_h_N",
        "e_NH3_g_h_S", "e_CH4_g_h_S",
        
        # Annual emissions (kg/yr)
        "e_NH3_kg_yr_N", "e_CH4_kg_yr_N", 
        "e_NH3_kg_yr_S","e_CH4_kg_yr_S")

# Tukey HSD
# Remove rows where any response variable has NA/NaN/Inf
emission_clean <- emission_and_ratio %>%
        filter(if_all(all_of(vars), ~ !is.na(.) & is.finite(.)))

# Then perform Tukey HSD
result_HSD_summary <- HSD_table(data = emission_clean,
                                response_vars = vars,
                                group_var = "analyzer")

# Save the result
write_excel_csv(result_HSD_summary, "20250408_20250414_HSD_table.csv")


######## Daily Mean ± SD Trend Plots ########
# Concentrations only
c_plot <- emiconplot(
        data = emission_reshaped,
        x = "day",
        y = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3"),
        var_type_filter = "concentration measured")

# Ratios only
r_plot <- emiconplot(
        data = emission_reshaped,
        x = "day",
        y = c("r_CH4/CO2", "r_NH3/CO2"),
        var_type_filter = "concentration ratio percentage")

# Delta variables only
d_plot <- emiconplot(
        data = emission_reshaped,
        x = "day",
        y = c("delta_CO2", "delta_CH4", "delta_NH3"),
        var_type_filter = "concentration delta")

# Ventilation only
q_plot <- emiconplot(
        data = emission_reshaped %>% filter(analyzer != "FTIR.4"),
        x = "day",
        y = c("Q_vent"),
        var_type_filter = "ventilation rate")

# Emissions only
e_g_plot <- emiconplot(
        data = emission_reshaped %>% filter(analyzer != "FTIR.4"),
        x = "day",
        y = c("e_CH4_g_h", "e_NH3_g_h"),
        var_type_filter = "emission gram per hour")

e_kg_plot <- emiconplot(
        data = emission_reshaped %>% filter(analyzer != "FTIR.4"),
        x = "day",
        y = c("e_CH4_kg_yr", "e_NH3_kg_yr"),
        var_type_filter = "emission per year")

# Named list of plots
dailyplots <- list(
        c_plot = c_plot,
        r_plot = r_plot,
        d_plot = d_plot,
        q_plot = q_plot,
        e_g_plot = e_g_plot,
        e_kg_plot = e_kg_plot)

# Corresponding size settings (width, height)
plot_sizes <- list(
        c_plot = c(18, 8.5),
        r_plot = c(18, 6.5),
        d_plot = c(14, 6.5),
        q_plot = c(14, 4.4),
        e_g_plot = c(14, 6.5),
        e_kg_plot = c(14, 6.5))

# Save each plot using its specific size
for (plot_name in names(dailyplots)) {
        size <- plot_sizes[[plot_name]]
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = dailyplots[[plot_name]],
                width = size[1], height = size[2], dpi = 300)
}


########## Heat maps (CV) #############
q_heatmap <- emiheatmap(data = q_stat_sum, 
                        response_vars = c("Q_vent"))

e_heatmap <- emiheatmap(data = e_stat_sum, 
                        response_vars = c("e_CH4", "e_NH3"))

# Named list of plots
dailyplots <- list(
        c_heatmap = c_heatmap,
        r_heatmap = r_heatmap,
        q_heatmap = q_heatmap,
        e_heatmap = e_heatmap)

# coresponding size settings (width, height)
plot_sizes <- list(
        c_heatmap = c(18, 8.5),
        r_heatmap = c(18, 6.5),
        q_heatmap = c(14, 4.4),
        e_heatmap = c(14, 6.5))

# Save each plot using its specific size
for (plot_name in names(dailyplots)) {
        size <- plot_sizes[[plot_name]]
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = dailyplots[[plot_name]],
                width = size[1], height = size[2], dpi = 300)
}


########## Box plots (HSD) ventilation and emission rates ##############
e_boxplot <- emiboxplot(data = emission_reshaped,
                         response_vars = c("e_CH4", "e_NH3"),
                         group_var = "analyzer")

q_boxplot <- emiboxplot(data = emission_reshaped,
                         response_vars = c("Q_Vent_rate"),
                         group_var = "analyzer")

# Save e_boxplot
ggsave(filename = "e_boxplot.pdf",
       plot = e_boxplot,
       width = 13.5, height = 8.5, dpi = 300)

# Save q_boxplot
ggsave(filename = "q_boxplot.pdf",
       plot = q_boxplot,
       width = 13.5, height = 5.8, dpi = 300)



######### Linear mix modelling ##########
#Data processing
# Load animal, temperature, and wind data
animal_temp <- read.csv("20250408-15_LVAT_Animal_Temp_data.csv")
T_RH_HOBO <- read.csv("T_RH_08_04_2025_To_30_06_2025.csv")
USA_Mst <- read.csv("USA_Mst_5_min_2025_04_08.csv")
USA_Trv <- read.csv("USA_Trv_5_min_2025_04_08.csv")

T_RH_HOBO <- T_RH_HOBO %>%
        mutate(date = trimws(date),
               time = trimws(time))

env_anm_data <- T_RH_HOBO %>%
        mutate(
                # Combine date and time as character string
                datetime_str = paste(date, time),
                # Parse using strptime with explicit format
                DATE.TIME = as.POSIXct(strptime(datetime_str, format = "%d/%m/%y %H:%M:%S")),
                # Floor to hour
                DATE.TIME = floor_date(DATE.TIME, "hour")
        ) %>%
        group_by(DATE.TIME) %>%
        summarise(
                T_outside = mean(T_outside, na.rm = TRUE),
                RH_out = mean(RH_out, na.rm = TRUE),
                T_inside = mean(T_inside, na.rm = TRUE),
                RH_inside = mean(RH_inside, na.rm = TRUE),
                .groups = "drop"
        )

library(lme4)

# Fit the model
q_model <- lmer(value ~ 1 + (1 | analyzer), data = emission_reshaped %>% filter(variable == "Q_Vent_rate"))

# Extract variance components
var_comp <- as.data.frame(VarCorr(q_model))

# Between-analyzer variance
sigma_R2 <- var_comp$vcov[var_comp$grp == "analyzer"]

# Within-analyzer variance (residual)
sigma_r2 <- var_comp$vcov[var_comp$grp == "Residual"]

# Calculate ICC
icc <- sigma_R2 / (sigma_R2 + sigma_r2)

cat("Between-analyzer variance (sigma_R^2):", sigma_R2, "\n")
cat("Within-analyzer variance (sigma_r^2):", sigma_r2, "\n")
cat("Intraclass correlation coefficient (ICC):", round(icc, 3), "\n")
