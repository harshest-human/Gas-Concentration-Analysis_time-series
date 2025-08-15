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
# Function to remove outliers based on IQR
remove_outliers <- function(df) {
        # Convert all columns except DATE.TIME to numeric
        df[ , -which(names(df) == "DATE.TIME")] <- lapply(df[ , -which(names(df) == "DATE.TIME")], function(x) as.numeric(as.character(x)))
        
        # Identify numeric columns
        num_cols <- sapply(df, is.numeric)
        
        # Replace outliers with NA
        df[num_cols] <- lapply(df[num_cols], function(x) {
                Q1 <- quantile(x, 0.25, na.rm = TRUE)
                Q3 <- quantile(x, 0.75, na.rm = TRUE)
                IQR <- Q3 - Q1
                x[x < (Q1 - 1.5*IQR) | x > (Q3 + 1.5*IQR)] <- NA
                return(x)
        })
        
        return(df)
}

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
        
        df <- df %>%
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
        
        # Calculate CH4/CO2 and NH3/CO2 ratios (in %)
        gases_to_compare <- c("CH4", "NH3")
        locations <- c("in", "N", "S")
        
        for (gas in gases_to_compare) {
                for (loc in locations) {
                        col_gas <- paste0(gas, "_mgm3_", loc)
                        col_CO2 <- paste0("CO2_mgm3_", loc)
                        ratio_col <- paste0("r_", gas, "/CO2_", loc)
                        
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
        
# Development of ratio and relative error calculator function
relerror <- function(df) {
        library(dplyr)
        library(tidyr)
        library(stringr)
        
        meta_cols <- c("DATE.TIME", "analyzer")
        
        # Identify numeric measurement columns
        measure_cols <- df %>%
                select(-all_of(meta_cols)) %>%
                select(where(is.numeric)) %>%
                names()
        
        # Step 1: Calculate baseline per DATE.TIME
        baseline <- df %>%
                group_by(DATE.TIME) %>%
                summarise(across(all_of(measure_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
                mutate(analyzer = "baseline") %>%
                select(names(df))
        
        # Step 2: Combine baseline with original data
        combined <- bind_rows(df, baseline)
        
        # Step 3: Calculate relative errors
        combined <- combined %>%
                group_by(DATE.TIME) %>%
                mutate(across(all_of(measure_cols),
                              .fns = list(err = ~ 100 * (.x - mean(.x[analyzer == "baseline"], na.rm = TRUE)) /
                                                  mean(.x[analyzer == "baseline"], na.rm = TRUE)),
                              .names = "err_{.col}")) %>%
                ungroup()
        
        # Step 4: Pivot to long format
        err_cols <- paste0("err_", measure_cols)
        
        df_long <- combined %>%
                pivot_longer(
                        cols = all_of(measure_cols),
                        names_to = "var",
                        values_to = "value"
                ) %>%
                pivot_longer(
                        cols = all_of(err_cols),
                        names_to = "err_var",
                        values_to = "err"
                ) %>%
                filter(str_remove(err_var, "^err_") == var) %>%
                select(-err_var) %>%
                mutate(
                        location = case_when(
                                str_detect(var, "_in$") ~ "Ringline inside",
                                str_detect(var, "_N$")  ~ "North background",
                                str_detect(var, "_S$")  ~ "South background",
                                TRUE ~ NA_character_
                        ),
                        var = str_remove(var, "_(in|N|S)$")
                ) %>%
                select(DATE.TIME, location, analyzer, var, value, err) %>%
                arrange(DATE.TIME, var, analyzer, location)
        
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
        
        # Keep only desired variables
        df <- df %>% filter(var %in% response_vars)
        
        # Create time grouping column
        df <- df %>%
                mutate(
                        time_grp = if (time_group == "hour") {
                                lubridate::hour(DATE.TIME)
                        } else {
                                as.Date(DATE.TIME)
                        }
                )
        
        # Group and summarise
        summary_df <- df %>%
                group_by(across(all_of(c("time_grp", group_vars))), var) %>%
                summarise(
                        n    = n(),
                        mean = round(mean(value, na.rm = TRUE), 2),
                        sd   = round(sd(value, na.rm = TRUE), 2),
                        cv   = round(DescTools::CoefVar(value, na.rm = TRUE) * 100, 2),
                        min  = round(min(value, na.rm = TRUE), 2),
                        max  = round(max(value, na.rm = TRUE), 2),
                        .groups = "drop"
                )
        
        # Rename the time column based on time_group
        time_col_name <- if (time_group == "hour") "hour" else "day"
        summary_df <- summary_df %>% rename(!!time_col_name := time_grp)
        
        return(summary_df)
}

# Development of function HSD_table
HSD_table <- function(data, group_var = "analyzer") {
        
        # Helper function to get labeled Tukey results for one variable
        get_tukey_labels <- function(var_name) {
                df_var <- data %>% filter(var == var_name)
                formula <- as.formula(paste("value ~", group_var))
                tukey_res <- rstatix::tukey_hsd(df_var, formula)
                
                # Add significance label
                tukey_res <- tukey_res %>%
                        mutate(
                                label = case_when(
                                        is.na(p.adj)    ~ "-",
                                        p.adj <= 0.001  ~ "***",
                                        p.adj <= 0.01   ~ "**",
                                        p.adj <= 0.05   ~ "*",
                                        p.adj > 0.05    ~ "ns"
                                )
                        ) %>%
                        select(group1, group2, label) %>%
                        rename(!!var_name := label)
                
                return(tukey_res)
        }
        
        # Unique variables
        vars <- unique(data$var)
        
        # Apply Tukey for all variables
        results_list <- lapply(vars, get_tukey_labels)
        
        # Join all by group1 and group2
        combined <- Reduce(function(x, y) full_join(x, y, by = c("group1", "group2")), results_list)
        
        # Arrange and rename grouping columns
        combined <- combined %>%
                arrange(group1, group2) %>%
                rename(analyzer.1 = group1, analyzer.2 = group2)
        
        return(combined)
}

# Development of function emission and concentration plot
emiconplot <- function(data, y = NULL, location_filter = NULL, plot_err = FALSE) {
        library(dplyr)
        library(ggplot2)
        library(scales)
        library(stringr)
        
        # Standardize column names
        if ("var" %in% names(data) && !"variable" %in% names(data)) {
                data <- data %>% rename(variable = var)
        }
        
        # Filter by location if specified
        if (!is.null(location_filter)) {
                data <- data %>% filter(location %in% location_filter)
        }
        
        # Filter by variable if specified
        if (!is.null(y)) {
                data <- data %>% filter(variable %in% y)
        }
        
        # Choose value column
        value_col <- if (plot_err) "err" else "value"
        
        # Summary statistics by DATE.TIME, analyzer, variable
        summary_data <- data %>%
                group_by(DATE.TIME, analyzer, location, variable) %>%
                summarise(
                        mean_val = mean(.data[[value_col]], na.rm = TRUE),
                        sd_val = sd(.data[[value_col]], na.rm = TRUE),
                        .groups = "drop"
                )
        
        # Smart facet labels (expressions)
        facet_labels_expr <- c(
                "CO2_mgm3"     = "c[CO2]~'(mg '*m^-3*')'",
                "CH4_mgm3"     = "c[CH4]~'(mg '*m^-3*')'",
                "NH3_mgm3"     = "c[NH3]~'(mg '*m^-3*')'",
                "r_CH4/CO2"    = "c[CH4]/c[CO2]~'('*'%'*')'",
                "r_NH3/CO2"    = "c[NH3]/c[CO2]~'('*'%'*')'",
                "Q_vent"       = "'Q Ventilation rate'~' ('*m^3*' h^-1)'",
                "e_CH4_g_h"    = "e[CH4]~'(g '*h^-1*')'",
                "e_NH3_g_h"    = "e[NH3]~'(g '*h^-1*')'",
                "e_CH4_kg_yr"  = "e[CH4]~'(kg '*yr^-1*')'",
                "e_NH3_kg_yr"  = "e[NH3]~'(kg '*yr^-1*')'",
                "delta_CO2"    = "Delta*c[CO2]~'(mg '*m^-3*')'",
                "delta_CH4"    = "Delta*c[CH4]~'(mg '*m^-3*')'",
                "delta_NH3"    = "Delta*c[NH3]~'(mg '*m^-3*')'"
        )
        
        # Add facet_label column as factor to preserve order
        summary_data <- summary_data %>%
                mutate(facet_label = factor(
                        facet_labels_expr[as.character(variable)],
                        levels = facet_labels_expr[y]
                ))
        
        # Colors & shapes for analyzers
        analyzer_colors <- c(
                "FTIR.1" = "#1b9e77", "FTIR.2" = "#d95f02", "FTIR.3" = "#7570b3",
                "FTIR.4" = "#e7298a", "CRDS.1" = "#66a61e", "CRDS.2" = "#e6ab02",
                "CRDS.3" = "#a6761d", "baseline" = "black"
        )
        analyzer_shapes <- c(
                "FTIR.1" = 0, "FTIR.2" = 1, "FTIR.3" = 2,
                "FTIR.4" = 5, "CRDS.1" = 15, "CRDS.2" = 19,
                "CRDS.3" = 17, "baseline" = 4
        )
        
        # X-axis breaks
        summary_data$DATE.TIME <- as.POSIXct(summary_data$DATE.TIME)
        x_breaks <- seq(
                from = min(summary_data$DATE.TIME, na.rm = TRUE),
                to = max(summary_data$DATE.TIME, na.rm = TRUE),
                by = "6 hours"
        )
        
        # Y-axis label
        ylab_to_use <- if (plot_err) "Relative error (%)" else "Mean ± SD"
        
        # Plot
        p <- ggplot(summary_data, aes(x = DATE.TIME, y = mean_val,
                                      color = analyzer, shape = analyzer, group = analyzer)) +
                geom_line(size = 0.5, alpha = 0.6) +
                geom_point(size = 2) +
                geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                              width = 0.2) +
                facet_grid(facet_label ~ ., scales = "free_y", switch = "y",
                           labeller = label_parsed) +
                scale_color_manual(values = analyzer_colors) +
                scale_shape_manual(values = analyzer_shapes) +
                scale_y_continuous(breaks = pretty_breaks(n = 9)) +
                scale_x_datetime(breaks = x_breaks, date_labels = "%Y-%m-%d %H:%M") +
                labs(x = "DATE.TIME", y = ylab_to_use, title = unique(summary_data$location)) +
                theme_classic() +
                theme(
                        text = element_text(size = 10),
                        axis.text = element_text(size = 10),
                        axis.title = element_text(size = 10),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                        strip.text.y.left = element_text(size = 10, vjust = 0.5),
                        panel.border = element_rect(color = "black", fill = NA),
                        legend.position = "bottom",
                        legend.title = element_blank(),
                        plot.title = element_text(hjust = 0.5)
                ) +
                guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
        
        print(p)
        return(p)
}

# Development of function HSD_boxplot
emiboxplot <- function(data, y = NULL, group_var = "analyzer", var_type_filter = NULL, x = "day") {
        library(dplyr)
        library(ggplot2)
        library(rstatix)
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
        
        # Filter by var_type if specified
        if (!is.null(var_type_filter)) {
                data <- data %>% filter(var_type %in% var_type_filter)
        }
        
        # Filter by requested y variables
        if (!is.null(y)) {
                data <- data %>% filter(variable %in% y)
        }
        
        # Remove outliers per group_var and variable
        df <- data %>%
                group_by(.data[[group_var]], variable) %>%
                filter(!is_outlier(value)) %>%
                ungroup()
        
        # Clean variable names for facets and preserve order
        df <- df %>%
                mutate(
                        variable_clean = variable %>%
                                str_remove("^e_") %>%
                                str_remove("^delta_") %>%
                                str_remove("_g_h$") %>%
                                str_remove("_kg_yr$") %>%
                                str_remove("_mgm3$") %>%
                                str_replace("^Q_vent.*", "Q"),
                        variable_clean = factor(variable_clean, levels = unique(variable_clean))
                )
        
        # Dynamic y-axis labels by var_type
        ylab_map <- list(
                "concentration mgm3" = expression(paste("Concentration (mg ", m^{-3}, ")")),
                "concentration ppm" = expression(paste("Concentration (ppm)")),
                "concentration delta" = expression(paste(Delta, "Concentration")),
                "concentration ratio percentage" = "Concentration ratio (%)",
                "ventilation rate" = expression(paste("Ventilation rate (m"^{3}, " h"^{-1}, ")")),
                "emission g per hour" = expression(paste("Emission (g h"^{-1}, ")")),
                "emission kg per year" = expression(paste("Emission (kg yr"^{-1}, ")")),
                "concentration mgm3 relative error" = "Relative error (%)",
                "concentration ppm relative error" = "Relative error (%)",
                "delta relative error" = "Relative error (%)",
                "ratio percentage relative error" = "Relative error (%)",
                "ventilation rate error" = "Relative error (%)",
                "emission g per hour relative error" = "Relative error (%)",
                "emission kg per year relative error" = "Relative error (%)"
        )
        categories_present <- unique(df$var_type)
        if (length(categories_present) == 1 && categories_present %in% names(ylab_map)) {
                ylab_to_use <- ylab_map[[categories_present]]
        } else {
                ylab_to_use <- "Value"
        }
        
        # Dynamic x-axis label
        xlab_to_use <- switch(
                x,
                "hour" = "Hour (aggregated by days)",
                "day" = "Day (aggregated by hours)",
                x
        )
        
        # Colors
        analyzer_colors <- c(
                "FTIR.1" = "#1b9e77", "FTIR.2" = "#d95f02", "FTIR.3" = "#7570b3",
                "FTIR.4" = "#e7298a", "CRDS.1" = "#66a61e", "CRDS.2" = "#e6ab02",
                "CRDS.3" = "#a6761d"
        )
        
        # Force x-axis factor
        if (x == "hour") {
                df[[x]] <- factor(df[[x]], levels = sprintf("%02d", 0:23))
        } else {
                df[[x]] <- factor(df[[x]], levels = unique(df[[x]]))
        }
        
        # Compute boxplot stats manually after outlier removal
        box_stats <- df %>%
                group_by(.data[[x]], variable_clean, .data[[group_var]], location) %>%
                summarise(
                        ymin = min(value),
                        lower = quantile(value, 0.25),
                        middle = median(value),
                        upper = quantile(value, 0.75),
                        ymax = max(value),
                        .groups = "drop"
                )
        
        # Plot using stat="identity" to show whiskers as min/max after outlier removal
        p <- ggplot() +
                geom_boxplot(
                        data = box_stats,
                        aes(
                                x = .data[[x]],
                                ymin = ymin, lower = lower, middle = middle,
                                upper = upper, ymax = ymax,
                                fill = .data[[group_var]]
                        ),
                        stat = "identity",
                        color = "black",
                        alpha = 0.8,
                        outlier.shape = NA
                ) +
                facet_grid(variable_clean ~ location, scales = "free_y", switch = "y") +
                scale_fill_manual(values = analyzer_colors) +
                guides(fill = guide_legend(nrow = 1)) +
                theme_classic() +
                theme(
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                        axis.text.y = element_text(size = 12),
                        axis.title = element_text(size = 12),
                        strip.text = element_text(size = 12),
                        legend.position = "bottom",
                        panel.border = element_rect(colour = "black", fill = NA),
                        axis.ticks.length = unit(0.2, "cm")
                ) +
                scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
                labs(y = ylab_to_use, x = xlab_to_use)
        
        print(p)
        return(p)
}

# Development of function emiheatmap
emiheatmap <- function(stat_df, y = NULL, group_var = "analyzer", time_group = c("hour", "day")) {
        library(dplyr)
        library(ggplot2)
        
        time_group <- match.arg(time_group)
        
        df <- stat_df
        
        # Filter by y variables if provided
        if (!is.null(y)) {
                df <- df %>% filter(var %in% y)
        }
        
        # Determine the column to use for x-axis
        x_col <- if (time_group == "hour") {
                if(!"hour" %in% colnames(df)) stop("Column 'hour' not found in stat_df")
                "hour"
        } else {
                if(!"time" %in% colnames(df)) stop("Column 'time' not found in stat_df for daily plot")
                "time"
        }
        
        # Convert hour to "HH:MM" format
        if (time_group == "hour") {
                df$hour_label <- sprintf("%02d:00", df$hour)
                x_col <- "hour_label"
        }
        
        # Ensure group_var is a factor
        df[[group_var]] <- factor(df[[group_var]], levels = sort(unique(df[[group_var]])))
        
        # Ensure var is a factor for faceting
        df$var <- factor(df$var, levels = unique(df$var))
        
        # Fixed color scale 0-200% with custom gradient
        legend_breaks <- seq(0, 150, by = 10)
        custom_colors <- c(
                "darkgreen","green4", "green2","yellow", "orange", "red", "darkred")
        
        # Plot
        p <- ggplot(df, aes_string(x = x_col, y = group_var, fill = "cv")) +
                geom_tile(color = "white") + 
                facet_grid(location ~ var, scales = "free_y") +
                scale_fill_gradientn(
                        colors = custom_colors,
                        limits = c(0, 150),
                        breaks = seq(0, 150, by = 30),
                        name = "CV (%)"
                ) +
                labs(
                        x = ifelse(time_group == "hour", "Hour of Day", "Date"),
                        y = NULL,
                        title = NULL
                ) +
                theme_minimal(base_size = 12) +
                theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                      axis.text.y = element_text(size = 10),
                      strip.text.x = element_blank(),
                      strip.text = element_text(size = 12),
                      panel.grid = element_blank(),
                      legend.position = "bottom",
                      legend.key.width = unit(2, "cm"),
                      plot.title = element_blank())
        
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

# Apply to each dataset
ATB_FTIR_clean <- remove_outliers(ATB_FTIR)
LUFA_FTIR_clean <- remove_outliers(LUFA_FTIR)
MBBM_FTIR_clean <- remove_outliers(MBBM_FTIR)
ANECO_FTIR_clean <- remove_outliers(ANECO_FTIR)
ATB_CRDS_clean <- remove_outliers(ATB_CRDS)
UB_CRDS_clean <- remove_outliers(UB_CRDS)
LUFA_CRDS_clean <- remove_outliers(LUFA_CRDS)

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
emission_result <- emission_result %>% select(-hour, -n_dairycows, -m_weight, -p_pregnancy_day, 
                                              -Y1_milk_prod, -temp_inside, -a, -h_min,
                                              -phi, -t_factor, -phi_T_cor, -A_cor,
                                              -hpu_T_A_cor_per_cow, -hpu_T_A_cor_all)%>%
        mutate(across(where(is.numeric), ~ round(.x, 2))) 

# Write csv 
write_excel_csv(emission_result, "20250408-15_ringversuche_emission_result.csv")

######## Reshape the data #########
emission_reshaped <-  relerror(emission_result) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2))) 

# Write csv
write_excel_csv(emission_reshaped, "20250408-15_ringversuche_emission_reshaped.csv")

######## ANOVA and HSD Summary ########
result_HSD_summary <- HSD_table(data = emission_reshaped)

write_excel_csv(result_HSD_summary, "result_HSD_summary.csv")

######## Normailty check ########
# Filter the gases of interest
gas_subset <- emission_reshaped %>%
        filter(var %in% c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3"))

# Histogram + density plot
ggplot(gas_subset, aes(x = value)) +
        geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
        geom_density(color = "red", size = 1) +
        facet_grid(var ~ analyzer, scales = "free") +
        theme_minimal() +
        labs(title = "Histogram and Density of Gas Measurements", x = "Value", y = "Density")

# Q-Q plot
ggplot(gas_subset, aes(sample = value)) +
        stat_qq() +
        stat_qq_line(color = "red") +
        facet_grid(var ~ analyzer, scales = "free") +
        theme_minimal() +
        labs(title = "Q-Q Plot of Gas Measurements")

 ######## Daily Mean ± SD Trend Plots ########
# Concentrations only
c_r_in_mgm3_plot <- emiconplot(
        data = emission_reshaped,
        y = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3", 
              "r_CH4/CO2", "r_NH3/CO2"),
        location_filter = "Ringline inside")

c_r_N_mgm3_plot <- emiconplot(
        data = emission_reshaped,
        y = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3", 
              "r_CH4/CO2", "r_NH3/CO2"),
        location_filter = "North background")

c_r_S_mgm3_plot <- emiconplot(
        data = emission_reshaped,
        y = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3", 
              "r_CH4/CO2", "r_NH3/CO2"),
        location_filter = "South background")


# Delta variables only
delta_N_plot <- emiconplot(
        data = emission_reshaped,
        y = c("delta_CO2", "delta_CH4", "delta_NH3"),
        location_filter = "North background")

delta_S_plot <- emiconplot(
        data = emission_reshaped,
        y = c("delta_CO2", "delta_CH4", "delta_NH3"),
        location_filter = "South background")

# Ventilation only
q_vent_N_plot <- emiconplot(
        data = emission_reshaped %>%
                filter(var == "Q_vent") %>%
                filter(value >= quantile(value, 0.02, na.rm = TRUE) &
                               value <= quantile(value, 0.98, na.rm = TRUE)),
        y = c("Q_vent"),
        location_filter = "North background")

q_vent_S_plot <- emiconplot(
        data = emission_reshaped %>%
                filter(var == "Q_vent") %>%
                filter(value >= quantile(value, 0.02, na.rm = TRUE) &
                                value <= quantile(value, 0.98, na.rm = TRUE)),
        y = c("Q_vent"),
        location_filter = "South background")

# Emissions only
emission_gph_N_plot <- emiconplot(
        data = bind_rows(
                # Filter CH4 to 1–99% quantiles
                emission_reshaped %>%
                        filter(var == "e_CH4_g_h") %>%
                        filter(value >= quantile(value, 0.01, na.rm = TRUE) &
                                       value <= quantile(value, 0.99, na.rm = TRUE)),
                
                # Filter NH3 to 0.1–99.9% quantiles
                emission_reshaped %>%
                        filter(var == "e_NH3_g_h") %>%
                        filter(value >= quantile(value, 0.001, na.rm = TRUE) &
                                       value <= quantile(value, 0.999, na.rm = TRUE))),
        y = c("e_CH4_g_h", "e_NH3_g_h"),
        location_filter = "North background")


emission_gph_S_plot <- emiconplot(
        data = bind_rows(
                # Filter CH4 to 1–99% quantiles
                emission_reshaped %>%
                        filter(var == "e_CH4_g_h") %>%
                        filter(value >= quantile(value, 0.01, na.rm = TRUE) &
                                       value <= quantile(value, 0.99, na.rm = TRUE)),
                
                # Filter NH3 to 0.1–99.9% quantiles
                emission_reshaped %>%
                        filter(var == "e_NH3_g_h") %>%
                        filter(value >= quantile(value, 0.001, na.rm = TRUE) &
                                       value <= quantile(value, 0.999, na.rm = TRUE))),
        y = c("e_CH4_g_h", "e_NH3_g_h"),
        location_filter = "South background")


# Named list of plots (matching new names)
dailyplots <- list(
        c_mgm3_plot = c_mgm3_plot,
        c_ppm_plot   = c_ppm_plot,
        ratio_plot   = ratio_plot,
        delta_plot   = delta_plot,
        q_vent_rate_plot = q_vent_rate_plot,
        emission_gph_plot = emission_gph_plot,
        emission_kgpy_plot = emission_kgpy_plot
)

# Corresponding size settings (width, height)
plot_sizes <- list(
        c_mgm3_plot = c(17, 8.5),
        c_ppm_plot  = c(17, 8.5),
        ratio_plot  = c(16.5, 6),
        delta_plot  = c(12, 8.5),
        q_vent_rate_plot = c(12, 4),
        emission_gph_plot = c(14, 6.5),
        emission_kgpy_plot = c(14, 6.5)
)

# Save each plot using its specific size
for (plot_name in names(dailyplots)) {
        size <- plot_sizes[[plot_name]]
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = dailyplots[[plot_name]],
                width = size[1],
                height = size[2],
                dpi = 300
        )
}

########### Calculate Statistical Summary (mean, SD and CV) #########
# Concentration and relative error
c_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3",
                          "err_CO2_mgm3", "err_CH4_mgm3", "err_NH3_mgm3"),
        var_type_filter = c("concentration mgm3", "concentration mgm3 relative error"),
        group_vars = c("analyzer", "location"),
        time_group = "hour")

# Delta concentrations and relative error
d_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("delta_CO2", "delta_CH4", "delta_NH3",
                          "err_delta_CO2", "err_delta_CH4", "err_delta_NH3"),
        var_type_filter = c("concentration delta", "concentration delta relative error"),
        group_vars = c("analyzer", "location"),
        time_group = "hour")

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
        time_group = "hour")

# Ventilation rates and error
q_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("Q_vent", "err_Q_vent"),
        var_type_filter = c("ventilation rate", "ventilation rate error"),
        group_vars = c("analyzer", "location"),
        time_group = "hour")

# Emissions per hour / year and relative error
e_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("e_CH4_g_h", "e_CH4_kg_yr", "e_NH3_g_h", "e_NH3_kg_yr",
                          "err_e_CH4_g_h", "err_e_CH4_kg_yr", "err_e_NH3_g_h", "err_e_NH3_kg_yr"),
        var_type_filter = c("emission gram per hour", "emission per year",
                            "emission per hour relative error", "emission per year relative error"),
        group_vars = c("analyzer", "location"),
        time_group = "hour")


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

########## Heat maps (CV) #############
min(q_stat_sum$cv[
        q_stat_sum$var == "Q_vent" & 
                q_stat_sum$location == "North background" & 
                q_stat_sum$analyzer != "FTIR.4"], na.rm = TRUE)

max(q_stat_sum$cv[
        q_stat_sum$var == "Q_vent" & 
                q_stat_sum$location == "North background" & 
                q_stat_sum$analyzer != "FTIR.4"], na.rm = TRUE)

min(q_stat_sum$cv[
        q_stat_sum$var == "Q_vent" & 
                q_stat_sum$location == "South background" & 
                q_stat_sum$analyzer != "FTIR.4"], na.rm = TRUE)

max(q_stat_sum$cv[
        q_stat_sum$var == "Q_vent" & 
                q_stat_sum$location == "South background" & 
                q_stat_sum$analyzer != "FTIR.4"], na.rm = TRUE)

q_heatmap <- emiheatmap(q_stat_sum %>% filter(analyzer != "FTIR.4"),
                        y = "Q_vent",
                        time_group = "hour")

e_CH4_heatmap <- emiheatmap(e_stat_sum %>% filter(analyzer != "FTIR.4"),
                       y = "e_CH4_g_h",
                       time_group = "hour")

e_NH3_heatmap <- emiheatmap(e_stat_sum %>% filter(analyzer != "FTIR.4"),
                            y = "e_NH3_g_h",
                            time_group = "hour")

# Named list of plots
dailyplots <- list(
        q_heatmap = q_heatmap,
        e_CH4_heatmap = e_CH4_heatmap,
        e_NH3_heatmap = e_NH3_heatmap)

# coresponding size settings (width, height)
plot_sizes <- list(
        q_heatmap = c(10, 4.4),
        e_CH4_heatmap = c(10, 4.4),
        e_NH3_heatmap = c(10, 4.4))

# Save each plot using its specific size
for (plot_name in names(dailyplots)) {
        size <- plot_sizes[[plot_name]]
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = dailyplots[[plot_name]],
                width = size[1], height = size[2], dpi = 300)
}

########## Box plots (HSD) ventilation and emission rates ##############
q_boxplot <- emiboxplot(
        data = emission_reshaped,
        y = "Q_vent",
        group_var = "analyzer")

e_CH4_boxplot <- emiboxplot(
        data = emission_reshaped,
        y = "e_CH4_g_h",
        group_var = "analyzer")

e_NH3_boxplot <- emiboxplot(
        data = emission_reshaped,
        y = "e_NH3_g_h",
        group_var = "analyzer")


# Save e_boxplot
ggsave(filename = "q_e_boxplot.pdf",
       plot = q_e_boxplot,
       width = 10.5, height = 8.5, dpi = 300)


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
