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
                        t_factor = 1 + 4e-5 * (20 - temp_in)^3,
                        phi_T_cor = phi * t_factor,
                        
                        # Relative animal activity correction
                        A_cor = 1 - a * (sin((2*pi/24) * (hour + 6 - h_min))),
                        
                        # Heat production per cow corrected for T and activity
                        hpu_T_A_cor_per_cow = phi_T_cor * A_cor,
                        
                        #PCO2
                        PCO2 = (0.185 * hpu_T_A_cor_per_cow) * 1000,
                        
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
                        Q_vent_N = ifelse(delta_CO2_N != 0, PCO2 / delta_CO2_N, NA_real_),
                        Q_vent_S = ifelse(delta_CO2_S != 0, PCO2 / delta_CO2_S, NA_real_),
                        
                        # Instantaneous emissions (g/h) divided by 1000 to convert mg to g
                        e_NH3_gh_N = (delta_NH3_N * Q_vent_N / 1000) * n_dairycows_in, 
                        e_CH4_gh_N = (delta_CH4_N * Q_vent_N / 1000) * n_dairycows_in,
                        e_NH3_gh_S = (delta_NH3_S * Q_vent_S / 1000) * n_dairycows_in,
                        e_CH4_gh_S = (delta_CH4_S * Q_vent_S / 1000) * n_dairycows_in,
                        
                        # Annual emissions (kg/year) divided by 1000 to convert g to kg
                        e_NH3_ghLU_N = (e_NH3_gh_N * 500) / (n_dairycows_in * m_weight),
                        e_CH4_ghLU_N = (e_CH4_gh_N * 500) / (n_dairycows_in * m_weight),
                        e_NH3_ghLU_S = (e_NH3_gh_S * 500) / (n_dairycows_in * m_weight),
                        e_CH4_ghLU_S = (e_CH4_gh_S * 500) / (n_dairycows_in * m_weight)
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
        
# Development of pivot longer function
reshaper <- function(df) {
        library(dplyr)
        library(tidyr)
        library(stringr)
        library(lubridate)
        
        meta_cols <- c("DATE.TIME", "analyzer")
        
        # Identify numeric measurement columns
        measure_cols <- df %>%
                select(-all_of(meta_cols)) %>%
                select(where(is.numeric)) %>%
                names()
        
        # Pivot to long format
        df_long <- df %>%
                pivot_longer(
                        cols = all_of(measure_cols),
                        names_to = "var",
                        values_to = "value"
                ) %>%
                mutate(
                        location = case_when(
                                str_detect(var, "_in$") ~ "Barn inside",
                                str_detect(var, "_N$")  ~ "North background",
                                str_detect(var, "_S$")  ~ "South background",
                                TRUE ~ NA_character_
                        ),
                        var = str_remove(var, "_(in|N|S)$"),
                        DATE.TIME = as.POSIXct(DATE.TIME),
                        day  = factor(as.Date(DATE.TIME)),
                        hour = factor(format(DATE.TIME, "%H:%M"))
                ) %>%
                select(DATE.TIME, day, hour, location, analyzer, var, value) %>%
                arrange(DATE.TIME, var, analyzer, location) %>%
                # Map special analyzers
                mutate(analyzer = case_when(
                        var %in% c("temp", "RH")                           ~ "HOBO",
                        var %in% c("wd_mst", "ws_mst", "wd_trv", "ws_trv") ~ "USA",
                        var %in% c("n_dairycows")                          ~ "RGB",
                        TRUE                                               ~ analyzer
                ))
        
        # ---- Add baseline per DATE.TIME, location, var ----
        baseline_df <- df_long %>%
                group_by(DATE.TIME, location, var) %>%
                summarise(
                        value = mean(value, na.rm = TRUE),
                        day   = first(day),    # copy day
                        hour  = first(hour),   # copy hour
                        .groups = "drop"
                ) %>%
                mutate(analyzer = "baseline")
        
        df_long <- bind_rows(df_long, baseline_df) %>%
                mutate(analyzer = factor(analyzer,
                                         levels = c("FTIR.1","FTIR.2","FTIR.3","FTIR.4",
                                                    "CRDS.1","CRDS.2","CRDS.3",
                                                    "HOBO","USA","RGB","baseline"))) %>%
                arrange(DATE.TIME, location, var)
        
        return(df_long)
}

# Development of function stat_table to calculate baseline, Mean, SD, CV, RE, ICC 
stat_error_table <- function(df_long, time.group = "hour", analyzer.levels = NULL) {
        library(dplyr)
        library(DescTools)
        
        # ---- Ensure day and hour columns exist ----
        if (!("day" %in% names(df_long))) 
                df_long <- df_long %>% mutate(day = as.Date(DATE.TIME))
        if (!("hour" %in% names(df_long))) 
                df_long <- df_long %>% mutate(hour = format(DATE.TIME, "%H"))
        
        # ---- Compute % error relative to baseline ----
        df_long <- df_long %>%
                group_by(DATE.TIME, location, var) %>%
                mutate(
                        baseline_value = value[analyzer == "baseline"],
                        pct_err = 100 * (value - baseline_value) / baseline_value
                ) %>%
                ungroup() %>%
                select(-baseline_value)
        
        # ---- Mean, SD, CV + Mean pct_err calculation ----
        result_df <- df_long %>%
                group_by(across(all_of(time.group)), analyzer, location, var) %>%
                summarise(
                        mean_value = mean(value, na.rm = TRUE),
                        sd         = sd(value, na.rm = TRUE),
                        cv         = DescTools::CoefVar(value, na.rm = TRUE) * 100,
                        pct_err    = mean(pct_err, na.rm = TRUE),
                        .groups = "drop"
                )
        
        # ---- Apply analyzer levels if provided ----
        if (!is.null(analyzer.levels)) {
                result_df <- result_df %>%
                        mutate(analyzer = factor(analyzer, levels = analyzer.levels)) %>%
                        filter(analyzer %in% analyzer.levels)
        }
        
        return(result_df %>% arrange(across(all_of(time.group)), location, var, analyzer))
}

icc_table <- function(df_long, vars = NULL, time.group = "day") {
        library(dplyr)
        library(tidyr)
        library(irr)
        
        df <- df_long
        
        # Ensure day/hour columns exist
        if (!"day" %in% names(df)) df <- df %>% mutate(day = as.Date(DATE.TIME))
        if (!"hour" %in% names(df)) df <- df %>% mutate(hour = format(DATE.TIME, "%H"))
        
        # Filter variables if provided
        if (!is.null(vars)) df <- df %>% filter(var %in% vars)
        
        results <- df %>%
                group_by(across(all_of(c(time.group, "location", "var")))) %>%
                summarise(
                        icc_result = list({
                                # Pivot wider: analyzers as columns
                                wide <- pivot_wider(cur_data(), names_from = analyzer, values_from = value)
                                
                                # Keep only analyzer columns for ICC
                                analyzer_cols <- setdiff(names(wide), c(time.group, "location", "var", "DATE.TIME", "hour"))
                                mat <- wide[, analyzer_cols, drop = FALSE]
                                
                                # Only compute ICC if we have >=2 analyzers and >=2 rows
                                if (ncol(mat) > 1 & nrow(mat) > 1) {
                                        suppressWarnings(irr::icc(mat, model = "twoway", type = "consistency", unit = "average"))
                                } else {
                                        NA
                                }
                        }),
                        .groups = "drop"
                ) %>%
                rowwise() %>%
                mutate(
                        ICC = if (is.list(icc_result)) round(icc_result$value, 2) else NA_real_,
                        F   = if (is.list(icc_result)) round(icc_result$Fvalue, 2) else NA_real_,
                        p   = if (is.list(icc_result) && !is.na(icc_result$p.value)) {
                                if (icc_result$p.value < 0.01) {
                                        formatC(icc_result$p.value, format = "e", digits = 2)
                                } else {
                                        as.character(round(icc_result$p.value, 2))
                                }
                        } else NA_character_,
                        signif = if (is.list(icc_result) && !is.na(icc_result$p.value)) case_when(
                                icc_result$p.value < 0.001 ~ "***",
                                icc_result$p.value < 0.01  ~ "**",
                                icc_result$p.value < 0.05  ~ "*",
                                TRUE                       ~ "ns"
                        ) else NA_character_
                ) %>%
                select(all_of(time.group), location, var, ICC, F, p, signif) %>%
                ungroup()
        
        return(results)
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
emiconplot <- function(data, y = NULL, location_filter = NULL, plot_err = FALSE, x = "DATE.TIME") {
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
        
        # Ensure the selected x column exists
        if (!x %in% names(data)) {
                stop("x column '", x, "' not found in data. Available columns: ", paste(names(data), collapse = ", "))
        }
        
        # Choose value column
        value_col <- if (plot_err) "err" else "value"
        
        # Summary statistics by grouping
        summary_data <- data %>%
                group_by(.data[[x]], analyzer, location, variable) %>%
                summarise(
                        mean_val = mean(.data[[value_col]], na.rm = TRUE),
                        sd_val   = sd(.data[[value_col]], na.rm = TRUE),
                        .groups = "drop"
                )
        
        # Smart facet labels
        facet_labels_expr <- c(
                "CO2_mgm3"     = "c[CO2]~'(mg '*m^-3*')'",
                "CH4_mgm3"     = "c[CH4]~'(mg '*m^-3*')'",
                "NH3_mgm3"     = "c[NH3]~'(mg '*m^-3*')'",
                "r_CH4/CO2"    = "c[CH4]/c[CO2]~'('*'%'*')'",
                "r_NH3/CO2"    = "c[NH3]/c[CO2]~'('*'%'*')'",
                "delta_CO2"    = "Delta*c[CO2]~'(mg '*m^-3*')'",
                "delta_CH4"    = "Delta*c[CH4]~'(mg '*m^-3*')'",
                "delta_NH3"    = "Delta*c[NH3]~'(mg '*m^-3*')'",
                "Q_vent"       = "Q~'('*m^3~h^-1~LU^-1*')'",
                "e_CH4_gh"     = "e[CH4]~'(g '*h^-1*')'",
                "e_NH3_gh"     = "e[NH3]~'(g '*h^-1*')'",
                "e_CH4_ghLU"   = "e[CH4]~'(g '*h^-1~LU^-1*')'",
                "e_NH3_ghLU"   = "e[NH3]~'(g '*h^-1~LU^-1*')'",
                "temp"         = "Temperature " ~ "(°C)",
                "RH"           = "Relative Humidity " ~ "(%)",
                "ws_mst"       = "Wind Speed Mast " ~ "(m " *s^-1* ")",
                "wd_mst"       = "Wind Diection Mast " ~ "(°)",
                "wd_trv"       = "Wind Direction Traverse " ~ "(°)",
                "ws_trv"       = "Wind Speed Traverse " ~ "(m " *s^-1* ")",
                "n_dairycows"  = "'Number of Cows'"
        )
        
        missing_vars <- setdiff(y, names(facet_labels_expr))
        if (length(missing_vars) > 0) {
                stop("These variables are missing facet labels: ", paste(missing_vars, collapse = ", "))
        }
        
        summary_data <- summary_data %>%
                mutate(facet_label = factor(
                        facet_labels_expr[as.character(variable)],
                        levels = facet_labels_expr[y]
                ))
        
        # Known colors & shapes
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
        
        # Fill all analyzers in data with defaults if missing
        all_analyzers <- unique(summary_data$analyzer)
        analyzer_colors_full <- setNames(rep("black", length(all_analyzers)), all_analyzers)
        analyzer_shapes_full <- setNames(rep(16, length(all_analyzers)), all_analyzers)
        analyzer_colors_full[names(analyzer_colors)] <- analyzer_colors
        analyzer_shapes_full[names(analyzer_shapes)] <- analyzer_shapes
        
        # Labels
        ylab_to_use <- if (plot_err) "Relative error (%)" else "Mean ± SD"
        xlab_to_use <- x
        
        # Base ggplot
        p <- ggplot(summary_data, aes(x = .data[[x]], y = mean_val,
                                      color = analyzer, shape = analyzer, group = analyzer))
        
        # Add geoms depending on x
        if (x %in% c("day", "hour")) {
                # Points + errorbars with dodge, no line
                p <- p +
                        geom_point(size = 2, position = position_dodge(width = 0.5)) +
                        geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                                      width = 0.4, position = position_dodge(width = 0.5))
        } else {
                # Line + points + errorbars
                p <- p +
                        geom_line(size = 0.5, alpha = 0.6) +
                        geom_point(size = 2) +
                        geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                                      width = 0.4)
        }
        
        # Apply full color & shape scales
        p <- p +
                scale_color_manual(values = analyzer_colors_full) +
                scale_shape_manual(values = analyzer_shapes_full)
        
        # Facets, scales, theme
        p <- p +
                facet_grid(facet_label ~ ., scales = "free_y", switch = "y",
                           labeller = label_parsed) +
                scale_y_continuous(breaks = pretty_breaks(n = 7)) +
                labs(x = xlab_to_use, y = ylab_to_use, title = unique(summary_data$location)) +
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
        
        # Add datetime scale only when x = DATE.TIME
        if (x == "DATE.TIME") {
                x_breaks <- seq(
                        from = min(summary_data$DATE.TIME, na.rm = TRUE),
                        to   = max(summary_data$DATE.TIME, na.rm = TRUE),
                        by   = "6 hours"
                )
                p <- p + scale_x_datetime(breaks = x_breaks, date_labels = "%Y-%m-%d %H:%M")
        }
        
        print(p)
        return(p)
}

# Development of function for correlogram
emicorrgram <- function(data, target_variable, locations = NULL) {
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(Hmisc)
        
        # Step 1: Filter for the target variable
        filtered <- data %>%
                filter(var == target_variable) %>%
                select(DATE.TIME, location, analyzer, value)
        
        # Step 2: Filter by user-selected locations if provided
        if (!is.null(locations)) {
                filtered <- filtered %>% filter(location %in% locations)
        }
        
        # Step 3: Remove duplicates by averaging
        filtered <- filtered %>%
                group_by(DATE.TIME, location, analyzer) %>%
                summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
        
        # Step 4: Pivot wider (each analyzer becomes a column)
        pivoted <- filtered %>%
                pivot_wider(
                        names_from = analyzer,
                        values_from = value,
                        values_fn = function(x, ...) mean(x, na.rm = TRUE)
                )
        
        # Step 5: Keep only numeric columns for correlation
        subdata_num <- pivoted %>%
                select(-DATE.TIME, -location) %>%
                mutate(across(everything(), as.numeric))
        
        # Check there are at least 2 analyzers
        if (ncol(subdata_num) < 2) {
                stop("Need at least 2 analyzers with numeric data to compute correlations")
        }
        
        # Step 6: Compute correlation and p-value
        cor_res <- Hmisc::rcorr(as.matrix(subdata_num), type = "pearson")
        cor_mat <- cor_res$r
        p_mat <- cor_res$P
        col_means <- colMeans(subdata_num, na.rm = TRUE)
        
        # Step 7: Prepare long-format data for plotting
        cor_long <- expand.grid(Var1 = colnames(cor_mat),
                                Var2 = colnames(cor_mat),
                                stringsAsFactors = FALSE) %>%
                mutate(
                        correlation = as.vector(cor_mat),
                        pvalue = as.vector(p_mat),
                        mean_value = col_means[Var1],
                        p_text = ifelse(pvalue > 0.05, "ns", "")
                ) %>%
                # Keep only lower triangle
                filter(as.numeric(factor(Var1)) > as.numeric(factor(Var2)))
        
        # Step 8: Title labels (parsed)
        title_labels_expr <- c(
                "CO2_mgm3"     = "c[CO2]~'(mg '*m^-3*')'",
                "CH4_mgm3"     = "c[CH4]~'(mg '*m^-3*')'",
                "NH3_mgm3"     = "c[NH3]~'(mg '*m^-3*')'",
                "r_CH4/CO2"    = "c[CH4]/c[CO2]~'('*'%'*')'",
                "r_NH3/CO2"    = "c[NH3]/c[CO2]~'('*'%'*')'",
                "delta_CO2"    = "Delta*c[CO2]~'(mg '*m^-3*')'",
                "delta_CH4"    = "Delta*c[CH4]~'(mg '*m^-3*')'",
                "delta_NH3"    = "Delta*c[NH3]~'(mg '*m^-3*')'",
                "Q_vent"       = "Q~'('*m^3~h^-1~LU^-1*')'",
                "e_CH4_gh"     = "e[CH4]~'(g '*h^-1*')'",
                "e_NH3_gh"     = "e[NH3]~'(g '*h^-1*')'",
                "e_CH4_ghLU"   = "e[CH4]~'(g '*h^-1~LU^-1*')'",
                "e_NH3_ghLU"   = "e[NH3]~'(g '*h^-1~LU^-1*')'",
                "temp"         = "Temperature " ~ "(°C)",
                "RH"           = "Relative Humidity " ~ "(%)",
                "ws_mst"       = "Wind Speed Mast " ~ "(m " *s^-1* ")",
                "wd_mst"       = "Wind Diection Mast " ~ "(°)",
                "wd_trv"       = "Wind Direction Traverse " ~ "(°)",
                "ws_trv"       = "Wind Speed Traverse " ~ "(m " *s^-1* ")",
                "n_dairycows"  = "'Number of Cows'"
        )
        
        # Build title: variable label + location
        title_label <- if (target_variable %in% names(title_labels_expr)) {
                title_labels_expr[[target_variable]]
        } else {
                target_variable
        }
        
        location_name <- if (!is.null(locations)) {
                paste(locations, collapse = ", ")
        } else {
                unique(filtered$location)
        }
        
        # Parse expression for ggplot
        full_title_expr <- parse(text = paste0(title_label, " ~ '(' ~ '", location_name, "' ~ ')'"))
        
        # Step 9: Plot
        p <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = correlation)) +
                geom_tile(color = "white") +
                geom_text(aes(label = paste0(
                        round(correlation, 2),
                        ifelse(pvalue <= 0.001, "\n***",
                               ifelse(pvalue <= 0.01, "\n**",
                                      ifelse(pvalue <= 0.05, "\n*", "\nns")))
                )),
                size = 4, color = "black") +
                scale_y_discrete(position = "right") + 
                scale_fill_gradientn(
                        colors = c("darkred","red","orange2","gold2","yellow","greenyellow","green1","green3","darkgreen"),
                        limits = c(0,1),
                        name = "PCC"
                ) +
                labs(title = full_title_expr, x = NULL, y = NULL) +
                theme_classic() +
                theme(
                        plot.title = element_text(hjust = 0.5),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        axis.text.y = element_text(angle = 45, hjust = 1),
                        legend.position = "bottom",
                        panel.border = element_rect(colour = "black", fill = NA)
                )
        
        print(p)
        return(p)
}

# Development of function cv emiheatmap
emiheatmap <- function(stat_df, vars = NULL, time.group = "hour", location.filter = NULL) {
        library(dplyr)
        library(ggplot2)
        
        df <- stat_df
        
        # Filter by variable(s) if requested
        if (!is.null(vars) && "var" %in% colnames(df)) {
                df <- df %>% filter(var %in% vars)
        }
        
        # Filter by location if requested
        if (!is.null(location.filter) && "location" %in% colnames(df)) {
                df <- df %>% filter(location %in% location.filter)
        }
        
        # Create hour/day if needed
        if (!"hour" %in% colnames(df) & "DATE.TIME" %in% colnames(df)) {
                df <- df %>% mutate(hour = format(DATE.TIME, "%H"))
        }
        if (!"day" %in% colnames(df) & "DATE.TIME" %in% colnames(df)) {
                df <- df %>% mutate(day = format(DATE.TIME, "%Y-%m-%d"))
        }
        
        # Ensure analyzer is factor for y-axis
        df$analyzer <- factor(df$analyzer, levels = unique(df$analyzer))
        
        # Title expressions (same as in emicorrgram)
        title_labels_expr <- c(
                "CO2_mgm3"     = "c[CO2]~'(mg '*m^-3*')'",
                "CH4_mgm3"     = "c[CH4]~'(mg '*m^-3*')'",
                "NH3_mgm3"     = "c[NH3]~'(mg '*m^-3*')'",
                "r_CH4/CO2"    = "c[CH4]/c[CO2]~'('*'%'*')'",
                "r_NH3/CO2"    = "c[NH3]/c[CO2]~'('*'%'*')'",
                "delta_CO2"    = "Delta*c[CO2]~'(mg '*m^-3*')'",
                "delta_CH4"    = "Delta*c[CH4]~'(mg '*m^-3*')'",
                "delta_NH3"    = "Delta*c[NH3]~'(mg '*m^-3*')'",
                "Q_vent"       = "Q~'('*m^3~h^-1~LU^-1*')'",
                "e_CH4_gh"     = "e[CH4]~'(g '*h^-1*')'",
                "e_NH3_gh"     = "e[NH3]~'(g '*h^-1*')'",
                "e_CH4_ghLU"   = "e[CH4]~'(g '*h^-1~LU^-1*')'",
                "e_NH3_ghLU"   = "e[NH3]~'(g '*h^-1~LU^-1*')'",
                "temp"         = "Temperature " ~ "(°C)",
                "RH"           = "Relative Humidity " ~ "(%)",
                "ws_mst"       = "Wind Speed Mast " ~ "(m " *s^-1* ")",
                "wd_mst"       = "Wind Diection Mast " ~ "(°)",
                "wd_trv"       = "Wind Direction Traverse " ~ "(°)",
                "ws_trv"       = "Wind Speed Traverse " ~ "(m " *s^-1* ")",
                "n_dairycows"  = "'Number of Cows'"
        )
        
        # Pick target variable(s) for title
        target_var <- if (!is.null(vars)) vars[1] else unique(df$var)[1]
        
        # Build title label
        title_label <- if (target_var %in% names(title_labels_expr)) {
                title_labels_expr[[target_var]]
        } else {
                target_var
        }
        
        # Add location info
        location_name <- if (!is.null(location.filter)) {
                paste(location.filter, collapse = ", ")
        } else if ("location" %in% colnames(df)) {
                paste(unique(df$location), collapse = ", ")
        } else {
                ""
        }
        
        # Parse full title expression
        full_title_expr <- parse(text = paste0(title_label, " ~ '(' ~ '", location_name, "' ~ ')'"))
        
        # Plot heatmap
        p <- ggplot(df, aes_string(x = time.group, y = "analyzer", fill = "cv")) +
                geom_tile(color = "white") +
                geom_text(aes(label = paste0(round(cv, 1))), size = 3, color = "black") +
                scale_fill_gradientn(
                        colors = c("darkgreen","green4","green2","yellow","orange","red","darkred"),
                        limits = c(0, 150),
                        breaks = seq(0, 150, by = 10),
                        name = "CV (%)",
                        na.value = "grey80"
                ) +
                labs(x = time.group, y = "Analyzer", title = full_title_expr) +
                theme_classic() +
                theme(axis.text.x = element_text(
                        size = 12, 
                        angle = 45, 
                        hjust = 1),
                      axis.text.y = element_text(size = 12),
                      panel.grid = element_blank(),
                      legend.position = "bottom",
                      legend.key.width = unit(2, "cm"),
                      panel.border = element_rect(color = "black", fill = NA),
                      plot.title = element_text(hjust = 0.5)
                )
        
        print(p)
        return(p)
}

# Development of function Bland-altman Plot
bland_altman_plot <- function(data, var_filter, analyzer_pair, location_filter = NULL, x = "DATE.TIME") {
        library(dplyr)
        library(ggplot2)
        library(scales)
        
        # mapping for variable labels
        var_labels <- c(
                "CO2_mgm3"     = "c[CO2]~'(mg '*m^-3*')'",
                "CH4_mgm3"     = "c[CH4]~'(mg '*m^-3*')'",
                "NH3_mgm3"     = "c[NH3]~'(mg '*m^-3*')'",
                "r_CH4/CO2"    = "c[CH4]/c[CO2]~'('*'%'*')'",
                "r_NH3/CO2"    = "c[NH3]/c[CO2]~'('*'%'*')'",
                "delta_CO2"    = "Delta*c[CO2]~'(mg '*m^-3*')'",
                "delta_CH4"    = "Delta*c[CH4]~'(mg '*m^-3*')'",
                "delta_NH3"    = "Delta*c[NH3]~'(mg '*m^-3*')'",
                "Q_vent"       = "Q~'('*m^3~h^-1~LU^-1*')'",
                "e_CH4_gh"     = "e[CH4]~'(g '*h^-1*')'",
                "e_NH3_gh"     = "e[NH3]~'(g '*h^-1*')'",
                "e_CH4_ghLU"   = "e[CH4]~'(g '*h^-1~LU^-1*')'",
                "e_NH3_ghLU"   = "e[NH3]~'(g '*h^-1~LU^-1*')'",
                "temp"         = "'Temperature (°C)'",
                "RH"           = "'Relative Humidity (%)'",
                "ws_mst"       = "'Wind Speed Mast (m s^-1)'",
                "wd_mst"       = "'Wind Direction Mast (°)'",
                "wd_trv"       = "'Wind Direction Traverse (°)'",
                "ws_trv"       = "'Wind Speed Traverse (m s^-1)'",
                "n_dairycows"  = "'Number of Cows'"
        )
        
        var_label_expr <- parse(text = var_labels[[var_filter]])[[1]]
        
        if ("var" %in% names(data) && !"variable" %in% names(data)) {
                data <- data %>% rename(variable = var)
        }
        
        df <- data %>%
                filter(variable == var_filter)
        
        if (!is.null(location_filter)) {
                df <- df %>% filter(location %in% location_filter)
        }
        
        df <- df %>%
                filter(analyzer %in% analyzer_pair)
        
        df_wide <- df %>%
                select(all_of(c(x, "location", "analyzer", "value"))) %>%
                tidyr::pivot_wider(names_from = analyzer, values_from = value)
        
        a1 <- analyzer_pair[1]
        a2 <- analyzer_pair[2]
        
        df_ba <- df_wide %>%
                mutate(
                        mean_val = ( .data[[a1]] + .data[[a2]] ) / 2,
                        diff_val = .data[[a1]] - .data[[a2]]
                )
        
        bias <- mean(df_ba$diff_val, na.rm = TRUE)
        sd_diff <- sd(df_ba$diff_val, na.rm = TRUE)
        loa_upper <- bias + 1.96 * sd_diff
        loa_lower <- bias - 1.96 * sd_diff
        
        subtitle_expr <- if (!is.null(location_filter)) {
                bquote(.(var_label_expr) ~ "|" ~ .(location_filter))
        } else {
                var_label_expr
        }
        
        p <- ggplot(df_ba, aes(x = mean_val, y = diff_val)) +
                geom_point(alpha = 0.6) +
                geom_hline(yintercept = bias, color = "blue", linetype = "dashed", size = 1) +
                geom_hline(yintercept = loa_upper, color = "red", linetype = "dotted", size = 1) +
                geom_hline(yintercept = loa_lower, color = "red", linetype = "dotted", size = 1) +
                scale_x_continuous(
                        breaks = pretty_breaks(n = 8),
                        labels = number_format(accuracy = 0.1)
                ) +
                scale_y_continuous(
                        breaks = pretty_breaks(n = 8),
                        labels = number_format(accuracy = 0.1)
                ) +
                ggtitle(
                        bquote(),
                        subtitle = subtitle_expr
                ) +
                labs(
                        x = bquote("Mean of" ~ .(a1) ~ "and" ~ .(a2)),
                        y = bquote("Difference (" ~ .(a1) - .(a2) ~ ")")
                ) +
                theme_classic() +
                theme(
                        plot.subtitle = element_text(hjust = 0.5),
                        plot.title = element_text(hjust = 0.5)
                )
        
        print(p)
        return(p)
}

# Development of Linear mixed model function
linear_mixed_model <- function(var_name, data) {
        library(lme4)
        library(broom.mixed)
        library(dplyr)
        library(stringr)
        library(ggplot2)
        
        # Ensure day column exists and is factor
        data <- data %>%
                mutate(day = factor(as.Date(DATE.TIME)))
        
        # Detect location from variable name
        if (str_detect(var_name, "_in$")) {
                location <- "in"
        } else if (str_detect(var_name, "_N$")) {
                location <- "N"
        } else if (str_detect(var_name, "_S$")) {
                location <- "S"
        } else {
                stop("Cannot detect location from variable name. Use _in, _N, or _S.")
        }
        
        # Predictors
        temp_var <- ifelse(location == "in", "temp_in", "temp_N")
        rh_var   <- ifelse(location == "in", "RH_in", "RH_N")
        cows_var <- "n_dairycows_in"
        wd_var <- ifelse(location %in% c("N","S"), "wd_mst_S", NULL)
        ws_var <- ifelse(location %in% c("N","S"), "ws_mst_S", NULL)
        
        predictors <- c(temp_var, rh_var, wd_var, ws_var, cows_var)
        predictors <- predictors[predictors %in% colnames(data)]
        
        formula <- as.formula(
                paste(var_name, "~", paste(predictors, collapse = " + "), "+ (1 | analyzer) + (1 | day)")
        )
        
        # Fit the model
        model <- lmer(formula, data = data, REML = TRUE)
        
        # Fixed effects tibble
        fixed_tbl <- broom.mixed::tidy(model, effects = "fixed")
        
        # Random effects
        ranefs <- ranef(model)
        analyzer_effects <- ranefs$analyzer %>%
                mutate(analyzer = rownames(.)) %>%
                rename(intercept_adj = `(Intercept)`) %>%
                mutate(percent_dev = 100 * intercept_adj / mean(fixed_tbl$estimate[fixed_tbl$term=="(Intercept)"]))
        
        day_effects <- ranefs$day %>%
                mutate(day = rownames(.)) %>%
                rename(intercept_adj = `(Intercept)`) %>%
                mutate(percent_dev = 100 * intercept_adj / mean(fixed_tbl$estimate[fixed_tbl$term=="(Intercept)"]))
        
        # Print results
        cat("\nFixed Effects:\n")
        print(fixed_tbl)
        cat("\nRandom Effects - Analyzer:\n")
        print(analyzer_effects)
        cat("\nRandom Effects - Day:\n")
        print(day_effects)
        
        # Heatmap for analyzer deviations
        heatmap_data <- analyzer_effects %>%
                left_join(day_effects, by = character()) %>%
                select(analyzer, day, percent_dev = percent_dev.y)
        
        ggplot(heatmap_data, aes(x = day, y = analyzer, fill = percent_dev)) +
                geom_tile(color = "black") +
                scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
                labs(title = paste("Analyzer Deviation (%) for", var_name), x = "Day", y = "Analyzer") +
                theme_classic() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        invisible(list(
                formula = formula,
                model_summary = summary(model),
                fixed_effects = fixed_tbl,
                random_effects_analyzer = analyzer_effects,
                random_effects_day = day_effects
        ))
}

######## Import Data #########
# Read all gas data
ATB_FTIR <- read.csv("20250408-14_ATB_wide_FTIR.1.csv")
LUFA_FTIR <- read.csv("20250408-14_LUFA_wide_FTIR.2.csv")
MBBM_FTIR <- read.csv("20250408-14_MBBM_wide_FTIR.4.csv")
ANECO_FTIR <- read.csv("20250408-14_ANECO_wide_FTIR.4.csv")
ATB_CRDS <- read.csv("20250408-14_ATB_wide_CRDS.1.csv")
UB_CRDS <- read.csv("20250408-14_UB_wide_CRDS.2.csv")
LUFA_CRDS <- read.csv("20250408-14_LUFA_wide_CRDS.3.csv")

# Define your start and end time as POSIXct
start_time <- "2025-04-08 12:00:00"
end_time   <- "2025-04-14 12:00:00"
                         
# Combine all data set
gas_data <- bind_rows(LUFA_FTIR, ANECO_FTIR, MBBM_FTIR, ATB_FTIR, ATB_CRDS, LUFA_CRDS, UB_CRDS) %>%
        mutate(DATE.TIME = ymd_hms(DATE.TIME, tz = "UTC")) %>%
        filter(DATE.TIME >= ymd_hms(start_time) &
                       DATE.TIME <= ymd_hms(end_time)) %>%
        select(DATE.TIME, analyzer, everything(), -lab) %>%
        arrange(DATE.TIME) %>% distinct()

# Count observations per analyzer
gas_data %>% summarise(across(where(is.numeric), ~ sum(!is.na(.)), .names = "{.col}"), .by = analyzer)
                
# Read animal data and fix time zone
animal_data <- read.csv("20250408-15_LVAT_Animal_data.csv") %>%
        mutate(DATE.TIME = ymd_hms(DATE.TIME, tz = "UTC")) %>%
        filter(DATE.TIME >= ymd_hms(start_time) &
                       DATE.TIME <= ymd_hms(end_time)) %>%
        group_by(DATE.TIME) %>%
        summarise(across(everything(), ~ first(.x)), .groups = "drop") %>%
        rename(n_dairycows_in = n_dairycows) 	

# Read weather data and fix time zone
T_RH_HOBO <- read.csv("T_RH_08_04_2025_To_30_06_2025.csv", stringsAsFactors = FALSE) %>%
        mutate(DATE.TIME = paste(date, time)) %>%
        mutate(DATE.TIME = ymd_hms(DATE.TIME, tz = "UTC")) %>%
        select(DATE.TIME, everything(), -date, -time)%>%
        rename(temp_N = temp_outside,
               RH_N   = RH_out,
               temp_in  = temp_inside,
               RH_in    = RH_inside)

USA_Mst <- read.csv("2025_04_08_USA_Mst_5_min_avg.csv", stringsAsFactors = FALSE) %>%
        mutate(DATE.TIME = ymd_hms(DATE.TIME, tz = "UTC")) %>%
        filter(DATE.TIME >= ymd_hms(start_time) &
                       DATE.TIME <= ymd_hms(end_time)) %>%
        mutate(DATE.TIME = floor_date(DATE.TIME, "hour")) %>%
        group_by(DATE.TIME) %>%
        summarise(wd_mst_S = mean(wd, na.rm = TRUE),
                  ws_mst_S = mean(ws, na.rm = TRUE),
                  .groups = "drop")

USA_Trv <- read.csv("2025_04_08_USA_Trv_5_min_avg.csv")%>%
        mutate(DATE.TIME = ymd_hms(DATE.TIME, tz = "UTC")) %>%
        filter(DATE.TIME >= ymd_hms(start_time) &
                       DATE.TIME <= ymd_hms(end_time)) %>%
        mutate(DATE.TIME = floor_date(DATE.TIME, "hour")) %>%
        group_by(DATE.TIME) %>%
        summarise(wd_trv_N = mean(wd_direction, na.rm = TRUE),
                  ws_trv_N = mean(wd_speed, na.rm = TRUE),
                  .groups = "drop")

weather_data <- T_RH_HOBO %>%
        full_join(USA_Mst, by = "DATE.TIME") %>%
        full_join(USA_Trv, by = "DATE.TIME") %>%
        arrange(DATE.TIME) 

# Join with gas data
input_combined <- gas_data %>%
        left_join(animal_data, by = "DATE.TIME") %>%
        left_join(weather_data, by = "DATE.TIME") %>%
        arrange(DATE.TIME) %>%
        distinct()

# Write csv
write_excel_csv(input_combined, "20250408-15_ringversuche_input_combined_data.csv")

######## Computation of emissions, ventilation rates and ratios #########
# Calculate emissions and ratios using the function
emission_result <- indirect.CO2.balance(input_combined)

# Remove constants and forumal inputs
emission_result <- emission_result %>% select(-hour, -m_weight, -p_pregnancy_day, 
                                              -Y1_milk_prod, -a, -h_min,
                                              -phi, -t_factor, -phi_T_cor, -A_cor,
                                              -hpu_T_A_cor_per_cow, -PCO2)

# Write csv 
write_excel_csv(emission_result, "20250408-15_ringversuche_emission_result.csv")

######## Reshape Data #########
emission_reshaped <-  reshaper(emission_result)  %>%
        mutate(across(where(is.numeric), ~ round(.x, 2)))
        
# Write csv
write_excel_csv(emission_reshaped, "20250408-15_ringversuche_emission_reshaped.csv")

######## ANOVA and HSD Summary ########
result_HSD_summary <- HSD_table(data = emission_reshaped %>%
                                        filter(analyzer %in% c("FTIR.1",
                                                               "FTIR.2",
                                                               "FTIR.3",
                                                               "FTIR.4",
                                                               "CRDS.1",
                                                               "CRDS.2",
                                                               "CRDS.3")))

write_excel_csv(result_HSD_summary, "result_HSD_summary.csv")

######## Absolute Concentration Trend plots ########
# Concentrations only
c_r_in_mgm3_plot <- emiconplot(
        data = emission_reshaped,
        y = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3", 
              "r_CH4/CO2", "r_NH3/CO2"),
        location_filter = "Barn inside")

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

# Combine all plots in a named list 
all_plots <- list(
        c_r_in_mgm3_plot = c_r_in_mgm3_plot,
        c_r_N_mgm3_plot  = c_r_N_mgm3_plot,
        c_r_S_mgm3_plot  = c_r_S_mgm3_plot)

# Define custom sizes (width, height) for each plot
plot_sizes <- list(
        c_r_in_mgm3_plot   = c(8, 8),
        c_r_N_mgm3_plot    = c(8, 8),
        c_r_S_mgm3_plot    = c(8, 8)
)

# Loop through and save each plot with its custom size
for (plot_name in names(all_plots)) {
        size <- plot_sizes[[plot_name]]
        
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = all_plots[[plot_name]],
                width = size[1],
                height = size[2],
                dpi = 300
        )
}


######## Daily Mean ± SD Trend Plots ########
# North background
d_q_e_day_N_plot <- emiconplot(
        data = emission_reshaped,
        x = "day",
        y = c("delta_CO2", "delta_CH4", "delta_NH3",
              "Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"),
        location_filter = "North background")

# South background
d_q_e_day_S_plot <- emiconplot(
        data = emission_reshaped,
        x = "day",
        y = c("delta_CO2", "delta_CH4", "delta_NH3",
              "Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"),
        location_filter = "South background")

# Combine all plots in a named list 
all_plots <- list(
        d_q_e_day_N_plot = d_q_e_day_N_plot,
        d_q_e_day_S_plot = d_q_e_day_S_plot
)

# Define custom sizes (width, height) for each plot
plot_sizes <- list(
        d_q_e_day_N_plot   = c(9, 9),
        d_q_e_day_S_plot   = c(9, 9)
)

# Loop through and save each plot with its custom size
for (plot_name in names(all_plots)) {
        size <- plot_sizes[[plot_name]]
        
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = all_plots[[plot_name]],
                width = size[1],
                height = size[2],
                dpi = 300
        )
}

######## Hourly Mean ± SD Trend Plots ########
# North background
d_q_e_hour_N_plot <- emiconplot(
        data = emission_reshaped,
        x = "hour",
        y = c("delta_CO2", "delta_CH4", "delta_NH3",
              "Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"),
        location_filter = "North background")

# South background
d_q_e_hour_S_plot <- emiconplot(
        data = emission_reshaped,
        x = "hour",
        y = c("delta_CO2", "delta_CH4", "delta_NH3",
              "Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"),
        location_filter = "South background")


# Combine all plots in a named list 
all_plots <- list(
        d_q_e_hour_N_plot = d_q_e_hour_N_plot,
        d_q_e_hour_S_plot = d_q_e_hour_S_plot
)

# Define custom sizes (width, height) for each plot
plot_sizes <- list(
        d_q_e_hour_N_plot   = c(9, 9),
        d_q_e_hour_S_plot   = c(9, 9)
)

# Loop through and save each plot with its custom size
for (plot_name in names(all_plots)) {
        size <- plot_sizes[[plot_name]]
        
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = all_plots[[plot_name]],
                width = size[1],
                height = size[2],
                dpi = 300
        )
}

######## Weather Plots ########
# Temperature and Relative Humidity 
temp_cows_in_plot <- emiconplot(
        data = emission_reshaped %>% filter(analyzer != "baseline"),
        x = "hour",
        y = c("temp", "RH", "n_dairycows"),
        location_filter = "Barn inside")

temp_cows_N_plot <- emiconplot(
        data = emission_reshaped %>% filter(analyzer != "baseline"),
        x = "hour",
        y = c("temp", "RH"),
        location_filter = "North background")

# Wind speed and direction
weather_N_plot <- emiconplot(
        data = emission_reshaped %>% filter(analyzer != "baseline"),
        x = "hour",
        y = c("wd_trv", "ws_trv"),
        location_filter = "North background")

weather_S_plot <- emiconplot(
        data = emission_reshaped %>% filter(analyzer != "baseline"),
        x = "hour",
        y = c("wd_mst", "ws_mst"),
        location_filter = "South background")


######## Bland-Altman plot ############
#Compare FTIR.1 vs CRDS.1 for delta_CH4 at "North background"
delta_CH4_N_A_blandplot <- bland_altman_plot(data = emission_reshaped,
                                             var_filter = "delta_CH4",
                                             analyzer_pair = c("FTIR.1", "CRDS.1"),
                                             location_filter = "North background")

#Compare FTIR.1 vs CRDS.1 for delta_CH4 at "South background"
delta_CH4_S_A_blandplot <- bland_altman_plot(data = emission_reshaped,
                                             var_filter = "delta_CH4",
                                             analyzer_pair = c("FTIR.1", "CRDS.1"),
                                             location_filter = "South background")

#Compare FTIR.1 vs CRDS.1 for delta_CO2 at "North background"
delta_CO2_N_A_blandplot <- bland_altman_plot(data = emission_reshaped,
                                             var_filter = "delta_CO2",
                                             analyzer_pair = c("FTIR.1", "CRDS.1"),
                                             location_filter = "North background")

#Compare FTIR.1 vs CRDS.1 for delta_CO2 at "South background"
delta_CO2_S_A_blandplot <- bland_altman_plot(data = emission_reshaped,
                                             var_filter = "delta_CO2",
                                             analyzer_pair = c("FTIR.1", "CRDS.1"),
                                             location_filter = "South background")

#Compare FTIR.1 vs CRDS.1 for delta_NH3 at "North background"
delta_NH3_N_A_blandplot <- bland_altman_plot(data = emission_reshaped,
                                             var_filter = "delta_NH3",
                                             analyzer_pair = c("FTIR.1", "CRDS.1"),
                                             location_filter = "North background")

#Compare FTIR.1 vs CRDS.1 for delta_NH3 at "South background"
delta_NH3_S_A_blandplot <- bland_altman_plot(data = emission_reshaped,
                                             var_filter = "delta_NH3",
                                             analyzer_pair = c("FTIR.1", "CRDS.1"),
                                             location_filter = "South background")

#save plots
# Add some margins to each plot
plots_list_A <- list(
        delta_CH4_N_A_blandplot, delta_CO2_N_A_blandplot, delta_NH3_N_A_blandplot,
        delta_CH4_S_A_blandplot, delta_CO2_S_A_blandplot, delta_NH3_S_A_blandplot)

plots_list_A <- lapply(plots_list_A, function(p) p + theme(plot.margin = margin(10, 10, 10, 10)))

# Combine with 3 columns × 2 rows
plots_A <- wrap_plots(plots_list_A, ncol = 3, nrow = 2)

ggsave("BlandAltman_AnalyzerA.pdf", plot = plots_A, width = 12, height = 6, units = "in", dpi = 300)

#Compare FTIR.2 vs CRDS.2 for delta_CH4 at "North background"
delta_CH4_N_B_blandplot <- bland_altman_plot(data = emission_reshaped,
                                             var_filter = "delta_CH4",
                                             analyzer_pair = c("FTIR.2", "CRDS.2"),
                                             location_filter = "North background")

#Compare FTIR.2 vs CRDS.2 for delta_CH4 at "South background"
delta_CH4_S_B_blandplot <- bland_altman_plot(data = emission_reshaped,
                                             var_filter = "delta_CH4",
                                             analyzer_pair = c("FTIR.2", "CRDS.2"),
                                             location_filter = "South background")

#Compare FTIR.2 vs CRDS.2 for delta_CO2 at "North background"
delta_CO2_N_B_blandplot <- bland_altman_plot(data = emission_reshaped,
                                             var_filter = "delta_CO2",
                                             analyzer_pair = c("FTIR.2", "CRDS.2"),
                                             location_filter = "North background")

#Compare FTIR.2 vs CRDS.2 for delta_CO2 at "South background"
delta_CO2_S_B_blandplot <- bland_altman_plot(data = emission_reshaped,
                                             var_filter = "delta_CO2",
                                             analyzer_pair = c("FTIR.2", "CRDS.2"),
                                             location_filter = "South background")

#Compare FTIR.2 vs CRDS.2 for delta_NH3 at "North background"
delta_NH3_N_B_blandplot <- bland_altman_plot(data = emission_reshaped,
                                             var_filter = "delta_NH3",
                                             analyzer_pair = c("FTIR.2", "CRDS.2"),
                                             location_filter = "North background")

#Compare FTIR.2 vs CRDS.2 for delta_NH3 at "South background"
delta_NH3_S_B_blandplot <- bland_altman_plot(data = emission_reshaped,
                                             var_filter = "delta_NH3",
                                             analyzer_pair = c("FTIR.2", "CRDS.2"),
                                             location_filter = "South background")

# Save plots
# Add some margins to each plot
plots_list_B <- list(
        delta_CH4_N_B_blandplot, delta_CO2_N_B_blandplot, delta_NH3_N_B_blandplot,
        delta_CH4_S_B_blandplot, delta_CO2_S_B_blandplot, delta_NH3_S_B_blandplot)

plots_list_B <- lapply(plots_list_B, function(p) p + theme(plot.margin = margin(10, 10, 10, 10)))

# Combine with 3 columns × 2 rows
plots_B <- wrap_plots(plots_list_B, ncol = 3, nrow = 2)

ggsave("BlandAltman_AnalyzerB.pdf", plot = plots_B, width = 12, height = 6, units = "in", dpi = 300)
######## Correlograms #######
# North background
d_CO2_N_corrgram <- emicorrgram(
        emission_reshaped,
        target_variable = "delta_CO2", 
        locations = "North background")

d_CH4_N_corrgram <- emicorrgram(
        emission_reshaped,
        target_variable = "delta_CH4",
        locations = "North background")

d_NH3_N_corrgram <- emicorrgram(
        emission_reshaped,
        target_variable = "delta_NH3",
        locations = "North background")

q_N_corrgram <- emicorrgram(
        emission_reshaped,
        target_variable = "Q_vent", 
        locations = "North background")

e_CH4_N_corrgram <- emicorrgram(
        emission_reshaped,
        target_variable = "e_CH4_ghLU", 
        locations = "North background")

e_NH3_N_corrgram <- emicorrgram(
        emission_reshaped,
        target_variable = "e_NH3_ghLU", 
        locations = "North background")

# South background
d_CO2_S_corrgram <- emicorrgram(
        emission_reshaped,
        target_variable = "delta_CO2", 
        locations = "South background")

d_CH4_S_corrgram <- emicorrgram(
        emission_reshaped,
        target_variable = "delta_CH4",
        locations = "South background")

d_NH3_S_corrgram <- emicorrgram(
        emission_reshaped,
        target_variable = "delta_NH3",
        locations = "South background")

q_S_corrgram <- emicorrgram(
        emission_reshaped,
        target_variable = "Q_vent", 
        locations = "South background")

e_CH4_S_corrgram <- emicorrgram(
        emission_reshaped,
        target_variable = "e_CH4_ghLU", 
        locations = "South background")

e_NH3_S_corrgram <- emicorrgram(
        emission_reshaped,
        target_variable = "e_NH3_ghLU", 
        locations = "South background")

# Combine all plots in a named list 
all_plots <- list(
        q_N_corrgram      = q_N_corrgram,
        q_S_corrgram      = q_S_corrgram,
        e_CH4_N_corrgram  = e_CH4_N_corrgram,
        e_CH4_S_corrgram  = e_CH4_S_corrgram,
        e_NH3_N_corrgram  = e_NH3_N_corrgram,
        e_NH3_S_corrgram  = e_NH3_S_corrgram
)

# Define custom sizes (width, height) for each plot
plot_sizes <- list(
        q_N_corrgram       = c(5, 6),
        q_S_corrgram       = c(5, 6),
        e_CH4_N_corrgram   = c(5, 6),
        e_CH4_S_corrgram   = c(5, 6),
        e_NH3_N_corrgram   = c(5, 6),
        e_NH3_S_corrgram   = c(5, 6)
)

# Loop through and save each plot with its custom size
for (plot_name in names(all_plots)) {
        size <- plot_sizes[[plot_name]]
        
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = all_plots[[plot_name]],
                width = size[1],
                height = size[2],
                dpi = 300
        )
}

######## Statistical Summary #########
emission_day_stat <- stat_error_table(emission_reshaped,
                                      time.group = "day",
                                      analyzer.level = c("FTIR.1", "FTIR.2", "FTIR.3", "FTIR.4",
                                                         "CRDS.1", "CRDS.2", "CRDS.3", "baseline")) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2)))

emission_hour_stat <- stat_error_table(emission_reshaped,
                                  time.group = "hour",
                                  analyzer.level = c("FTIR.1", "FTIR.2", "FTIR.3", "FTIR.4",
                                                     "CRDS.1", "CRDS.2", "CRDS.3", "baseline")) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2)))

# Delta gases daily summary
delta_day_summary <- emission_day_stat %>%
        filter(var %in% c("delta_CO2", "delta_CH4", "delta_NH3")) %>%
        group_by(day, analyzer, location, var) %>%
        summarise(mean_value = mean(mean_value, na.rm = TRUE),
                  sd         = mean(sd, na.rm = TRUE),
                  cv         = mean(cv, na.rm = TRUE),
                  mean_pct_err = mean(pct_err, na.rm = TRUE),
                  .groups = "drop")

# Q ventilation daily summary
q_day_summary <- emission_day_stat %>%
        filter(var == "Q_vent") %>%
        group_by(day, analyzer, location, var) %>%
        summarise(mean_value = mean(mean_value, na.rm = TRUE),
                  sd         = mean(sd, na.rm = TRUE),
                  cv         = mean(cv, na.rm = TRUE),
                  mean_pct_err = mean(pct_err, na.rm = TRUE),
                  .groups = "drop")

# CH4 emission daily summary
e_CH4_day_summary <- emission_day_stat %>%
        filter(var == "e_CH4_ghLU") %>%
        group_by(day, analyzer, location, var) %>%
        summarise(mean_value = mean(mean_value, na.rm = TRUE),
                  sd         = mean(sd, na.rm = TRUE),
                  cv         = mean(cv, na.rm = TRUE),
                  mean_pct_err = mean(pct_err, na.rm = TRUE),
                  .groups = "drop")

# NH3 emission daily summary
e_NH3_day_summary <- emission_day_stat %>%
        filter(analyzer == "baseline",
               var == "e_NH3_ghLU") %>%
        group_by(day, analyzer, location, var) %>%
        summarise(mean_value = mean(mean_value, na.rm = TRUE),
                  sd         = mean(sd, na.rm = TRUE),
                  cv         = mean(cv, na.rm = TRUE),
                  mean_pct_err = mean(pct_err, na.rm = TRUE),
                  .groups = "drop")


# Delta gases daily summary
delta_hour_summary <- emission_hour_stat %>%
        filter(var %in% c("delta_CO2", "delta_CH4", "delta_NH3")) %>%
        group_by(hour, analyzer, location, var) %>%
        summarise(mean_value = mean(mean_value, na.rm = TRUE),
                  sd         = mean(sd, na.rm = TRUE),
                  cv         = mean(cv, na.rm = TRUE),
                  mean_pct_err = mean(pct_err, na.rm = TRUE),
                  .groups = "drop")

# Q ventilation daily summary
q_hour_summary <- emission_hour_stat %>%
        filter(var == "Q_vent") %>%
        group_by(hour, analyzer, location, var) %>%
        summarise(mean_value = mean(mean_value, na.rm = TRUE),
                  sd         = mean(sd, na.rm = TRUE),
                  cv         = mean(cv, na.rm = TRUE),
                  mean_pct_err = mean(pct_err, na.rm = TRUE),
                  .groups = "drop")

# CH4 emission daily summary
e_CH4_hour_summary <- emission_hour_stat %>%
        filter(var == "e_CH4_ghLU") %>%
        group_by(hour, analyzer, location, var) %>%
        summarise(mean_value = mean(mean_value, na.rm = TRUE),
                  sd         = mean(sd, na.rm = TRUE),
                  cv         = mean(cv, na.rm = TRUE),
                  mean_pct_err = mean(pct_err, na.rm = TRUE),
                  .groups = "drop")

# NH3 emission daily summary
e_NH3_hour_summary <- emission_hour_stat %>%
        filter(analyzer == "baseline",
               var == "e_NH3_ghLU") %>%
        group_by(hour, analyzer, location, var) %>%
        summarise(mean_value = mean(mean_value, na.rm = TRUE),
                  sd         = mean(sd, na.rm = TRUE),
                  cv         = mean(cv, na.rm = TRUE),
                  mean_pct_err = mean(pct_err, na.rm = TRUE),
                  .groups = "drop")

# Save Day CSVs
readr::write_excel_csv(delta_day_summary, "delta_day_summary.csv")
readr::write_excel_csv(q_day_summary, "q_day_summary.csv")
readr::write_excel_csv(e_CH4_day_summary, "e_CH4_day_summary.csv")
readr::write_excel_csv(e_NH3_day_summary, "e_NH3_day_summary.csv")

# Save Hour CSVs
readr::write_excel_csv(delta_hour_summary, "delta_hour_summary.csv")
readr::write_excel_csv(q_hour_summary, "q_hour_summary.csv")
readr::write_excel_csv(e_CH4_hour_summary, "e_CH4_hour_summary.csv")
readr::write_excel_csv(e_NH3_hour_summary, "e_NH3_hour_summary.csv")

# Only ICC 
delta_icc_day <- icc_table(emission_reshaped  %>% filter(analyzer != "baseline"),
                           time.group = "day",
                           vars = c("delta_CO2", "delta_CH4", "delta_NH3")) 

q_icc_day <- icc_table(emission_reshaped  %>% filter(analyzer != "baseline"),
                       time.group = "day",
                       vars = "Q_vent")

e_icc_day <- icc_table(emission_reshaped%>% filter(analyzer != c("FTIR.4", "baseline")),
                       time.group = "day", 
                       vars = c("e_CH4_ghLU", "e_NH3_ghLU"))

# Save Day CSVs
readr::write_excel_csv(delta_icc_day, "delta_icc_day.csv")
readr::write_excel_csv(q_icc_day, "q_icc_day.csv")
readr::write_excel_csv(e_icc_day, "e_icc_day.csv")

######## Heatmap #########
# Daily heatmaps
q_S_day_heatmap  <- emiheatmap(emission_day_stat,
                               vars = "Q_vent",
                               time.group = "day",
                               location.filter = "South background")

q_N_day_heatmap  <- emiheatmap(emission_day_stat,
                               vars = "Q_vent",
                               time.group = "day",
                               location.filter = "North background")

# Hourly heatmaps
q_N_hour_heatmap <- emiheatmap(emission_hour_stat,
                               vars = "Q_vent",
                               time.group = "hour",
                               location.filter = "North background")

q_S_hour_heatmap <- emiheatmap(emission_hour_stat,
                               vars = "Q_vent",
                               time.group = "hour",
                               location.filter = "South background")

# Combine all plots in a named list 
all_plots <- list(
        q_S_day_heatmap  = q_S_day_heatmap,
        q_N_day_heatmap  = q_N_day_heatmap,
        q_N_hour_heatmap = q_N_hour_heatmap,
        q_S_hour_heatmap = q_S_hour_heatmap
)

# Define custom sizes (width, height) for each plot
plot_sizes <- list(
        q_S_day_heatmap    = c(6, 4),
        q_N_day_heatmap    = c(6, 4),
        q_N_hour_heatmap   = c(12, 4),
        q_S_hour_heatmap   = c(12, 4)
)

# Loop through and save each plot with its custom size
for (plot_name in names(all_plots)) {
        size <- plot_sizes[[plot_name]]
        
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = all_plots[[plot_name]],
                width = size[1],
                height = size[2],
                dpi = 300
        )
}

######## Linear mix modelling ##########
colnames(emission_result)

#emission_result <- emission_result %>% filter( analyzer != "FTIR.4")

# Delta CH4 North
linear_mixed_model("delta_CH4_N", emission_result)

# Delta CO2 North
linear_mixed_model("delta_CO2_N", emission_result)

# Delta NH3 North
linear_mixed_model("delta_NH3_N", emission_result)

# Delta CH4 South
linear_mixed_model("delta_CH4_S", emission_result)

# Delta CO2 South
linear_mixed_model("delta_CO2_S", emission_result)

# Delta NH3 South
linear_mixed_model("delta_NH3_S", emission_result)

# Q_vent emission rate North
linear_mixed_model("Q_vent_N", emission_result)

# Q_vent emission rate South
linear_mixed_model("Q_vent_S", emission_result)

# NH3 emission rate North
linear_mixed_model("e_NH3_ghLU_N", emission_result)

# NH3 emission rate South
linear_mixed_model("e_NH3_ghLU_S", emission_result)

# CH4 emission rate North
linear_mixed_model("e_CH4_ghLU_N", emission_result)

# CH4 emission rate South
linear_mixed_model("e_CH4_ghLU_S", emission_result)

