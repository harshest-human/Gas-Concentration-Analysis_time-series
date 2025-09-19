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

######## Development of Statistics functions #######
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

# Development of function stat_table to calculate summary statistics and RE
stat_table <- function(df_long, analyzer.levels = NULL) {
        library(dplyr)
        library(DescTools)
        library(lubridate)
        
        df_long <- df_long %>%
                # Add day and hour columns
                mutate(
                        day  = factor(format(as.Date(DATE.TIME), "%Y.%m.%d")),
                        hour = format(floor_date(DATE.TIME, unit = "hour"), "%H:%M")
                )
        
        # Filter analyzers if requested
        if (!is.null(analyzer.levels)) {
                df_long <- df_long %>%
                        mutate(analyzer = factor(analyzer, levels = analyzer.levels)) %>%
                        filter(analyzer %in% analyzer.levels | analyzer == "baseline")
        }
        
        # Calculate consortium metrics (exclude baseline row)
        consortium_metrics <- df_long %>%
                filter(analyzer != "baseline") %>%
                group_by(DATE.TIME, day, hour, location, var) %>%
                summarise(
                        baseline_mean = mean(value, na.rm = TRUE),
                        baseline_sd   = sd(value, na.rm = TRUE),
                        baseline_cv   = abs(CoefVar(value, na.rm = TRUE)) * 100,
                        .groups = "drop"
                )
        
        # Join consortium metrics to all rows (including baseline)
        df_long <- df_long %>%
                left_join(consortium_metrics, by = c("DATE.TIME", "day", "hour", "location", "var")) %>%
                # Compute relative error vs consortium mean
                mutate(
                        re = if_else(analyzer == "baseline", 0, 100 * (value - baseline_mean)/baseline_mean)
                ) %>%
                arrange(DATE.TIME, location, var, analyzer) %>%
                mutate(across(where(is.numeric), ~ round(.x, 2)))
        
        return(df_long)
}

# Development of function stat_table to calculate RE
pct_err <- function(df) {
        df %>%
                group_by(DATE.TIME, location, var) %>%
                mutate(
                        baseline_value = value[analyzer == "baseline"],
                        pct_err = 100 * (value - baseline_value) / baseline_value
                ) %>%
                ungroup() %>%
                select(-baseline_value)
}

# Development of function HSD_table and Pairwise difference in mean values
pairwise_diff <- function(data, vars, locations, group_var = "analyzer") {
        library(dplyr)
        library(rstatix)
        library(tidyr)
        library(ggplot2)
        
        df_plot_all <- list()
        
        # Loop over vars × locations
        for(v in vars) {
                for(loc in locations) {
                        df_var <- data %>% filter(var == v, location == loc)
                        analyzers <- sort(unique(df_var[[group_var]]))
                        n <- length(analyzers)
                        
                        # Skip if not enough data
                        counts <- df_var %>%
                                group_by(.data[[group_var]]) %>%
                                summarise(n = n(), .groups = "drop")
                        if(any(counts$n < 2)) next
                        
                        # Tukey HSD
                        tuk <- tryCatch(
                                rstatix::tukey_hsd(df_var, formula = as.formula(paste("value ~", group_var))) %>%
                                        mutate(label = case_when(
                                                is.na(p.adj)    ~ "-",
                                                p.adj <= 0.001  ~ "***",
                                                p.adj <= 0.01   ~ "**",
                                                p.adj <= 0.05   ~ "*",
                                                p.adj > 0.05    ~ "ns"
                                        )),
                                error = function(e) data.frame(group1 = character(), group2 = character(), label = character())
                        )
                        
                        # Percent error (lower triangle)
                        pct_err_df <- expand.grid(analyzer1 = analyzers, analyzer2 = analyzers, stringsAsFactors = FALSE) %>%
                                rowwise() %>%
                                mutate(
                                        pct_err = if(analyzer1 == analyzer2) NA else {
                                                val1 <- mean(df_var$value[df_var[[group_var]] == analyzer1], na.rm = TRUE)
                                                val2 <- mean(df_var$value[df_var[[group_var]] == analyzer2], na.rm = TRUE)
                                                100 * (val1 - val2)/val2
                                        }
                                ) %>% ungroup()
                        
                        # Fill matrix
                        mat <- matrix(NA, nrow = n, ncol = n, dimnames = list(analyzers, analyzers))
                        
                        # Lower triangle -> percent error
                        for(i in 2:n){
                                for(j in 1:(i-1)){
                                        mat[i,j] <- paste0(round(pct_err_df$pct_err[pct_err_df$analyzer1 == analyzers[i] &
                                                                                            pct_err_df$analyzer2 == analyzers[j]], 1), "%")
                                }
                        }
                        
                        # Upper triangle -> Tukey labels
                        for(i in 1:nrow(tuk)){
                                a1 <- tuk$group1[i]
                                a2 <- tuk$group2[i]
                                if(a1 %in% analyzers & a2 %in% analyzers){
                                        row_idx <- match(a1, analyzers)
                                        col_idx <- match(a2, analyzers)
                                        mat[row_idx, col_idx] <- tuk$label[i]
                                }
                        }
                        
                        # Diagonal
                        diag(mat) <- "-"
                        
                        # Convert to long format
                        df_long <- as.data.frame(mat, stringsAsFactors = FALSE) %>%
                                tibble::rownames_to_column("analyzer1") %>%
                                pivot_longer(-analyzer1, names_to = "analyzer2", values_to = "label") %>%
                                mutate(var = v, location = loc)
                        
                        df_plot_all[[paste(v, loc)]] <- df_long
                }
        }
        
        # Combine all
        df_final <- bind_rows(df_plot_all)
        
        # Plot
        p <- ggplot(df_final, aes(x = analyzer2, y = analyzer1)) +
                geom_tile(fill = "white", color = "black") +
                geom_text(aes(label = label), size = 3) +
                scale_y_discrete(limits = rev) +
                facet_grid(var ~ location, scales = "free") +
                theme_minimal() +
                theme(
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        strip.text = element_text(face = "bold"),
                        axis.text.x = element_text(angle = 45, hjust = 1)
                )
        
        return(p)
}

# Development of a Linear Mixed-Effects Model to test for interactions
linear_mixed_model <- function(var_name, data) {
        library(lme4)
        library(broom.mixed)
        library(dplyr)
        library(stringr)
        
        # Ensure 'day' column exists and is a factor
        data <- data %>%
                mutate(day = factor(as.Date(DATE.TIME)))
        
        # Detect location from variable name to select correct predictors
        if (str_detect(var_name, "_in$")) {
                location <- "in"
        } else if (str_detect(var_name, "_N$")) {
                location <- "N"
        } else if (str_detect(var_name, "_S$")) {
                location <- "S"
        } else {
                stop("Cannot detect location from variable name. Use _in, _N, or _S.")
        }
        
        # Define predictor variables based on location
        temp_var <- ifelse(location == "in", "temp_in", "temp_N")
        rh_var   <- ifelse(location == "in", "RH_in", "RH_N")
        cows_var <- "n_dairycows_in"
        wd_var <- ifelse(location %in% c("N","S"), "wd_mst_S", NULL)
        ws_var <- ifelse(location %in% c("N","S"), "ws_mst_S", NULL)
        
        predictors <- c(temp_var, rh_var, wd_var, ws_var, cows_var)
        predictors <- predictors[predictors %in% colnames(data)]
        
        # --- KEY CHANGE: MODIFIED FORMULA ---
        # 'analyzer' is now a fixed effect, interacting (*) with the environmental predictors.
        # This allows the model to estimate a different slope for each predictor for each analyzer.
        # '(1 | day)' remains as a random intercept for day-to-day variability.
        formula <- as.formula(
                paste(var_name, "~ (", paste(predictors, collapse = " + "), ") * analyzer + (1 | day)")
        )
        
        # Fit the model
        model <- lmer(formula, data = data, REML = TRUE)
        
        # Tidy the fixed effects table (this now includes the key interaction terms)
        fixed_tbl <- broom.mixed::tidy(model, effects = "fixed")
        
        # Tidy the random effects for 'day'
        ranefs <- ranef(model)
        day_effects <- ranefs$day %>%
                mutate(day = rownames(.)) %>%
                rename(intercept_adj = `(Intercept)`) %>%
                mutate(percent_dev = 100 * intercept_adj / mean(fixed_tbl$estimate[fixed_tbl$term=="(Intercept)"]))
        
        # Print results
        cat("\nFormula used:\n")
        print(formula)
        cat("\nFixed Effects (including interaction terms):\n")
        print(fixed_tbl, n = 50) # Print more rows to see all interactions
        cat("\nRandom Effects - Day:\n")
        print(day_effects)
        
        # Return a list of results (plotting is removed)
        invisible(list(
                formula = formula,
                model_summary = summary(model),
                fixed_effects = fixed_tbl,
                random_effects_day = day_effects
        ))
}

######## Development of Plotting functions ########
# Boxplot Function
emiboxplot <- function(data, y = NULL, location_filter = NULL, plot_err = FALSE) {
        library(dplyr)
        library(ggplot2)
        library(scales)
        library(stringr)
        
        if ("var" %in% names(data) && !"variable" %in% names(data)) {
                data <- data %>% rename(variable = var)
        }
        
        if (!is.null(location_filter)) {
                data <- data %>% filter(location %in% location_filter)
        }
        
        if (!is.null(y)) {
                data <- data %>% filter(variable %in% y)
        }
        
        value_col <- if (plot_err) "pct_err" else "value"
        if (!value_col %in% names(data)) stop("Column '", value_col, "' not found in data.")
        
        # Remove outliers using IQR method
        data <- data %>%
                group_by(location, analyzer, variable, day) %>%
                mutate(
                        Q1 = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
                        Q3 = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
                        IQR = Q3 - Q1,
                        is_outlier = (.data[[value_col]] < (Q1 - 1.5 * IQR)) | (.data[[value_col]] > (Q3 + 1.5 * IQR))
                ) %>%
                filter(!is_outlier) %>%
                ungroup()
        
        # Facet labels
        facet_labels_value <- c(
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
        
        facet_labels_err <- c(
                "CO2_mgm3"     = "c[CO2]",
                "CH4_mgm3"     = "c[CH4]",
                "NH3_mgm3"     = "c[NH3]",
                "r_CH4/CO2"    = "c[CH4]/c[CO2]",
                "r_NH3/CO2"    = "c[NH3]/c[CO2]",
                "delta_CO2"    = "Delta*c[CO2]",
                "delta_CH4"    = "Delta*c[CH4]",
                "delta_NH3"    = "Delta*c[NH3]",
                "Q_vent"       = "Q",
                "e_CH4_gh"     = "e[CH4]",
                "e_NH3_gh"     = "e[NH3]",
                "e_CH4_ghLU"   = "e[CH4]",
                "e_NH3_ghLU"   = "e[NH3]",
                "temp"         = "Temperature",
                "RH"           = "Relative~Humidity",
                "ws_mst"       = "Wind~Speed~Mast",
                "wd_mst"       = "Wind~Direction~Mast",
                "wd_trv"       = "Wind~Direction~Traverse",
                "ws_trv"       = "Wind~Speed~Traverse",
                "n_dairycows"  = "'Number of Cows'"
        )
        
        facet_labels_expr <- if (plot_err) facet_labels_err else facet_labels_value
        
        data <- data %>%
                mutate(facet_label = factor(
                        facet_labels_expr[as.character(variable)],
                        levels = facet_labels_expr[y]
                ))
        
        analyzer_colors <- c(
                "FTIR.1" = "#1b9e77", "FTIR.2" = "#d95f02", "FTIR.3" = "#7570b3",
                "FTIR.4" = "#e7298a", "CRDS.1" = "#66a61e", "CRDS.2" = "#e6ab02",
                "CRDS.3" = "#a6761d", "baseline" = "darkgrey"
        )
        
        all_analyzers <- unique(data$analyzer)
        analyzer_colors_full <- setNames(rep("black", length(all_analyzers)), all_analyzers)
        analyzer_colors_full[names(analyzer_colors)] <- analyzer_colors
        
        # No x and y axis labels
        xlab_to_use <- NULL
        ylab_to_use <- NULL
        
        p <- ggplot(data, aes(x = factor(day), y = .data[[value_col]], fill = analyzer)) +
                geom_boxplot(outlier.shape = NA) +
                scale_fill_manual(values = analyzer_colors_full) +
                facet_grid(facet_label ~ location, scales = "free_y", switch = "y",
                           labeller = labeller(
                                   facet_label = label_parsed,
                                   location = label_value
                           )) +
                scale_y_continuous(
                        breaks = scales::pretty_breaks(n = 6),
                        labels = scales::label_number(accuracy = 1, big.mark = "")
                ) +
                labs(x = xlab_to_use, y = ylab_to_use) +
                theme_classic() +
                theme(
                        text = element_text(size = 12),
                        axis.text = element_text(size = 12),
                        axis.title = element_text(size = 12),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                        strip.text.y.left = element_text(size = 12, vjust = 0.5),
                        strip.text.x = element_text(size = 12),
                        panel.border = element_rect(color = "black", fill = NA),
                        legend.position = "bottom",
                        legend.title = element_blank(),
                        plot.title = element_text(hjust = 0.5)
                ) +
                guides(fill = guide_legend(nrow = 1))
        
        print(p)
        return(p)
}

# Trend Plot Function 
emitrendplot <- function(data, y = NULL, location_filter = NULL, plot_err = FALSE, x = "DATE.TIME") {
        library(dplyr)
        library(ggplot2)
        library(scales)
        library(stringr)
        
        # Standardize variable name
        if ("var" %in% names(data) && !"variable" %in% names(data)) {
                data <- data %>% rename(variable = var)
        }
        
        # Filter by location and variable
        if (!is.null(location_filter)) {
                data <- data %>% filter(location %in% location_filter)
        }
        if (!is.null(y)) {
                data <- data %>% filter(variable %in% y)
        }
        
        # Choose value column
        value_col <- if (plot_err && "pct_err" %in% names(data)) "pct_err" else "value"
        if (!value_col %in% names(data)) {
                stop("Column '", value_col, "' not found in data.")
        }
        
        # Compute summary (mean ± SD)
        summary_data <- data %>%
                group_by(.data[[x]], analyzer, location, variable) %>%
                summarise(
                        mean_val = mean(.data[[value_col]], na.rm = TRUE),
                        sd_val = sd(.data[[value_col]], na.rm = TRUE),
                        .groups = "drop"
                ) %>%
                filter(!is.na(mean_val))
        
        if (nrow(summary_data) == 0) stop("No data available for the selected variables and plot_err setting.")
        
        # Facet labels
        facet_labels_value <- c(
                "CO2_mgm3"    = "c[CO2]~'(mg '*m^-3*')'",
                "CH4_mgm3"    = "c[CH4]~'(mg '*m^-3*')'",
                "NH3_mgm3"    = "c[NH3]~'(mg '*m^-3*')'",
                "r_CH4/CO2"   = "c[CH4]/c[CO2]~'('*'%'*')'",
                "r_NH3/CO2"   = "c[NH3]/c[CO2]~'('*'%'*')'",
                "delta_CO2"   = "Delta*c[CO2]~'(mg '*m^-3*')'",
                "delta_CH4"   = "Delta*c[CH4]~'(mg '*m^-3*')'",
                "delta_NH3"   = "Delta*c[NH3]~'(mg '*m^-3*')'",
                "Q_vent"      = "Q~'('*m^3~h^-1~LU^-1*')'",
                "e_CH4_gh"    = "e[CH4]~'(g '*h^-1*')'",
                "e_NH3_gh"    = "e[NH3]~'(g '*h^-1*')'",
                "e_CH4_ghLU"  = "e[CH4]~'(g '*h^-1~LU^-1*')'",
                "e_NH3_ghLU"  = "e[NH3]~'(g '*h^-1~LU^-1*')'",
                "temp"        = "Temperature~(degree*C)",
                "RH"          = "Relative~Humidity~('%')",
                "ws_mst"      = "Wind~Speed~Mast~(m~s^-1)",
                "wd_mst"      = "Wind~Direction~Mast~(degree)",
                "wd_trv"      = "Wind~Direction~Traverse~(degree)",
                "ws_trv"      = "Wind~Speed~Traverse~(m~s^-1)",
                "n_dairycows" = "'Number of Cows'"
        )
        
        facet_labels_err <- c(
                "CO2_mgm3"    = "c[CO2]",
                "CH4_mgm3"    = "c[CH4]",
                "NH3_mgm3"    = "c[NH3]",
                "r_CH4/CO2"   = "c[CH4]/c[CO2]",
                "r_NH3/CO2"   = "c[NH3]/c[CO2]",
                "delta_CO2"   = "Delta*c[CO2]",
                "delta_CH4"   = "Delta*c[CH4]",
                "delta_NH3"   = "Delta*c[NH3]",
                "Q_vent"      = "Q",
                "e_CH4_gh"    = "e[CH4]",
                "e_NH3_gh"    = "e[NH3]",
                "e_CH4_ghLU"  = "e[CH4]",
                "e_NH3_ghLU"  = "e[NH3]",
                "temp"        = "Temperature",
                "RH"          = "Relative~Humidity",
                "ws_mst"      = "Wind~Speed~Mast",
                "wd_mst"      = "Wind~Direction~Mast",
                "wd_trv"      = "Wind~Direction~Traverse",
                "ws_trv"      = "Wind~Speed~Traverse",
                "n_dairycows" = "'Number of Cows'"
        )
        
        facet_labels_expr <- if (plot_err) facet_labels_err else facet_labels_value
        
        summary_data <- summary_data %>%
                mutate(facet_label = factor(
                        facet_labels_expr[as.character(variable)],
                        levels = facet_labels_expr[y]
                ))
        
        # Analyzer aesthetics
        analyzer_colors <- c(
                "FTIR.1"  = "#1b9e77", "FTIR.2"  = "#d95f02", "FTIR.3"  = "#7570b3",
                "FTIR.4"  = "#e7298a", "CRDS.1"  = "#66a61e", "CRDS.2"  = "#e6ab02",
                "CRDS.3"  = "#a6761d", "baseline" = "black"
        )
        analyzer_shapes <- c(
                "FTIR.1"  = 0, "FTIR.2"  = 1, "FTIR.3"  = 2, "FTIR.4" = 5,
                "CRDS.1"  = 15, "CRDS.2" = 19, "CRDS.3" = 17, "baseline" = 4
        )
        all_analyzers <- unique(summary_data$analyzer)
        analyzer_colors_full <- setNames(rep("black", length(all_analyzers)), all_analyzers)
        analyzer_shapes_full <- setNames(rep(16, length(all_analyzers)), all_analyzers)
        analyzer_colors_full[names(analyzer_colors)] <- analyzer_colors
        analyzer_shapes_full[names(analyzer_shapes)] <- analyzer_shapes
        
        # Base plot
        p <- ggplot(summary_data, aes(x = .data[[x]], y = mean_val,
                                      color = analyzer, shape = analyzer, group = analyzer))
        
        if (x %in% c("day", "hour")) {
                p <- p +
                        geom_point(size = 2, position = position_dodge(width = 0.5)) +
                        geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                                      width = 0.4, position = position_dodge(width = 0.5))
        } else {
                p <- p +
                        geom_line(size = 0.5, alpha = 0.6, na.rm = TRUE) +
                        geom_point(size = 2, na.rm = TRUE) +
                        geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                                      width = 0.4, na.rm = TRUE)
        }
        
        # Add styling/scales
        p <- p +
                scale_color_manual(values = analyzer_colors_full) +
                scale_shape_manual(values = analyzer_shapes_full) +
                facet_grid(facet_label ~ location, scales = "free_y", switch = "y",
                           labeller = labeller(
                                   facet_label = label_parsed,
                                   location = label_value
                           )) +
                scale_y_continuous(
                        breaks = scales::pretty_breaks(n = 6),
                        labels = scales::label_number(accuracy = 1, big.mark = "")
                ) +
                labs(x = NULL, y = NULL) +
                theme_classic() +
                theme(
                        text = element_text(size = 12),
                        axis.text = element_text(size = 12),
                        axis.title = element_text(size = 12),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                        strip.text.y.left = element_text(size = 12, vjust = 0.5),
                        strip.text.x = element_text(size = 12),
                        panel.border = element_rect(color = "black", fill = NA),
                        legend.position = "bottom",
                        legend.title = element_blank(),
                        plot.title = element_text(hjust = 0.5)
                ) +
                guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
        
        # Time axis formatting
        if (x == "DATE.TIME") {
                x_breaks <- seq(from = min(summary_data$DATE.TIME, na.rm = TRUE),
                                to = max(summary_data$DATE.TIME, na.rm = TRUE),
                                by = "6 hours")
                p <- p + scale_x_datetime(breaks = x_breaks, date_labels = "%Y-%m-%d %H:%M")
        }
        
        print(p)
        return(p)
}

# Development of function for correlogram
emicorrgram <- function(data, target_variables, locations = NULL) {
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(Hmisc)
        library(purrr)
        
        # Parsed labels for facets
        facet_labels_expr <- c(
                "CO2_mgm3"     = "c[CO2]",
                "CH4_mgm3"     = "c[CH4]",
                "NH3_mgm3"     = "c[NH3]",
                "delta_CO2"    = "Delta*c[CO2]",
                "delta_CH4"    = "Delta*c[CH4]",
                "delta_NH3"    = "Delta*c[NH3]",
                "Q_vent"       = "Q",
                "e_CH4_gh"     = "e[CH4]",
                "e_NH3_gh"     = "e[NH3]",
                "e_CH4_ghLU"   = "e[CH4]",
                "e_NH3_ghLU"   = "e[NH3]"
        )
        
        # Helper function to compute correlation for one variable & location
        compute_corr <- function(var_sel, loc_sel) {
                filtered <- data %>%
                        filter(var == var_sel) %>%
                        select(DATE.TIME, location, analyzer, value)
                
                if (!is.null(loc_sel)) {
                        filtered <- filtered %>% filter(location %in% loc_sel)
                }
                
                filtered <- filtered %>%
                        group_by(DATE.TIME, location, analyzer) %>%
                        summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
                
                pivoted <- filtered %>%
                        pivot_wider(names_from = analyzer, values_from = value,
                                    values_fn = function(x, ...) mean(x, na.rm = TRUE))
                
                subdata_num <- pivoted %>%
                        select(-DATE.TIME, -location) %>%
                        mutate(across(everything(), as.numeric))
                
                if (ncol(subdata_num) < 2) return(NULL)
                
                cor_res <- Hmisc::rcorr(as.matrix(subdata_num), type = "pearson")
                cor_mat <- cor_res$r
                p_mat <- cor_res$P
                
                cor_long <- expand.grid(
                        Var1 = colnames(cor_mat),
                        Var2 = colnames(cor_mat),
                        stringsAsFactors = FALSE
                ) %>%
                        mutate(
                                correlation = as.vector(cor_mat),
                                pvalue = as.vector(p_mat),
                                target_variable = var_sel,
                                location = loc_sel
                        ) %>%
                        filter(as.numeric(factor(Var1)) > as.numeric(factor(Var2)))
                
                return(cor_long)
        }
        
        # Loop over variables and locations
        corr_df <- map_df(target_variables, function(tv) {
                map_df(locations, function(loc) compute_corr(tv, loc))
        })
        
        # Add parsed facet labels
        corr_df <- corr_df %>%
                mutate(
                        facet_label = factor(
                                facet_labels_expr[target_variable],
                                levels = facet_labels_expr[target_variables]
                        )
                )
        
        # Plot
        p <- ggplot(corr_df, aes(x = Var1, y = Var2, fill = correlation)) +
                geom_tile(color = "white") +
                geom_text(aes(label = paste0(
                        round(correlation, 2),
                        ifelse(pvalue <= 0.001, "\n***",
                               ifelse(pvalue <= 0.01, "\n**",
                                      ifelse(pvalue <= 0.05, "\n*", "\nns")))
                )), size = 3, color = "white") +
                scale_y_discrete(position = "right") +
                scale_fill_gradientn(
                        colors = c("darkred","red","orange2","gold2","yellow",
                                   "greenyellow","green1","green3","darkgreen"),
                        limits = c(0,1),
                        name = "PCC"
                ) +
                facet_grid(location ~ facet_label, scales = "free", space = "free", switch = "y",
                           labeller = labeller(facet_label = label_parsed)) +
                theme_classic() +
                theme(
                        axis.title = element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
                        legend.position = "bottom",
                        strip.text = element_text(size = 12),
                        panel.border = element_rect(colour = "black", fill = NA)
                )
        
        print(p)
        return(p)
}

# Development of function cv emiheatmap
emiheatmap <- function(data, vars, locations = c("North background", "South background")) {
        library(ggplot2)
        library(dplyr)
        library(lubridate)
        library(tidyr)
        
        # --- 1. Filter data for variable and locations ---
        df_filtered <- data %>%
                filter(var %in% vars, location %in% locations)
        
        if (nrow(df_filtered) == 0) {
                warning("No data found for specified variable/location.")
                return(ggplot() + theme_void() + ggtitle("No data available"))
        }
        
        # --- 2. Extract day and hour ---
        df_filtered <- df_filtered %>%
                mutate(
                        day  = factor(format(as.Date(DATE.TIME), "%Y.%m.%d"), 
                                      levels = sort(unique(format(as.Date(DATE.TIME), "%Y.%m.%d")))),
                        hour = factor(format(floor_date(DATE.TIME, unit = "hour"), "%H:%M"),
                                      levels = sprintf("%02d:00", 0:23),   # ensures 00 to 23
                                      ordered = TRUE)
                )
        
        # --- 3. Prepare facet labels ---
        title_labels_expr <- c(
                "delta_CO2"    = "Delta*c[CO2]",
                "delta_CH4"    = "Delta*c[CH4]",
                "delta_NH3"    = "Delta*c[NH3]",
                "Q_vent"       = "Q",
                "e_CH4_gh"     = "e[CH4]",
                "e_NH3_gh"     = "e[NH3]",
                "e_CH4_ghLU"   = "e[CH4]",
                "e_NH3_ghLU"   = "e[NH3]"
        )
        
        df_prepared <- df_filtered %>%
                mutate(
                        facet_label_str = ifelse(var %in% names(title_labels_expr),
                                                 paste0(title_labels_expr[var], " ~ '(", location, ")'"),
                                                 paste0(var, " ~ '(", location, ")'"))
                )
        
        # --- Ordered factor for facet wrap ---
        facet_order <- paste0(title_labels_expr[vars], " ~ '(", locations, ")'")
        df_prepared$facet_label <- factor(df_prepared$facet_label_str, levels = facet_order)
        
        # --- 4. Create heatmap ---
        p <- ggplot(df_prepared, aes(x = day, y = hour, fill = baseline_cv)) +
                geom_tile(color = "white") +
                geom_text(aes(label = ifelse(is.na(baseline_cv), "", round(baseline_cv, 1))),
                          size = 3, color = "black") +
                scale_fill_gradientn(
                        colors = c("darkgreen", "green4", "green2", "yellow", "orange", "red", "darkred"),
                        limits = c(0, 100),
                        breaks = seq(0, 100, by = 10),
                        name = "CV (%)",
                        na.value = "grey90"
                ) +
                facet_wrap(~ facet_label, scales = "free", ncol = length(unique(df_prepared$facet_label)),
                           labeller = label_parsed) +
                labs(x = NULL, y = NULL) +
                theme_classic() +
                theme(
                        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  
                        axis.text.y = element_text(size = 10),                  
                        strip.text = element_text(size = 10),
                        strip.background = element_rect(fill = "white", color = "black"),
                        panel.border = element_rect(color = "black", fill = NA),
                        legend.position = "bottom",
                        legend.key.width = unit(1.5, "cm"),
                        plot.margin = margin(t = 5, r = 10, b = 5, l = 25, unit = "pt"),
                        panel.spacing = unit(1.5, "cm")
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

# Remove constants and formla inputs
emission_result <- emission_result %>% select(-hour, -m_weight, -p_pregnancy_day, 
                                              -Y1_milk_prod, -a, -h_min,
                                              -phi, -t_factor, -phi_T_cor, -A_cor,
                                              -hpu_T_A_cor_per_cow, -PCO2)

# Write csv 
write_excel_csv(emission_result, "20250408-15_ringversuche_emission_result.csv")

######## Reshape Data #########
# Concentrations dataset
concentration_result <- indirect.CO2.balance(input_combined) %>%
        select("DATE.TIME", "analyzer", 
               "CO2_mgm3_in", "CH4_mgm3_in", "NH3_mgm3_in",
               "CO2_mgm3_S", "CH4_mgm3_S", "NH3_mgm3_S",
               "CO2_mgm3_N", "CH4_mgm3_N", "NH3_mgm3_N",
               "delta_CO2_N", "delta_CH4_N", "delta_NH3_N",
               "delta_CO2_S", "delta_CH4_S", "delta_NH3_S")

concentration_reshaped <- reshaper(concentration_result) 

concentration_reshaped <- pct_err(concentration_reshaped) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2)))

write_excel_csv(concentration_reshaped, "20250408-15_ringversuche_concentration_reshaped.csv")

#Emissions Dataset
emission_reshaped <-  reshaper(emission_result %>% filter(analyzer != "FTIR.4")) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2)))

# Write csv
write_excel_csv(emission_reshaped, "20250408-15_ringversuche_emission_reshaped.csv")

######## HSD and RE Matrix Plots ########
vars <- c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3")
locations <- c("Barn inside", "North background", "South background")

p <- pairwise_diff(concentration_reshaped, vars = vars, locations = locations)
print(p)

######## Absolute Concentration Plots ########
c_boxplot <- emiboxplot(
        data = concentration_reshaped,
        y = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3"),
        plot_err = FALSE)

c_trend_plot <- emitrendplot(
        data = concentration_reshaped,
        y = c("CO2_mgm3", "CH4_mgm3", "NH3_mgm3"),
        plot_err = FALSE)

all_plots_conc <- list(
        c_boxplot   = c_boxplot,
        c_trend_plot = c_trend_plot)

plot_sizes_conc <- list(
        c_boxplot   = c(14, 12),
        c_trend_plot = c(16, 12))

for (plot_name in names(all_plots_conc)) {
        size <- plot_sizes_conc[[plot_name]]
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = all_plots_conc[[plot_name]],
                width = size[1],
                height = size[2],
                dpi = 300)
}

######## Delta Concentration Plots ########
d_boxplot <- emiboxplot(
        data = concentration_reshaped,
        y = c("delta_CO2", "delta_CH4", "delta_NH3"),
        plot_err = FALSE)

d_trend_plot <- emitrendplot(
        data = concentration_reshaped,
        y = c("delta_CO2", "delta_CH4", "delta_NH3"),
        plot_err = FALSE)

all_plots_delta <- list(
        d_boxplot   = d_boxplot,
        d_trend_plot = d_trend_plot)

plot_sizes_delta <- list(
        d_boxplot   = c(12, 8),
        d_trend_plot = c(12, 8))

for (plot_name in names(all_plots_delta)) {
        size <- plot_sizes_delta[[plot_name]]
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = all_plots_delta[[plot_name]],
                width = size[1],
                height = size[2],
                dpi = 300)
}

######## Ventilation and Emission Rate Plots ########
q_e_boxplot <- emiboxplot(data = emission_reshaped,
                          y = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"),
                          plot_err = FALSE)

q_e_trend_plot <- emitrendplot(data = emission_reshaped,
                               y = c("Q_vent", "e_CH4_ghLU", "e_NH3_ghLU"),
                               plot_err = FALSE)


all_plots_vent <- list(
        q_e_boxplot   = q_e_boxplot,
        q_e_trend_plot = q_e_trend_plot
)

plot_sizes_vent <- list(
        q_e_boxplot   = c(12, 8),
        q_e_trend_plot = c(12, 8))

for (plot_name in names(all_plots_vent)) {
        size <- plot_sizes_vent[[plot_name]]
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = all_plots_vent[[plot_name]],
                width = size[1],
                height = size[2],
                dpi = 300
        )
}

######## Bland-Altman plot ############
# Lab A
e_CH4_N_A_blandplot <- bland_altman_plot(data = emission_reshaped,
                                         var_filter = "e_CH4_ghLU",
                                         analyzer_pair = c("FTIR.1", "CRDS.1"),
                                         location_filter = "North background")

e_CH4_S_A_blandplot <- bland_altman_plot(data = emission_reshaped,
                                         var_filter = "e_CH4_ghLU",
                                         analyzer_pair = c("FTIR.1", "CRDS.1"),
                                         location_filter = "South background")

e_NH3_N_A_blandplot <- bland_altman_plot(data = emission_reshaped,
                                         var_filter = "e_NH3_ghLU",
                                         analyzer_pair = c("FTIR.1", "CRDS.1"),
                                         location_filter = "North background")

e_NH3_S_A_blandplot <- bland_altman_plot(data = emission_reshaped,
                                         var_filter = "e_NH3_ghLU",
                                         analyzer_pair = c("FTIR.1", "CRDS.1"),
                                         location_filter = "South background")

q_N_A_blandplot <- bland_altman_plot(data = emission_reshaped,
                                     var_filter = "Q_vent",
                                     analyzer_pair = c("FTIR.1", "CRDS.1"),
                                     location_filter = "North background")

q_S_A_blandplot <- bland_altman_plot(data = emission_reshaped,
                                     var_filter = "Q_vent",
                                     analyzer_pair = c("FTIR.1", "CRDS.1"),
                                     location_filter = "South background")

# Save multiple plots as one pdf
plots_list_A <- list(
        e_CH4_N_A_blandplot, e_NH3_N_A_blandplot,
        e_CH4_S_A_blandplot, e_NH3_S_A_blandplot)

plots_list_A <- lapply(plots_list_A, function(p) p + theme(plot.margin = margin(10, 10, 10, 10)))

# Arrange in 1 row × 4 columns
plots_A <- wrap_plots(plots_list_A, ncol = 4, nrow = 1)

ggsave("e_BlandAltman_AnalyzerA.pdf", plot = plots_A,
       width = 16, height = 4, units = "in", dpi = 300)

# Save multiple plots as one pdf
plots_list_A_Q <- list(
        q_N_A_blandplot, q_S_A_blandplot)

plots_list_A_Q <- lapply(plots_list_A_Q, function(p) p + theme(plot.margin = margin(10, 10, 10, 10)))

plots_A_Q <- wrap_plots(plots_list_A_Q, ncol = 2, nrow = 1)

ggsave("q_BlandAltman_AnalyzerA.pdf", plot = plots_A_Q,
       width = 10, height = 4, units = "in", dpi = 300)


# Lab B
e_CH4_N_B_blandplot <- bland_altman_plot(data = emission_reshaped,
                                         var_filter = "e_CH4_ghLU",
                                         analyzer_pair = c("FTIR.2", "CRDS.2"),
                                         location_filter = "North background")

e_CH4_S_B_blandplot <- bland_altman_plot(data = emission_reshaped,
                                         var_filter = "e_CH4_ghLU",
                                         analyzer_pair = c("FTIR.2", "CRDS.2"),
                                         location_filter = "South background")

e_NH3_N_B_blandplot <- bland_altman_plot(data = emission_reshaped,
                                         var_filter = "e_NH3_ghLU",
                                         analyzer_pair = c("FTIR.2", "CRDS.2"),
                                         location_filter = "North background")

e_NH3_S_B_blandplot <- bland_altman_plot(data = emission_reshaped,
                                         var_filter = "e_NH3_ghLU",
                                         analyzer_pair = c("FTIR.2", "CRDS.2"),
                                         location_filter = "South background")

q_N_B_blandplot <- bland_altman_plot(data = emission_reshaped,
                                     var_filter = "Q_vent",
                                     analyzer_pair = c("FTIR.2", "CRDS.2"),
                                     location_filter = "North background")

q_S_B_blandplot <- bland_altman_plot(data = emission_reshaped,
                                     var_filter = "Q_vent",
                                     analyzer_pair = c("FTIR.2", "CRDS.2"),
                                     location_filter = "South background")


# Save multiple plots as one pdf
plots_list_B <- list(
        e_CH4_N_B_blandplot, e_NH3_N_B_blandplot,
        e_CH4_S_B_blandplot, e_NH3_S_B_blandplot)

plots_list_B <- lapply(plots_list_B, function(p) p + theme(plot.margin = margin(10, 10, 10, 10)))

plots_B <- wrap_plots(plots_list_B, ncol = 4, nrow = 1)

ggsave("e_BlandAltman_AnalyzerB.pdf", plot = plots_B,
       width = 16, height = 4, units = "in", dpi = 300)

# Save multiple plots as one pdf
plots_list_B_Q <- list(
        q_N_B_blandplot, q_S_B_blandplot)

plots_list_B_Q <- lapply(plots_list_B_Q, function(p) p + theme(plot.margin = margin(10, 10, 10, 10)))

plots_B_Q <- wrap_plots(plots_list_B_Q, ncol = 2, nrow = 1)

ggsave("q_BlandAltman_AnalyzerB.pdf", plot = plots_B_Q,
       width = 10, height = 4, units = "in", dpi = 300)


######## Correlograms #######
q_e_corrgram <- emicorrgram(emission_reshaped,
                            target_variables = c("e_CH4_ghLU", "e_NH3_ghLU", "Q_vent"),
                            locations = c("North background", "South background"))

ggsave("q_e_corrgram.pdf", plot = q_e_corrgram,
       width = 10, height = 6, units = "in", dpi = 300)

######## Heatmap #########
d_CO2_heatmap <- emiheatmap(
        data = concentration_stat,
        var = "delta_CO2",
        locations = c("North background", "South background"))

d_CH4_heatmap <- emiheatmap(
        data = concentration_stat,
        var = "delta_CH4",
        locations = c("North background", "South background"))


d_NH3_heatmap <- emiheatmap(
        data = concentration_stat,
        var = "delta_NH3",
        locations = c("North background", "South background"))


q_heatmap <- emiheatmap(
        data = emission_stat,
        var = "Q_vent",
        locations = c("North background", "South background"))

e_CH4_heatmap <- emiheatmap(
        data = emission_stat,
        var = "e_CH4_ghLU",
        locations = c("North background", "South background"))

e_NH3_heatmap <- emiheatmap(
        data = emission_stat,
        var = "e_NH3_ghLU",
        locations = c("North background", "South background"))

# 1. Combine all plots in a named list
all_plots <- list(
        d_CO2_heatmap  = d_CO2_heatmap,
        d_CH4_heatmap  = d_CH4_heatmap,
        d_NH3_heatmap  = d_NH3_heatmap,
        q_heatmap      = q_heatmap,
        e_CH4_heatmap  = e_CH4_heatmap,
        e_NH3_heatmap  = e_NH3_heatmap
)

# 2. Define custom sizes
plot_sizes <- list(
        d_CO2_heatmap  = c(8, 8),
        d_CH4_heatmap  = c(8, 8),
        d_NH3_heatmap  = c(8, 8),
        q_heatmap      = c(8, 8),
        e_CH4_heatmap  = c(8, 8),
        e_NH3_heatmap  = c(8, 8)
)

# 3. Save all plots
for (plot_name in names(all_plots)) {
        size <- plot_sizes[[plot_name]]
        
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = all_plots[[plot_name]],
                width = size[1],
                height = size[2],
                dpi = 300
        )
        
        message("Saved: ", paste0(plot_name, ".pdf"))
}

######## Environment plots #######
emission_reshaped <-  emission_reshaped  %>%
        mutate(location = if_else(
                var %in% c("temp", "wd_mst", "ws_mst", "RH"),"Weather", location))

weather_trendplot <- emitrendplot(
        data = emission_reshaped,
        y = c("temp", "wd_mst", "ws_mst", "RH"))

ggsave("weather_trendplot.pdf", plot = weather_trendplot,
       width = 10, height = 12, units = "in", dpi = 300)

######## Linear mix modelling ##########
colnames(emission_result)

emission_result <- emission_result %>% filter( analyzer != "FTIR.4")

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


