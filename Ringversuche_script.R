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

######## Development of functions #######
# Development of indirect.CO2.balance function
indirect.CO2.balance <- function(df) {
        df %>%
                mutate(
                        # Constants and time component
                        P_CO2_term = 0.185,
                        a = 0.22,
                        h_min = 2.9,
                        hour = hour(DATE.TIME),
                        
                        # Correction factors and animal-level CO2 production
                        A_cor = 1 - a * 3 * sin((2 * pi / 24) * (hour + 6 - h_min)),
                        Phi_tot = 5.6 * (m_weight ^ 0.75) + 22 * Y1_milk_prod + 1.6e-5 * (p_pregnancy_day ^ 3),
                        Phi_T_cor = Phi_tot * (1 + 4e-5 * (20 - Temperature)^3),
                        hpu_T_A_cor_all_animal = ((Phi_T_cor / 1000) * A_cor) * n_animals,
                        n_LU = (n_animals * m_weight) / 500,
                        P_CO2_T_A_all_animal = hpu_T_A_cor_all_animal * P_CO2_term,
                        
                        # Concentration differences (mg m3)
                        delta_NH3_N = NH3_in - NH3_N,
                        delta_CH4_N = CH4_in - CH4_N,
                        delta_CO2_N = CO2_in - CO2_N,
                        
                        delta_NH3_S = NH3_in - NH3_S,
                        delta_CH4_S = CH4_in - CH4_S,
                        delta_CO2_S = CO2_in - CO2_S,
                        
                        # Ventilation rate (m³/h)
                        Q_Vent_rate_N = P_CO2_T_A_all_animal / ((delta_CO2_N) * 1e-6),
                        Q_Vent_rate_S = P_CO2_T_A_all_animal / ((delta_CO2_S) * 1e-6),
                        
                        # Emissions (g/h)
                        e_NH3_N = (delta_NH3_N * Q_Vent_rate_N) / 1000,
                        e_CH4_N = (delta_CH4_N * Q_Vent_rate_N) / 1000,
                        e_NH3_S = (delta_NH3_S * Q_Vent_rate_S) / 1000,
                        e_CH4_S = (delta_CH4_S * Q_Vent_rate_S) / 1000,
                        
                        # NH3 to CO2 ratios (%)
                        NHCO_in = (NH3_in / CO2_in) * 100,
                        NHCO_N  = (NH3_N  / CO2_N) * 100,
                        NHCO_S  = (NH3_S  / CO2_S) * 100,
                        
                        # NH3 to CH4 ratios (%)
                        NHCH_in = (NH3_in / CH4_in) * 100,
                        NHCH_N  = (NH3_N  / CH4_N) * 100,
                        NHCH_S  = (NH3_S  / CH4_S) * 100,
                        
                        # CH4 to CO2 ratios (%)
                        CHCO_in = (CH4_in / CO2_in) * 100,
                        CHCO_N  = (CH4_N  / CO2_N) * 100,
                        CHCO_S  = (CH4_S  / CO2_S) * 100
                )
}

# Development of function relative error calculator
emirror <- function(df, vars) {
        # Columns to keep besides vars
        meta_cols <- c("DATE.TIME", "hour", "day", "lab", "analyzer")
        
        # Select only meta_cols + vars from df
        df_selected <- df %>% select(all_of(c(meta_cols, vars)))
        
        # Compute relative errors only on vars
        df_result <- df_selected %>%
                mutate(across(
                        all_of(vars),
                        ~ (. - mean(., na.rm = TRUE)) / mean(., na.rm = TRUE) * 100,
                        .names = "Err_{.col}"
                ))%>%
                mutate(across(where(is.numeric) & !any_of(c("hour")), ~ round(.x, 2)))
        
        return(df_result)
}

# Development of function reshape parameters
reparam <- function(data) {
        library(dplyr)
        library(tidyr)
        library(stringr)
        
        meta_cols <- c("DATE.TIME", "day", "hour", "lab", "analyzer")
        
        # Patterns for gases and ventilation columns
        gases <- c("CO2", "CH4", "NH3", "NHCO", "NHCH", "CHCO")
        vent_cols <- c("Q_Vent_rate_N", "Q_Vent_rate_S")
        
        # Pattern to select relevant columns (including prefixed ones)
        pattern <- paste0(
                "^(", paste(c(
                        # gases with suffixes
                        paste0(gases, "_(in|N|S)$"),
                        # prefixed variables
                        paste0("delta_(", gases[1:3], ")_(N|S)$"),
                        paste0("e_(", gases[3:2], ")_(N|S)$"),  # e_NH3, e_CH4 emissions
                        paste0("Q_Vent_rate_(N|S)$"),
                        # errors for gases, delta, emissions, ventilation
                        paste0("Err_(", gases, ")_(in|N|S)$"),
                        paste0("Err_delta_(", gases[1:3], ")_(N|S)$"),
                        paste0("Err_e_(", gases[3:2], ")_(N|S)$"),
                        paste0("Err_Q_Vent_rate_(N|S)$")
                ), collapse = "|"), ")$"
        )
        
        # Select columns matching pattern
        cols_to_use <- grep(pattern, names(data), value = TRUE)
        
        all_vars <- c(meta_cols, cols_to_use)
        
        df_long <- data %>%
                select(all_of(all_vars)) %>%
                pivot_longer(
                        cols = -all_of(meta_cols),
                        names_to = "variable",
                        values_to = "value"
                ) %>%
                mutate(
                        value = round(value, 2),
                        location_code = str_extract(variable, "(?<=_)(in|N|S)(?=(_|$))"),
                        location = case_when(
                                location_code == "in" ~ "Ringline inside",
                                location_code == "N" ~ "North background",
                                location_code == "S" ~ "South background",
                                TRUE ~ NA_character_
                        ),
                        
                        var_type = case_when(
                                str_starts(variable, "e_") ~ "emission",
                                str_starts(variable, "delta_") ~ "delta",
                                str_starts(variable, "Err_") ~ "error",
                                str_detect(variable, paste0("^(", paste(gases[c(4,5,6)], collapse = "|"), ")_")) ~ "ratio",
                                str_starts(variable, "Q_Vent_rate") ~ "ventilation",
                                TRUE ~ "concentration"
                        ),
                        
                        # Clean suffix only: remove the location suffix _in, _N, _S for clarity
                        variable = str_remove(variable, "_(in|N|S)$"),
                        
                        # Replace ratio names with more readable forms
                        variable = case_when(
                                variable == "NHCO" ~ "NH3/CO2",
                                variable == "NHCH" ~ "NH3/CH4",
                                variable == "CHCO" ~ "CH4/CO2",
                                variable == "Err_NHCO" ~ "Err_NH3/CO2",
                                variable == "Err_NHCH" ~ "Err_NH3/CH4",
                                variable == "Err_CHCO" ~ "Err_CH4/CO2",
                                TRUE ~ variable
                        )
                ) %>%
                select(all_of(meta_cols), location, var_type, variable, value)
        
        return(df_long)
}

# Development of function stat_table
stat_table <- function(data, response_vars, group_vars, var_type_filter = NULL) {
        require(dplyr)
        require(DescTools)
        
        df <- data
        
        if (!is.null(var_type_filter)) {
                df <- df %>% filter(var_type %in% var_type_filter)
        }
        
        df %>%
                filter(variable %in% response_vars) %>%
                group_by(across(all_of(c("hour", group_vars))), variable) %>%
                summarise(
                        n = n(),
                        mean = round(mean(value, na.rm = TRUE), 2),
                        sd = round(sd(value, na.rm = TRUE), 2),
                        cv = round(DescTools::CoefVar(value, na.rm = TRUE) * 100, 2),
                        .groups = "drop"
                )
}

# Function to summarize, pivot wider, and round
hour_sum <- function(df) {
        
        library(dplyr)
        library(tidyr)
        
        df %>%
                group_by(analyzer, location, variable) %>%
                summarise(
                        mean = mean(mean, na.rm = TRUE),
                        sd = mean(sd, na.rm = TRUE),
                        cv = mean(cv, na.rm = TRUE),
                        .groups = "drop"
                ) %>%
                pivot_wider(
                        names_from = variable,
                        values_from = c(mean, sd, cv),
                        names_glue = "{.value}_{variable}"
                ) %>%
                mutate(across(where(is.numeric), ~ round(.x, 2)))
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

# Development of function emission and concentration plot
emiconplot <- function(data, y = NULL, var_type_filter = NULL, x = "day") {
        library(dplyr)
        library(ggplot2)
        library(scales)
        
        required_cols <- c("location", "analyzer", "var_type", "variable", "value", x)
        missing_cols <- setdiff(required_cols, names(data))
        if (length(missing_cols) > 0) {
                stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))
        }
        
        if (!is.null(var_type_filter)) {
                data <- data %>% filter(var_type %in% var_type_filter)
        }
        
        if (!is.null(y)) {
                data <- data %>% filter(variable %in% y)
        }
        
        summary_data <- data %>%
                group_by(across(all_of(c(x, "location", "analyzer", "var_type", "variable")))) %>%
                summarise(
                        mean_val = mean(value, na.rm = TRUE),
                        se_val = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
                        .groups = "drop"
                )
        
        # Default y-axis labels
        ylab_map <- list(
                concentration = expression(paste("Mean ± SD (mg ", m^{-3}, ")")),
                ratio         = "Mean ± SD (%)",
                ventilation   = expression(paste("Ventilation rate (m"^{-3} * " h"^{-1}, ")")),
                emission      = expression(paste("Emission (g h"^{-1}, ")")),
                error         = "Relative errors Mean ± SD (%)"
        )
        
        categories_present <- unique(summary_data$var_type)
        if (length(categories_present) == 1 && categories_present %in% names(ylab_map)) {
                ylab_to_use <- ylab_map[[categories_present]]
        } else {
                ylab_to_use <- "Mean ± SD"
        }
        
        # X-axis label logic
        xlab_to_use <- switch(
                x,
                "hour" = "Hour (aggregated by days)",
                "day" = "Day (aggregated by hours)",
                x
        )
        
        analyzer_colors <- c(
                "FTIR.1" = "#1b9e77", "FTIR.2" = "#d95f02", "FTIR.3" = "#7570b3",
                "FTIR.4" = "#e7298a", "CRDS.1" = "#66a61e", "CRDS.2" = "#e6ab02",
                "CRDS.3" = "#a6761d"
        )
        analyzer_shapes <- c(
                "FTIR.1" = 0, "FTIR.2" = 1, "FTIR.3" = 2,
                "FTIR.4" = 5, "CRDS.1" = 15, "CRDS.2" = 19, "CRDS.3" = 17
        )
        
        p <- ggplot(summary_data, aes_string(x = x, y = "mean_val", color = "analyzer", shape = "analyzer", group = "analyzer")) +
                geom_line() +
                geom_point(size = 2) +
                geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), width = 0.2) +
                scale_color_manual(values = analyzer_colors) +
                scale_shape_manual(values = analyzer_shapes) +
                facet_grid(variable ~ location, scales = "free_y", switch = "y") +
                scale_y_continuous(breaks = pretty_breaks(n = 8)) +
                labs(
                        x = xlab_to_use,
                        y = ylab_to_use
                ) +
                guides(
                        color = guide_legend(nrow = 1),
                        shape = guide_legend(nrow = 1)
                ) +
                theme_classic() +
                theme(
                        text = element_text(size = 12),
                        axis.text = element_text(size = 12),
                        axis.title = element_text(size = 12),
                        strip.text = element_text(size = 12),
                        panel.border = element_rect(colour = "black", fill = NA),
                        axis.text.x = element_text(hjust = 1),
                        legend.position = "bottom",
                        legend.direction = "horizontal",
                        legend.box = "horizontal",
                        legend.box.just = "left",
                        legend.title = element_blank(),
                        legend.text = element_text(size = 12),
                        legend.key.width = unit(0.1, "lines")
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

######## Import Gas Data #########
# Load animal and temperature data
animal_temp <- read.csv("20250408-15_LVAT_Animal_Temp_data.csv")

# Read and convert DATE.TIME for FTIR data
ATB_FTIR <- read.csv("20250408-15_long_ATB_FTIR.1.csv")
LUFA_FTIR <- read.csv("20250408-15_long_LUFA_FTIR.2.csv")
MBBM_FTIR <- read.csv("20250408-15_long_MBBM_FTIR.3.csv")
ANECO_FTIR <- read.csv("20250408-15_long_ANECO_FTIR.4.csv")

# Read and convert DATE.TIME for CRDS data
ATB_CRDS <- read.csv("20250408-15_long_ATB_CRDS.1.csv")
UB_CRDS <- read.csv("20250408-15_long_UB_CRDS.2.csv")
LUFA_CRDS <- read.csv("20250408-15_long_LUFA_CRDS.3.csv")

# Combine all data set
gas_data <- bind_rows(LUFA_FTIR, ANECO_FTIR, MBBM_FTIR, ATB_FTIR, ATB_CRDS, LUFA_CRDS, UB_CRDS)
input_combined <-full_join(gas_data, animal_temp, by = "DATE.TIME") %>%
        mutate(day = as.Date(DATE.TIME)) %>% 
        select(DATE.TIME, day, hour, everything())%>%
        mutate(across(where(is.numeric) & !any_of(c("hour")), ~ round(.x, 2)))

# Write csv
input_combined <- input_combined %>% select(DATE.TIME, hour, everything())
write_excel_csv(input_combined, "20250408-15_ringversuche_input_combined_data.csv")

######## Computation of ratios, ventilation rates and emissions #########
# Convert DATE.TIME format
input_combined <- input_combined %>% filter(DATE.TIME >= "2025-04-08 12:00:00" & DATE.TIME <= "2025-04-14 10:00:00")

# Calculate emissions using the function
emission_combined  <- indirect.CO2.balance(input_combined)

# Write csv
write_excel_csv(emission_combined, "20250408-15_ringversuche_emission_combined_data.csv")

######## Computation of relative errors #########
# List all variables
vars <- c(
        # Concentrations
        "CO2_in", "CO2_N", "CO2_S",
        "CH4_in", "CH4_N", "CH4_S",
        "NH3_in", "NH3_N", "NH3_S",
        
        # Ratios
        "NHCO_in", "NHCO_N", "NHCO_S",
        "CHCO_in", "CHCO_N", "CHCO_S",
        
        # Deltas
        "delta_NH3_N", "delta_CH4_N", "delta_CO2_N",
        "delta_NH3_S", "delta_CH4_S", "delta_CO2_S",
        
        # Ventilation
        "Q_Vent_rate_N", "Q_Vent_rate_S",
        
        # Emissions
        "e_NH3_N", "e_CH4_N", "e_NH3_S", "e_CH4_S"
)

# Calculate errors using the emirror function
emission_error  <- emirror(emission_combined, vars = vars)

# Write csv
write_excel_csv(emission_error, "20250408-15_ringversuche_emission_error_data.csv")

######## Reshape the data #########
emission_reshaped <-  reparam(emission_error) %>%
        mutate(
                DATE.TIME = as.POSIXct(DATE.TIME),
                day = as.factor(as.Date(DATE.TIME)),
                hour = as.factor(hour)) 

write_excel_csv(emission_reshaped, "20250408-15_ringversuche_emission_reshaped.csv")

########### Calculate mean, sd and cv #########
c_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("CO2", "CH4", "NH3","Err_CO2", "Err_CH4", "Err_NH3"),
        var_type_filter = c("concentration","error"),
        group_vars = c("analyzer", "location"))

d_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("delta_CO2", "delta_CH4", "delta_NH3", "Err_delta_CO2", "Err_delta_CH4", "Err_delta_NH3"),
        var_type_filter = c("delta","error"),
        group_vars = c("analyzer", "location"))

r_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("NH3/CO2", "CH4/CO2", "Err_NH3/CO2", "Err_CH4/CO2"),
        var_type_filter = c("ratio","error"),
        group_vars = c("analyzer", "location"))

q_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("Q_Vent_rate", "Err_Q_Vent_rate"),
        var_type_filter = c("ventilation","error"),
        group_vars = c("analyzer", "location"))

e_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("e_CH4", "e_NH3", "Err_e_CH4", "Err_e_NH3"),
        var_type_filter = c("emission","error"),
        group_vars = c("analyzer", "location"))

# Apply function to each stat table and create *_hour_sum objects
c_hour_sum <- hour_sum(c_stat_sum)
d_hour_sum <- hour_sum(d_stat_sum)
r_hour_sum <- hour_sum(r_stat_sum)
q_hour_sum <- hour_sum(q_stat_sum)
e_hour_sum <- hour_sum(e_stat_sum)

# write csv
readr::write_excel_csv(c_hour_sum, "c_hour_summary.csv")
readr::write_excel_csv(d_hour_sum, "d_hour_summary.csv")
readr::write_excel_csv(r_hour_sum, "r_hour_summary.csv")
readr::write_excel_csv(e_hour_sum, "e_hour_summary.csv")
readr::write_excel_csv(q_hour_sum, "q_hour_summary.csv")

######## ANOVA and HSD Summary ########
# Tukey HSD
# Remove rows where any response variable has NA/NaN/Inf
emission_clean <- emission_combined %>%
        filter(if_all(all_of(vars), ~ !is.na(.) & is.finite(.)))

# Then perform Tukey HSD
result_HSD_summary <- HSD_table(data = emission_clean,
                                response_vars = vars,
                                group_var = "analyzer")

# Save the result
write_excel_csv(result_HSD_summary, "20250408_20250414_HSD_table.csv")


######## Hourly Mean ± SD Trend Plots ########
c_r_erbr_plot <- emiconplot(data = emission_reshaped,
                          x = "hour",
                          y = c("Err_CO2", "Err_CH4", "Err_NH3", "Err_NH3/CO2", "Err_CH4/CO2"),
                          var_type_filter = "error")

q_e_erbr_plot <- emiconplot(data = emission_reshaped,
                          x = "hour",
                          y = c("Err_Q_Vent_rate","Err_e_CH4", "Err_e_NH3"),
                          var_type_filter = "error")

#save plots
ggsave(filename = "c_r_erbr_plot.pdf",
       plot = c_r_erbr_plot,
       width = 18, height = 12,
       dpi = 100)

ggsave(filename = "q_e_erbr_plot.pdf",
       plot = q_e_erbr_plot,
       width = 13, height = 7,
       dpi = 100)

########## Heat maps (CV) #############
c_heatmap <- emiheatmap(data = c_stat_sum,
                        response_vars = c("CO2", "CH4", "NH3"))

r_heatmap <- emiheatmap(data = r_stat_sum, 
                        response_vars = c("NH3/CO2", "CH4/CO2"))

q_heatmap <- emiheatmap(data = q_stat_sum, 
                        response_vars = c("Q_Vent_rate"))

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


