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

######## Development of functions #######
# Development of function reshape parameters
reparam <- function(data) {
        library(dplyr)
        library(tidyr)
        library(stringr)
        
        meta_cols <- c("DATE.TIME", "day", "hour", "lab", "analyzer")
        gases <- c("CO2", "CH4", "NH3")
        
        gas_pattern <- paste0(
                "(",
                paste(c(
                        paste0("^", gases, "_(in|N|S)$"),
                        paste0("^delta_", gases, "_(N|S)(_ppm|_mgm3)?$"),
                        paste0("^e_", gases, "_(N|S)$")
                ), collapse = "|"),
                ")"
        )
        
        gas_cols <- grep(gas_pattern, names(data), value = TRUE)
        
        # Add ventilation columns
        vent_cols <- c("Q_Vent_rate_N", "Q_Vent_rate_S")
        all_cols <- c(meta_cols, gas_cols, vent_cols[vent_cols %in% names(data)])
        
        if (length(gas_cols) == 0 && !any(vent_cols %in% names(data))) {
                stop("No gas-related or ventilation columns matched. Please check column names or patterns.")
        }
        
        df <- data %>% select(all_of(all_cols))
        
        long_df <- df %>%
                pivot_longer(
                        cols = -all_of(meta_cols),
                        names_to = "variable",
                        values_to = "value"
                ) %>%
                mutate(
                        gas = case_when(
                                variable %in% vent_cols ~ "Q",
                                TRUE ~ str_extract(variable, "CO2|CH4|NH3")
                        ),
                        type = case_when(
                                str_starts(variable, "e_") ~ "emission",
                                str_starts(variable, "delta_") ~ "delta",
                                variable %in% vent_cols ~ "Ventilation rate",
                                TRUE ~ "concentration"
                        ),
                        unit = case_when(
                                variable %in% vent_cols ~ "m^3 h^-1",
                                type == "emission" ~ "g h^-1",
                                str_detect(variable, "_mgm3") ~ "mg m^-3",
                                str_detect(variable, "_ppm") ~ "ppm",
                                TRUE ~ "ppm"
                        ),
                        location_code = str_extract(variable, "(?<=_)(in|N|S)(?=(_|$))"),
                        location = case_when(
                                location_code == "in" ~ "Ringline inside",
                                location_code == "N" ~ "North outside",
                                location_code == "S" ~ "South outside",
                                TRUE ~ NA_character_
                        )
                ) %>%
                select(all_of(meta_cols), location, type, unit, gas, value)
        
        wide_gas_df <- long_df %>%
                pivot_wider(
                        names_from = gas,
                        values_from = value
                ) %>%
                relocate(CO2, CH4, NH3, Q, .after = unit) %>%
                mutate(
                        DATE.TIME = as.POSIXct(DATE.TIME, format = "%Y-%m-%d %H:%M:%S"),
                        day = as.factor(day),
                        hour = as.factor(hour),
                        lab = as.factor(lab),
                        analyzer = as.factor(analyzer),
                        location = as.factor(location),
                        type = as.factor(type),
                        unit = as.factor(unit),
                        CO2 = as.numeric(CO2),
                        CH4 = as.numeric(CH4),
                        NH3 = as.numeric(NH3),
                        Q = as.numeric(Q)
                )
        
        return(wide_gas_df)
}

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

# Development of function stat_table
stat_table <- function(data, response_vars, group_vars) {
        require(dplyr)
        require(DescTools)
        
        data %>%
                group_by(across(all_of(group_vars))) %>%
                summarise(
                        n = n(),
                        across(
                                all_of(response_vars),
                                list(
                                        mean = ~round(mean(., na.rm = TRUE), 2),
                                        sd   = ~round(sd(., na.rm = TRUE), 2),
                                        cv   = ~round(DescTools::CoefVar(., na.rm = TRUE) * 100, 2)
                                ),
                                .names = "{.fn}_{.col}"
                        ),
                        .groups = "drop"
                )
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
emiconplot <- function(data, variable, type_filter = NULL, unit_filter = NULL) {
        library(dplyr)
        library(ggplot2)
        library(scales)
        
        # Ensure 'day' is Date class
        if ("day" %in% names(data)) {
                if (!inherits(data$day, "Date")) {
                        data <- data %>% mutate(day = as.Date(as.character(day)))
                }
        } else {
                stop("Data must have a 'day' column for faceting.")
        }
        
        # Filter by type if specified
        if (!is.null(type_filter) && "type" %in% names(data)) {
                data <- data %>% filter(type == type_filter)
        }
        
        # Filter by unit if specified
        if (!is.null(unit_filter) && "unit" %in% names(data)) {
                data <- data %>% filter(unit == unit_filter)
        }
        
        # Check variable exists
        if (!variable %in% names(data)) {
                stop(paste("Variable", variable, "not found in data"))
        }
        
        # Summarize mean and SE by day, hour, location, analyzer
        summary_data <- data %>%
                group_by(day, hour, location, analyzer) %>%
                summarise(
                        mean_val = mean(.data[[variable]], na.rm = TRUE),
                        se_val = sd(.data[[variable]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[variable]]))),
                        .groups = "drop"
                )
        
        analyzer_colors <- c(
                "FTIR.1" = "#1b9e77",
                "FTIR.2" = "#d95f02",
                "FTIR.3" = "#7570b3", 
                "FTIR.4" = "#e7298a",
                "CRDS.1" = "#66a61e",
                "CRDS.2" = "#e6ab02",
                "CRDS.3" = "#a6761d"
        )
        
        y_label <- variable
        if (!is.null(type_filter)) y_label <- paste0(y_label, " (", type_filter, ")")
        if (!is.null(unit_filter)) y_label <- paste0(y_label, " [", unit_filter, "]")
        
        p <- ggplot(summary_data, aes(x = hour, y = mean_val, color = analyzer, group = analyzer)) +
                geom_line() +
                geom_point(size = 1) +
                scale_color_manual(values = analyzer_colors) +
                labs(
                        x = "Hour",
                        y = y_label,
                        color = "Analyzer"
                ) +
                scale_y_continuous(breaks = pretty_breaks(n = 10)) +
                facet_grid(day ~ location, scales = "free_x") +
                theme_bw()
        
        print(p)
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
        select(DATE.TIME, day, hour, everything())


# Write csv
input_combined <- input_combined %>% select(DATE.TIME, hour, everything())
write.csv(input_combined, "20250408-15_ringversuche_input_combined_data.csv", row.names = FALSE)

######## Computation of ventilation rates and emissions #########
# Convert DATE.TIME format
input_combined <- input_combined %>% filter(DATE.TIME >= "2025-04-08 12:00:00" & DATE.TIME <= "2025-04-14 10:00:00")

# Calculate emissions using the function
emission_combined  <- indirect.CO2.balance(input_combined)

# Write csv
write.csv(emission_combined, "20250408-15_ringversuche_emission_combined_data.csv", row.names = FALSE)

######## Statistical Analysis ########
# list variables
cvars <- c(
        "CO2_in", "CH4_in", "NH3_in",
        "CO2_N", "CH4_N", "NH3_N",
        "CO2_S", "CH4_S", "NH3_S"
)

dvars <- c(
        "delta_CO2_N_ppm", "delta_CH4_N_ppm", "delta_NH3_N_ppm",
        "delta_CO2_S_ppm", "delta_CH4_S_ppm", "delta_NH3_S_ppm"
)

evars <- c(
        "e_CH4_N", "e_NH3_N", "e_CH4_S", "e_NH3_S"
)

qvars <- c(
        "Q_Vent_rate_N", "Q_Vent_rate_S"
)

vars <- c(
        "CO2_in", "CH4_in", "NH3_in",
        "CO2_N", "CH4_N", "NH3_N",
        "CO2_S", "CH4_S", "NH3_S",
        "delta_CO2_N_ppm", "delta_CH4_N_ppm", "delta_NH3_N_ppm",
        "delta_CO2_S_ppm", "delta_CH4_S_ppm", "delta_NH3_S_ppm",
        "e_CH4_N", "e_NH3_N", "e_CH4_S", "e_NH3_S",
        "Q_Vent_rate_N", "Q_Vent_rate_S"
        )

# Tukey HSD
result_HSD_summary <- HSD_table(data = emission_combined, response_vars = vars, group_var = "analyzer")
write_excel_csv(result_HSD_summary, "20250408_20250414_HSD_table.csv")


######## Trend Visualization ########
# Reshape the data
emission_reshaped <-  reparam(emission_combined) 
write_excel_csv(emission_reshaped, "20250408-15_ringversuche_emission_reshaped.csv")

# Emission plots (g h^-1)
eNH3 <- emiconplot(
        data = emission_reshaped, 
        variable = "NH3",
        type_filter = "emission",
        unit_filter = "g h^-1")

eCH4 <- emiconplot(
        data = emission_reshaped, 
        variable = "CH4",
        type_filter = "emission",
        unit_filter = "g h^-1")


# Concentration plots (ppm)
cNH3 <- emiconplot(
        data = emission_reshaped, 
        variable = "NH3",
        type_filter = "concentration",
        unit_filter = "ppm")

cCH4 <- emiconplot(
        data = emission_reshaped, 
        variable = "CH4",
        type_filter = "concentration",
        unit_filter = "ppm")

cCO2 <- emiconplot(
        data = emission_reshaped, 
        variable = "CO2",
        type_filter = "concentration",
        unit_filter = "ppm")

# Delta plots (ppm)
dNH3 <- emiconplot(
        data = emission_reshaped, 
        variable = "NH3",
        type_filter = "delta",
        unit_filter = "ppm")

dCH4 <- emiconplot(
        data = emission_reshaped, 
        variable = "CH4",
        type_filter = "delta",
        unit_filter = "ppm")

dCO2 <- emiconplot(
        data = emission_reshaped, 
        variable = "CO2",
        type_filter = "delta",
        unit_filter = "ppm")

qVent <-  emiconplot(
        data = emission_reshaped, 
        variable = "Q",
        type_filter = "Ventilation rate",
        unit_filter = "m^3 h^-1")


# Create a named list of all your plots and desired file names
dailyplots <- list(
        eNH3   = eNH3,
        eCH4   = eCH4,
        cNH3   = cNH3,
        cCH4   = cCH4,
        cCO2   = cCO2,
        dNH3   = dNH3,
        dCH4   = dCH4,
        dCO2   = dCO2,
        qVent  = qVent)

# Save each plot using ggsave
for (plot_name in names(dailyplots)) {
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = dailyplots[[plot_name]],
                width = 14, height = 12, dpi = 600)
}


######## Stats Visualization ########
# Concentration long
c_long <- emission_reshaped %>%
        filter(type == "concentration", unit == "ppm") %>%
        select(DATE.TIME, analyzer, location, CO2, CH4, NH3) %>%
        pivot_longer(cols = c(CO2, CH4, NH3), names_to = "gas", values_to = "concentration") %>%
        drop_na(concentration)

# Calculate mean, sd and cv
c_stat_sum <- stat_table(c_long, response_vars = "concentration", group_var = c("analyzer", "location", "gas"))


c_boxplot <- ggplot(c_long, aes(x = analyzer, y = concentration, fill = analyzer)) +
        geom_boxplot(outliers = FALSE) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "yellow") +
        facet_grid(gas ~ location, scales = "free_y") + 
        labs(title = "CO2, CH4 and NH3 Mean Concentration by Analyzer and Location",
             y = "(Concentration) [ppm]",
             fill = "Analyzer") +
        theme_bw() + theme(legend.position = "top")

# Delta long
d_long <- emission_reshaped %>%
        filter(type == "delta", unit == "ppm") %>%
        select(DATE.TIME, analyzer, location, CO2, CH4, NH3) %>%
        pivot_longer(cols = c(CO2, CH4, NH3), names_to = "gas", values_to = "concentration") %>%
        drop_na(concentration)

# Calculate mean, sd, and cv for delta concentrations
d_stat_sum <- stat_table(d_long, response_vars = "concentration", group_var = c("analyzer", "location", "gas"))

# Boxplot for delta concentrations
d_boxplot <- ggplot(d_long, aes(x = analyzer, y = concentration, fill = analyzer)) +
        geom_boxplot(outliers = FALSE) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "yellow") +
        facet_grid(gas ~ location, scales = "free_y") + 
        labs(
                title = "Delta Concentrations of CO2, CH4, and NH3 by Analyzer and Location",
                y = "(Delta Concentration) [ppm]",
                fill = "Analyzer"
        ) +
        theme_bw() +
        theme(legend.position = "top")

# Emission long
e_long <- emission_reshaped %>%
        filter(type == "emission", unit == "g h^-1") %>%
        select(DATE.TIME, analyzer, location, CH4, NH3) %>%
        pivot_longer(cols = c(CH4, NH3), names_to = "gas", values_to = "emission") %>%
        drop_na(emission)

# Calculate mean, sd, and cv for emissions
e_stat_sum <- stat_table(e_long, response_vars = "emission", group_var = c("analyzer", "location", "gas"))

# Boxplot for emissions
e_boxplot <- ggplot(e_long, aes(x = analyzer, y = emission, fill = analyzer)) +
        geom_boxplot(outliers = FALSE) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "yellow") +
        facet_grid(gas ~ location, scales = "free_y") +
        labs(
                title = "Emissions of CH4 and NH3 by Analyzer and Location",
                fill = "Analyzer"
        ) +
        theme_bw() +
        theme(legend.position = "top")


# Ventilation rate long
q_long <- emission_reshaped %>%
        filter(type == "Ventilation rate", unit == "m^3 h^-1") %>%
        select(DATE.TIME, analyzer, location, Q) %>%
        drop_na(Q)

# Calculate mean, sd, and cv for ventilation rate
q_stat_sum <- stat_table(q_long, response_vars = "Q", group_var = c("analyzer", "location"))

# Boxplot for ventilation rate
q_boxplot <- ggplot(q_long, aes(x = analyzer, y = Q, fill = analyzer)) +
        geom_boxplot(outlier.shape = NA) + 
        stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "yellow") +
        facet_wrap(~ location) +
        labs(
                title = "Ventilation Rate (Q) by Analyzer and Location",
                y = "Ventilation Rate [m³/h]",
                fill = "Analyzer"
        ) +
        theme_bw() +
        theme(legend.position = "top")

# Write stat summary as csv
readr::write_excel_csv(c_stat_sum, "c_stat_summary.csv")
readr::write_excel_csv(d_stat_sum, "d_stat_summary.csv")
readr::write_excel_csv(e_stat_sum, "e_stat_summary.csv")
readr::write_excel_csv(q_stat_sum, "q_stat_summary.csv")


# Named list of plots with specific dimensions
statplots <- list(
        concentration_boxplot = list(plot = c_boxplot, width = 12, height = 10),
        delta_boxplot         = list(plot = d_boxplot, width = 12, height = 10),
        emission_boxplot      = list(plot = e_boxplot, width = 8, height = 6),
        ventilation_boxplot   = list(plot = q_boxplot, width = 8,  height = 4)
)

# Save plots with their specific sizes
for (plot_name in names(statplots)) {
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot     = statplots[[plot_name]]$plot,
                width    = statplots[[plot_name]]$width,
                height   = statplots[[plot_name]]$height,
                dpi      = 600
        )
}

########## Correlograms ##############
# Change pivot to wide
eNH3_matrix <- emission_combined %>%
        select(DATE.TIME, analyzer, e_NH3_N) %>%
        pivot_wider(names_from = analyzer, values_from = e_NH3_N)

# Remove DATE.Time column
eNH3_matrix <- eNH3_matrix %>% select(-DATE.TIME)

# Compute correlation matrix (pairwise complete to allow NAs)
eNH3_matrix <- cor(eNH3_matrix, use = "pairwise.complete.obs")

# Plot the correlogram
ggcorrplot(eNH3_matrix, method = "circle", type = "lower", lab = TRUE,
           title = "Correlation of e_NH3_N Across Analyzers")


# Change pivot to wide
eCH4_matrix <- emission_combined %>%
        select(DATE.TIME, analyzer, e_CH4_N) %>%
        pivot_wider(names_from = analyzer, values_from = e_CH4_N)

# Remove DATE.Time column
eCH4_matrix <- eCH4_matrix %>% select(-DATE.TIME)

# Compute correlation matrix (pairwise complete to allow NAs)
eCH4_matrix <- cor(eCH4_matrix, use = "pairwise.complete.obs")

# Plot the correlogram
ggcorrplot(eCH4_matrix, method = "circle", type = "lower", lab = TRUE,
           title = "Correlation of e_CH4_N Across Analyzers")



result <- stat_table(emission_reshaped, c("CO2", "CH4", "NH3", "Q"), c("type", "location"))

