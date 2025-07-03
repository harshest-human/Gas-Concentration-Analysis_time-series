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
input_combined <-full_join(gas_data, animal_temp, by = "DATE.TIME")

# Write csv
input_combined <- input_combined %>% select(DATE.TIME, hour, everything())
write.csv(input_combined, "20250408-15_ringversuche_input_combined_data.csv", row.names = FALSE)

######## Computation of ventilation rates and emissions #########
# Convert DATE.TIME format
input_combined <- input_combined %>% filter(DATE.TIME >= "2025-04-08 12:00:00" & DATE.TIME <= "2025-04-14 23:00:00")

# Development of indirect.CO2.balance function
indirect.CO2.balance <- function(df) {
        
        library(dplyr)
        library(lubridate)
        
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
                        
                        # Annual emissions (tons/year), assuming constant emission rates
                        e_NH3_N_per_year = e_NH3_N * 24 * 365 / 1000,
                        e_CH4_N_per_year = e_CH4_N * 24 * 365 / 1000,
                        e_CO2_N_per_year = e_CO2_N * 24 * 365 / 1000,
                        
                        e_NH3_S_per_year = e_NH3_S * 24 * 365 / 1000,
                        e_CH4_S_per_year = e_CH4_S * 24 * 365 / 1000,
                        e_CO2_S_per_year = e_CO2_S * 24 * 365 / 1000
                )
}


# Calculate emissions using the function
emission_combined  <- indirect.CO2.balance(input_combined)

# Write csv
write.csv(emission_combined, "20250408-15_ringversuche_emission_combined_data.csv", row.names = FALSE)


######## Data visualization ############
# Analyzer color and fill mappings
analyzer_colors <- c(
        "FTIR.1" = "#1b9e77",
        "FTIR.2" = "#d95f02",
        "FTIR.3" = "#e7298a", 
        "FTIR.4" = "#7570b3",
        "CRDS.1" = "#66a61e",
        "CRDS.2" = "#e6ab02",
        "CRDS.3" = "#a6761d"
)

# Background side labels
bg <- list(
        N = "(computed by North side)",
        S = "(computed by South side)",
        
)

# Per year suffix
py <- list(per_year = "per year")

# Function to generate plotmath labels
get_plot_label <- function(varname) {
        match <- regexec("^(e|delta|Q_Vent_rate|err)_?(CH4|CO2|NH3)?_?(N|S)?(?:_(ppm|mgm3|per_year))?$", varname)
        parts <- regmatches(varname, match)[[1]]
        
        if (length(parts) == 0) return(varname)
        
        type <- parts[2]
        gas <- parts[3]
        side <- parts[4]
        suffix <- parts[5]
        
        gas_expr <- switch(gas,
                           CH4 = quote(CH[4]),
                           CO2 = quote(CO[2]),
                           NH3 = quote(NH[3]),
                           NULL
        )
        
        label <- switch(type,
                        e = bquote(.(gas_expr)~"Emission (g·h"^{-1}*")"),
                        delta = if (!is.null(suffix) && suffix %in% c("ppm", "mgm3")) {
                                if (suffix == "ppm") {
                                        bquote(Delta*.(gas_expr)~"(ppm)")
                                } else if (suffix == "mgm3") {
                                        bquote(Delta*.(gas_expr)~"(mg·m"^{-3}*")")
                                }
                        },
                        Q_Vent_rate = bquote("Ventilation rate (m"^3*"/h)"),
                        err = bquote("Relative error of "*.(gas_expr)*" (%)"),
                        varname
        )
        
        # Add background side
        if (!is.null(side) && nzchar(side)) {
                label <- bquote(.(label)~.(bg[[side]]))
        }
        
        # Add per_year if present
        if (!is.null(suffix) && suffix == "per_year") {
                label <- bquote(.(label)~.(py[["per_year"]]))
        }
        
        return(label)
}

# Main plotting function
emicon.plot <- function(df, x, y) {
        # Capture and convert to strings
        x_var <- enquo(x)
        y_var <- enquo(y)
        x_str <- as_label(x_var)
        y_str <- as_label(y_var)
        
        # Get math label for y axis
        y_label <- get_plot_label(y_str)
        
        # Base plot
        p <- ggline(
                df, x = x_str, y = y_str,
                add = "mean_se",
                color = "analyzer",
                palette = analyzer_colors
        ) +
                labs(
                        x = x_str,
                        y = y_label,
                        color = "Analyzer",
                        fill = "Analyzer"
                ) +
                scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
                theme_light() +
                theme(
                        legend.position = "right",
                        axis.title.y = element_text(size = 10),
                        plot.title = element_blank(),
                        plot.caption = element_text(size = 10, hjust = 0.5)
                )
        
        # If x is datetime, add datetime scale and rotate labels
        if (inherits(df[[x_str]], "POSIXct") || inherits(df[[x_str]], "POSIXt")) {
                p <- p +
                        scale_x_datetime(date_breaks = "6 hours", date_labels = "%d.%m %H:%M") +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
        
        # If x = "hour" (numeric from 0 to 23), set breaks for every hour
        if (x_str == "hour") {
                p <- p + 
                        scale_x_continuous(breaks = 0:23) +
                        theme(axis.text.x = element_text(angle = 0))
        }
        
        return(p)
}


# DATE.TIME conversion
result <- emission_combined %>%
        mutate(
                DATE.TIME = as.POSIXct(DATE.TIME, format = "%Y-%m-%d %H:%M:%S"),
                analyzer = as.factor(analyzer)
        )

# Concentration plots in ppm
emicon.plot(df = result, x = hour, y = delta_NH3_N_ppm)
emicon.plot(df = result, x = hour, y = delta_NH3_S_ppm)
emicon.plot(df = result, x = hour, y = delta_CH4_N_ppm)
emicon.plot(df = result, x = hour, y = delta_CH4_S_ppm)
emicon.plot(df = result, x = hour, y = delta_CO2_N_ppm)
emicon.plot(df = result, x = hour, y = delta_CO2_S_ppm)

# Concentration plots in mgm3
emicon.plot(df = result, x = hour, y = delta_NH3_N_mgm3)
emicon.plot(df = result, x = hour, y = delta_NH3_S_mgm3)
emicon.plot(df = result, x = hour, y = delta_CH4_N_mgm3)
emicon.plot(df = result, x = hour, y = delta_CH4_S_mgm3)
emicon.plot(df = result, x = hour, y = delta_CO2_N_mgm3)
emicon.plot(df = result, x = hour, y = delta_CO2_S_mgm3)

# Ventilation rate plots
emicon.plot(df = result, x = hour, y = Q_Vent_rate_N)
emicon.plot(df = result, x = hour, y = Q_Vent_rate_S)

# Emission plots
emicon.plot(df = result, x = hour, y = e_CH4_N)
emicon.plot(df = result, x = hour, y = e_CH4_S)
emicon.plot(df = result, x = hour, y = e_NH3_N)
emicon.plot(df = result, x = hour, y = e_NH3_S)



# save plots
y_vars <- c("Q_Vent_rate_N", "Q_Vent_rate_S",
            "delta_NH3_N_ppm", "delta_NH3_S_ppm",
            "delta_CH4_N_ppm", "delta_CH4_S_ppm",
            "delta_CO2_N_ppm", "delta_CO2_S_ppm",
            "delta_NH3_N_mgm3", "delta_NH3_S_mgm3",
            "delta_CH4_N_mgm3", "delta_CH4_S_mgm3",
            "delta_CO2_N_mgm3", "delta_CO2_S_mgm3",
            "e_CH4_N", "e_CH4_S")

for (y_var in y_vars) {
        y_sym <- sym(y_var)  # convert string to symbol
        
        p <- emicon.plot(df = result, x = hour, y = !!y_sym)
        
        file_name <- paste0(y_var, ".png")
        
        ggsave(filename = file_name,
               plot = p,
               width = 8,
               height = 5,
               dpi = 300)
        
        cat("Plot saved as", file_name, "\n")
}


######## Statistical Data Analysis #######
# Paired t-test between North and South delta NH3 (in ppm)
t.test(result$delta_NH3_N_ppm, result$delta_NH3_S_ppm, paired = TRUE)
t.test(result$delta_CH4_N_ppm, result$delta_CH4_S_ppm, paired = TRUE)
t.test(result$delta_CO2_N_ppm, result$delta_CO2_S_ppm, paired = TRUE)

# Paired t-tests for ventilation rate
t.test(result$Q_Vent_rate_N, result$Q_Vent_rate_S, paired = TRUE)

# Paired t-tests for emission
t.test(result$e_NH3_N, result$e_NH3_S, paired = TRUE)
t.test(result$e_CH4_N, result$e_CH4_S, paired = TRUE)
t.test(result$e_CO2_N, result$e_CO2_S, paired = TRUE)

######## Relative percentage errors #######
err <- result %>%
        select(
                DATE.TIME, hour, analyzer,
                delta_NH3_S_ppm, delta_CH4_S_ppm, delta_CO2_S_ppm,
                e_NH3_S, e_CH4_S, e_CO2_S
        ) %>%
        group_by(DATE.TIME, hour, analyzer) %>%
        summarise(
                delta_NH3_S_ppm = mean(delta_NH3_S_ppm, na.rm = TRUE),
                delta_CH4_S_ppm = mean(delta_CH4_S_ppm, na.rm = TRUE),
                delta_CO2_S_ppm = mean(delta_CO2_S_ppm, na.rm = TRUE),
                e_NH3_S = mean(e_NH3_S, na.rm = TRUE),
                e_CH4_S = mean(e_CH4_S, na.rm = TRUE),
                e_CO2_S = mean(e_CO2_S, na.rm = TRUE),
                .groups = "drop"
        ) %>%
        pivot_wider(
                names_from = analyzer,
                values_from = c(delta_NH3_S_ppm, delta_CH4_S_ppm, delta_CO2_S_ppm,
                                e_NH3_S, e_CH4_S, e_CO2_S),
                names_sep = "_"
        ) %>%
        rename_with(~str_remove_all(.x, "_S_ppm|_S"))

err <- err %>%
        mutate(
                # NH3 emissions
                e_NH3_FTIR.1_err = 100 * (e_NH3_FTIR.1 - e_NH3_CRDS.1) / e_NH3_CRDS.1,
                e_NH3_FTIR.2_err = 100 * (e_NH3_FTIR.2 - e_NH3_CRDS.1) / e_NH3_CRDS.1,
                e_NH3_FTIR.3_err = 100 * (e_NH3_FTIR.3 - e_NH3_CRDS.1) / e_NH3_CRDS.1,
                e_NH3_FTIR.4_err = 100 * (e_NH3_FTIR.4 - e_NH3_CRDS.1) / e_NH3_CRDS.1,
                e_NH3_CRDS.1_err = 100 * (e_NH3_CRDS.1 - e_NH3_CRDS.1) / e_NH3_CRDS.1,
                e_NH3_CRDS.2_err = 100 * (e_NH3_CRDS.2 - e_NH3_CRDS.1) / e_NH3_CRDS.1,
                e_NH3_CRDS.3_err = 100 * (e_NH3_CRDS.3 - e_NH3_CRDS.1) / e_NH3_CRDS.1,
                
                # CH4 emissions
                e_CH4_FTIR.1_err = 100 * (e_CH4_FTIR.1 - e_CH4_CRDS.1) / e_CH4_CRDS.1,
                e_CH4_FTIR.2_err = 100 * (e_CH4_FTIR.2 - e_CH4_CRDS.1) / e_CH4_CRDS.1,
                e_CH4_FTIR.3_err = 100 * (e_CH4_FTIR.3 - e_CH4_CRDS.1) / e_CH4_CRDS.1,
                e_CH4_FTIR.4_err = 100 * (e_CH4_FTIR.4 - e_CH4_CRDS.1) / e_CH4_CRDS.1,
                e_CH4_CRDS.1_err = 100 * (e_CH4_CRDS.1 - e_CH4_CRDS.1) / e_CH4_CRDS.1,
                e_CH4_CRDS.2_err = 100 * (e_CH4_CRDS.2 - e_CH4_CRDS.1) / e_CH4_CRDS.1,
                e_CH4_CRDS.3_err = 100 * (e_CH4_CRDS.3 - e_CH4_CRDS.1) / e_CH4_CRDS.1,
                
                # NH3 deltas
                delta_NH3_FTIR.1_err = 100 * (delta_NH3_FTIR.1 - delta_NH3_CRDS.1) / delta_NH3_CRDS.1,
                delta_NH3_FTIR.2_err = 100 * (delta_NH3_FTIR.2 - delta_NH3_CRDS.1) / delta_NH3_CRDS.1,
                delta_NH3_FTIR.3_err = 100 * (delta_NH3_FTIR.3 - delta_NH3_CRDS.1) / delta_NH3_CRDS.1,
                delta_NH3_FTIR.4_err = 100 * (delta_NH3_FTIR.4 - delta_NH3_CRDS.1) / delta_NH3_CRDS.1,
                delta_NH3_CRDS.1_err = 100 * (delta_NH3_CRDS.1 - delta_NH3_CRDS.1) / delta_NH3_CRDS.1,
                delta_NH3_CRDS.2_err = 100 * (delta_NH3_CRDS.2 - delta_NH3_CRDS.1) / delta_NH3_CRDS.1,
                delta_NH3_CRDS.3_err = 100 * (delta_NH3_CRDS.3 - delta_NH3_CRDS.1) / delta_NH3_CRDS.1,
                
                # CH4 deltas
                delta_CH4_FTIR.1_err = 100 * (delta_CH4_FTIR.1 - delta_CH4_CRDS.1) / delta_CH4_CRDS.1,
                delta_CH4_FTIR.2_err = 100 * (delta_CH4_FTIR.2 - delta_CH4_CRDS.1) / delta_CH4_CRDS.1,
                delta_CH4_FTIR.3_err = 100 * (delta_CH4_FTIR.3 - delta_CH4_CRDS.1) / delta_CH4_CRDS.1,
                delta_CH4_FTIR.4_err = 100 * (delta_CH4_FTIR.4 - delta_CH4_CRDS.1) / delta_CH4_CRDS.1,
                delta_CH4_CRDS.1_err = 100 * (delta_CH4_CRDS.1 - delta_CH4_CRDS.1) / delta_CH4_CRDS.1,
                delta_CH4_CRDS.2_err = 100 * (delta_CH4_CRDS.2 - delta_CH4_CRDS.1) / delta_CH4_CRDS.1,
                delta_CH4_CRDS.3_err = 100 * (delta_CH4_CRDS.3 - delta_CH4_CRDS.1) / delta_CH4_CRDS.1,
                
                # CO2 deltas
                delta_CO2_FTIR.1_err = 100 * (delta_CO2_FTIR.1 - delta_CO2_CRDS.1) / delta_CO2_CRDS.1,
                delta_CO2_FTIR.2_err = 100 * (delta_CO2_FTIR.2 - delta_CO2_CRDS.1) / delta_CO2_CRDS.1,
                delta_CO2_FTIR.3_err = 100 * (delta_CO2_FTIR.3 - delta_CO2_CRDS.1) / delta_CO2_CRDS.1,
                delta_CO2_FTIR.4_err = 100 * (delta_CO2_FTIR.4 - delta_CO2_CRDS.1) / delta_CO2_CRDS.1,
                delta_CO2_CRDS.1_err = 100 * (delta_CO2_CRDS.1 - delta_CO2_CRDS.1) / delta_CO2_CRDS.1,
                delta_CO2_CRDS.2_err = 100 * (delta_CO2_CRDS.2 - delta_CO2_CRDS.1) / delta_CO2_CRDS.1,
                delta_CO2_CRDS.3_err = 100 * (delta_CO2_CRDS.3 - delta_CO2_CRDS.1) / delta_CO2_CRDS.1
        )


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


# Usage
err <- rm_outliers_IQR(err, c(
        "e_NH3_FTIR.1_err", "e_NH3_FTIR.2_err", "e_NH3_FTIR.3_err", "e_NH3_FTIR.4_err",
        "e_NH3_CRDS.1_err", "e_NH3_CRDS.2_err", "e_NH3_CRDS.3_err",
        "e_CH4_FTIR.1_err", "e_CH4_FTIR.2_err", "e_CH4_FTIR.3_err", "e_CH4_FTIR.4_err",
        "e_CH4_CRDS.1_err", "e_CH4_CRDS.2_err", "e_CH4_CRDS.3_err",
        "delta_NH3_FTIR.1_err", "delta_NH3_FTIR.2_err", "delta_NH3_FTIR.3_err", "delta_NH3_FTIR.4_err",
        "delta_NH3_CRDS.1_err", "delta_NH3_CRDS.2_err", "delta_NH3_CRDS.3_err",
        "delta_CH4_FTIR.1_err", "delta_CH4_FTIR.2_err", "delta_CH4_FTIR.3_err", "delta_CH4_FTIR.4_err",
        "delta_CH4_CRDS.1_err", "delta_CH4_CRDS.2_err", "delta_CH4_CRDS.3_err",
        "delta_CO2_FTIR.1_err", "delta_CO2_FTIR.2_err", "delta_CO2_FTIR.3_err", "delta_CO2_FTIR.4_err",
        "delta_CO2_CRDS.1_err", "delta_CO2_CRDS.2_err", "delta_CO2_CRDS.3_err"
))


err_long <- err %>%
        select(DATE.TIME, hour, matches("^(e_NH3|e_CH4|delta_NH3|delta_CH4|delta_CO2).*_err$")) %>%
        pivot_longer(
                cols = -c(DATE.TIME, hour),
                names_to = c("gas", "analyzer"),
                names_pattern = "(e_NH3|e_CH4|delta_NH3|delta_CH4|delta_CO2)_(FTIR\\.\\d|CRDS\\.\\d)_err"
        ) %>%
        pivot_wider(
                id_cols = c(DATE.TIME, hour, analyzer),
                names_from = gas,
                values_from = value,
                names_glue = "{gas}_err"
        ) %>%
        mutate(analyzer = factor(analyzer, levels = unique(analyzer)))


# Plot boxplots comparing relative percentage error by gas and analyzer
ggplot(err_long, aes(x = analyzer, y = e_NH3_err, fill = analyzer)) +
        geom_jitter(width = 0.2, alpha = 0.3, color = "gray") +
        geom_boxplot() + theme_light() +
        scale_y_continuous(breaks = seq(-100, 100, by = 10))

ggplot(err_long, aes(x = analyzer, y = e_CH4_err, fill = analyzer)) +
        geom_jitter(width = 0.2, alpha = 0.3, color = "gray") +
        geom_boxplot() + theme_light() +
        scale_y_continuous(breaks = seq(-100, 100, by = 10))

ggplot(err_long, aes(x = analyzer, y = delta_NH3_err, fill = analyzer)) +
        geom_jitter(width = 0.2, alpha = 0.3, color = "gray") +
        geom_boxplot() + theme_light() +
        scale_y_continuous(breaks = seq(-100, 100, by = 10))

ggplot(err_long, aes(x = analyzer, y = delta_CH4_err, fill = analyzer)) +
        geom_jitter(width = 0.2, alpha = 0.3, color = "gray") +
        geom_boxplot() + theme_light() +
        scale_y_continuous(breaks = seq(-100, 100, by = 10))

ggplot(err_long, aes(x = analyzer, y = delta_CO2_err, fill = analyzer)) +
        geom_boxplot() +theme_light()


#save plots 
y_errs <- c("e_NH3_err", "e_CH4_err", "delta_NH3_err", "delta_CH4_err", "delta_CO2_err")

for (y_err in y_errs) {
        p <- ggplot(err_long, aes_string(x = "analyzer", y = y_err, fill = "analyzer")) +
                geom_boxplot() +
                theme_light() +
                labs(title = paste("Boxplot of", y_err, "by Analyzer"),
                     x = "Analyzer",
                     y = "Relative % Error") +
                theme(legend.position = "none",
                      axis.text.x = element_text(angle = 45, hjust = 1))
        
        file_name <- paste0(y_err, ".png")
        
        ggsave(filename = file_name,
               plot = p,
               width = 8,
               height = 5,
               dpi = 300)
        
        cat("Plot saved as", file_name, "\n")
}

