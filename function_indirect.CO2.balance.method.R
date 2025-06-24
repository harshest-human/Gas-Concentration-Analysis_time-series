####### function to calculate emissions by indirect CO2 balancing method #########

indirect.CO2.balance.method  <- function(df) {
        
        library(dplyr)
        library(lubridate)
        library(stringr)
        
        # Detect gas concentration columns
        co2_in_col <- names(df)[str_detect(names(df), "CO2_in")]
        co2_n_col  <- names(df)[str_detect(names(df), "CO2_N")]
        co2_s_col  <- names(df)[str_detect(names(df), "CO2_S")]
        
        nh3_in_col <- names(df)[str_detect(names(df), "NH3_in")]
        nh3_n_col  <- names(df)[str_detect(names(df), "NH3_N")]
        nh3_s_col  <- names(df)[str_detect(names(df), "NH3_S")]
        
        ch4_in_col <- names(df)[str_detect(names(df), "CH4_in")]
        ch4_n_col  <- names(df)[str_detect(names(df), "CH4_N")]
        ch4_s_col  <- names(df)[str_detect(names(df), "CH4_S")]
        
        # Extract lab suffix (e.g., "LUFA") from column name
        lab_suffix <- str_extract(co2_in_col, "(?<=_in_).*")
        
        # Sanity check
        if (is.na(lab_suffix)) stop("Lab suffix not detected. Please ensure columns follow '*_in_*' format.")
        
        df %>%
                mutate(
                        # Constants
                        P_CO2_term = 0.185,
                        a = 0.22,
                        h_min = 2.9,
                        CO2_Molmass = 44.01,
                        NH3_Molmass = 17.031,
                        CH4_Molmass = 16.04,
                        R_gas_constant = 8.314472,
                        p_pressure_ref = 1013,
                        
                        # Time and correction
                        hour = hour(DATE.TIME),
                        A_corr = 1 - a * 3 * sin((2 * pi / 24) * (hour + 6 - h_min)),
                        
                        Phi_tot = 5.6 * (m_weight ^ 0.75) + 22 * Y1_milk_prod + 1.6e-5 * (p_pregnancy_day ^ 3),
                        Phi_T_corr = Phi_tot * (1 + 4e-5 * (20 - Temperature)^3),
                        hpu_T_A_corr_all_animal = ((Phi_T_corr / 1000) * A_corr) * n_animals,
                        n_LU = (n_animals * m_weight) / 500,
                        P_CO2_T_A_all_animal = hpu_T_A_corr_all_animal * P_CO2_term,
                        
                        # -------- NORTH --------
                        Q_Vent_rate_N = P_CO2_T_A_all_animal / ((.data[[co2_in_col]] - .data[[co2_n_col]]) * 1e-6),
                        
                        delta_NH3_N = (0.1 * NH3_Molmass * p_pressure_ref * (.data[[nh3_in_col]] - .data[[nh3_n_col]])) / 
                                ((Temperature + 273.15) * R_gas_constant),
                        !!paste0("emission_NH3_N_", lab_suffix) := (delta_NH3_N * Q_Vent_rate_N) / 1000,
                        !!paste0("emission_NH3_N_per_year_", lab_suffix) := get(paste0("emission_NH3_N_", lab_suffix)) * 24 * 365 / 1000,
                        
                        delta_CH4_N = (0.1 * CH4_Molmass * p_pressure_ref * (.data[[ch4_in_col]] - .data[[ch4_n_col]])) / 
                                ((Temperature + 273.15) * R_gas_constant),
                        !!paste0("emission_CH4_N_", lab_suffix) := (delta_CH4_N * Q_Vent_rate_N) / 1000,
                        !!paste0("emission_CH4_N_per_year_", lab_suffix) := get(paste0("emission_CH4_N_", lab_suffix)) * 24 * 365 / 1000,
                        
                        # -------- SOUTH --------
                        Q_Vent_rate_S = P_CO2_T_A_all_animal / ((.data[[co2_in_col]] - .data[[co2_s_col]]) * 1e-6),
                        
                        delta_NH3_S = (0.1 * NH3_Molmass * p_pressure_ref * (.data[[nh3_in_col]] - .data[[nh3_s_col]])) / 
                                ((Temperature + 273.15) * R_gas_constant),
                        !!paste0("emission_NH3_S_", lab_suffix) := (delta_NH3_S * Q_Vent_rate_S) / 1000,
                        !!paste0("emission_NH3_S_per_year_", lab_suffix) := get(paste0("emission_NH3_S_", lab_suffix)) * 24 * 365 / 1000,
                        
                        delta_CH4_S = (0.1 * CH4_Molmass * p_pressure_ref * (.data[[ch4_in_col]] - .data[[ch4_s_col]])) / 
                                ((Temperature + 273.15) * R_gas_constant),
                        !!paste0("emission_CH4_S_", lab_suffix) := (delta_CH4_S * Q_Vent_rate_S) / 1000,
                        !!paste0("emission_CH4_S_per_year_", lab_suffix) := get(paste0("emission_CH4_S_", lab_suffix)) * 24 * 365 / 1000
                )
}
