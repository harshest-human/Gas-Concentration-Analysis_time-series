####### function to calculate emissions by indirect CO2 balancing method #########
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
