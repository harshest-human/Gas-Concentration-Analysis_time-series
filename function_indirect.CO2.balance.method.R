####### function to calculate emissions by indirect CO2 balancing method #########

indirect.CO2.balance.method <- function(df) {
        
        library(dplyr)
        library(lubridate)
        
        df %>%
                mutate(
                        # Add constants
                        P_CO2_term = 0.185,            # CO2 per hpu (m³/h/hpu)
                        a = 0.22,                      # amplitude coefficient
                        h_min = 2.9,                   # time of minimum activity
                        CO2_Molmass = 44.01,           # CO2 molar mass (g/mol)
                        NH3_Molmass = 17.031,          # NH3 molar mass (g/mol)
                        CH4_Molmass = 16.04,           # CH4 molar mass (g/mol)
                        R_gas_constant = 8.314472,     # gas constant (J/mol·K)
                        p_pressure_ref = 1013,         # reference pressure (mbar)
                        
                        # Extract hour from timestamp
                        hour = hour(DATE.TIME),
                        
                        # Activity correction
                        A_corr = 1 - a * 3 * sin((2 * pi / 24) * (hour + 6 - h_min)),
                        
                        # Heat production per animal
                        Phi_tot = 5.6 * (m_weight ^ 0.75) + 22 * Y1_milk_prod + 1.6e-5 * (p_pregnancy_day ^ 3),
                        
                        # Temperature-corrected heat production
                        Phi_T_corr = Phi_tot * (1 + 4e-5 * (20 - Temperature)^3),
                        
                        # hpu for all animals corrected for T & A
                        hpu_T_A_corr_all_animal = ((Phi_T_corr / 1000) * A_corr) * n_animals,
                        
                        # Livestock units
                        n_LU = (n_animals * m_weight) / 500,
                        
                        # CO2 production (m³/h)
                        P_CO2_T_A_all_animal = hpu_T_A_corr_all_animal * P_CO2_term,
                        
                        # Ventilation rate (m³/h)
                        Q_Vent_rate = P_CO2_T_A_all_animal / ((CO2_in - CO2_out) * 1e-6),
                        
                        # NH3 calculations
                        delta_NH3 = (0.1 * NH3_Molmass * p_pressure_ref * (NH3_in - NH3_out)) / ((Temperature + 273.15) * R_gas_constant),
                        emission_NH3 = (delta_NH3 * Q_Vent_rate) / 1000,
                        emission_NH3_per_year = emission_NH3 * 24 * 365 / 1000,
                        
                        # CH4 calculations
                        delta_CH4 = (0.1 * CH4_Molmass * p_pressure_ref * (CH4_in - CH4_out)) / ((Temperature + 273.15) * R_gas_constant),
                        emission_CH4 = (delta_CH4 * Q_Vent_rate) / 1000,
                        emission_CH4_per_year = emission_CH4 * 24 * 365 / 1000
                )
}
