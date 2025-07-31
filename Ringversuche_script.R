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
library(networkD3)

######## Development of functions #######
# Development of indirect.CO2.balance function
indirect.CO2.balance <- function(df) {
        df %>%
                mutate(
                        P_CO2_term = 0.185,
                        a = 0.22,
                        h_min = 2.9,
                        hour = hour(DATE.TIME),
                        
                        A_corr = 1 - a * 3 * sin((2 * pi / 24) * (hour + 6 - h_min)),
                        Phi_tot = 5.6 * (m_weight ^ 0.75) + 22 * Y1_milk_prod + 1.6e-5 * (p_pregnancy_day ^ 3),
                        Phi_T_corr = Phi_tot * (1 + 4e-5 * (20 - Temperature)^3),
                        hpu_T_A_corr_all_animal = ((Phi_T_corr / 1000) * A_corr) * n_animals,
                        n_LU = (n_animals * m_weight) / 500,
                        P_CO2_T_A_all_animal = hpu_T_A_corr_all_animal * P_CO2_term,
                        
                        # Delta in ppm (original concentration differences)
                        delta_NH3_N = NH3_in - NH3_N,
                        delta_CH4_N = CH4_in - CH4_N,
                        delta_CO2_N = CO2_in - CO2_N,
                        
                        delta_NH3_S = NH3_in - NH3_S,
                        delta_CH4_S = CH4_in - CH4_S,
                        delta_CO2_S = CO2_in - CO2_S,
                        
                        # Ventilation rate North
                        Q_Vent_rate_N = P_CO2_T_A_all_animal / ((delta_CO2_N) * 1e-6),
                        # Ventilation rate South
                        Q_Vent_rate_S = P_CO2_T_A_all_animal / ((delta_CO2_S) * 1e-6),
                        
                        # Emissions North
                        e_NH3_N = (delta_NH3_N * Q_Vent_rate_N) / 1000,
                        e_CH4_N = (delta_CH4_N * Q_Vent_rate_N) / 1000,
                        
                        # Emissions South
                        e_NH3_S = (delta_NH3_S * Q_Vent_rate_S) / 1000,
                        e_CH4_S = (delta_CH4_S * Q_Vent_rate_S) / 1000,
                        )
}

# Development of function reshape parameters
reparam <- function(data) {
        library(dplyr)
        library(tidyr)
        library(stringr)
        
        meta_cols <- c("DATE.TIME", "day", "hour", "lab", "analyzer")
        gases <- c("CO2", "CH4", "NH3", "NHCO", "NHCH", "CHCO")
        vent_cols <- c("Q_Vent_rate_N", "Q_Vent_rate_S")
        
        gas_pattern <- paste0(
                "(",
                paste(c(
                        paste0("^", gases, "_(in|N|S)$"),
                        paste0("^delta_", c("CO2", "CH4", "NH3"), "_(N|S)$"),
                        paste0("^e_", c("CO2", "CH4", "NH3"), "_(N|S)$")
                ), collapse = "|"),
                ")"
        )
        
        gas_cols <- grep(gas_pattern, names(data), value = TRUE)
        all_vars <- c(meta_cols, gas_cols, vent_cols[vent_cols %in% names(data)])
        
        if(length(gas_cols) == 0 && !any(vent_cols %in% names(data))) {
                stop("No gas, emission, delta, ratio, or ventilation columns found in data.")
        }
        
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
                                str_starts(variable, "NHCO") ~ "ratio",
                                str_starts(variable, "NHCH") ~ "ratio",
                                str_starts(variable, "CHCO") ~ "ratio",
                                variable %in% vent_cols ~ "ventilation",
                                TRUE ~ "concentration"
                        ),
                        # Remove suffix (location code) from variable name:
                        variable = str_remove(variable, "_(in|N|S)$"),
                        # Add prefix c_ to base gases only for concentration var_type
                        variable = if_else(
                                var_type == "concentration" & variable %in% c("CO2", "CH4", "NH3"),
                                paste0("c_", variable),
                                variable
                        ),
                        # Rename ratios to human-readable form
                        variable = case_when(
                                variable == "NHCO" ~ "NH3/CO2",
                                variable == "NHCH" ~ "NH3/CH4",
                                variable == "CHCO" ~ "CH4/CO2",
                                TRUE ~ variable
                        )
                ) %>%
                select(all_of(meta_cols), location, var_type, variable, value)
        
        return(df_long)
}

# Development of function stat_table
stat_table <- function(data, response_vars, group_vars) {
        require(dplyr)
        require(DescTools)
        
        data %>%
                filter(variable %in% response_vars) %>%
                group_by(across(all_of(group_vars)), variable) %>%
                summarise(
                        n = n(),
                        mean = round(mean(value, na.rm = TRUE), 2),
                        sd = round(sd(value, na.rm = TRUE), 2),
                        cv = round(DescTools::CoefVar(value, na.rm = TRUE) * 100, 2),
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
emiconplot <- function(data, y = NULL, var_type_filter = NULL, x = "day") {
        library(dplyr)
        library(ggplot2)
        library(scales)
        
        # Check for required columns
        required_cols <- c("location", "analyzer", "var_type", "variable", "value", x)
        missing_cols <- setdiff(required_cols, names(data))
        if (length(missing_cols) > 0) {
                stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))
        }
        
        # Filter by var_type if requested
        if (!is.null(var_type_filter)) {
                data <- data %>% filter(var_type %in% var_type_filter)
        }
        
        # Filter by y (variable names) if requested
        if (!is.null(y)) {
                data <- data %>% filter(variable %in% y)
        }
        
        # Summarize mean and SE by grouping over x
        summary_data <- data %>%
                group_by(across(all_of(c(x, "location", "analyzer", "var_type", "variable")))) %>%
                summarise(
                        mean_val = mean(value, na.rm = TRUE),
                        se_val = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
                        .groups = "drop"
                )
        
        # y-axis labels by var_type
        ylab_map <- list(
                concentration = expression(paste("Mean ± SE (mg ", m^{-3}, ")")),
                ratio = "Mean ± SE (%)",
                ventilation = expression(paste("Ventilation rate (m"^3, "/s)")),
                emission = expression(paste("Emission (g h"^{-1}, ")"))
        )
        
        categories_present <- unique(summary_data$var_type)
        if (length(categories_present) == 1 && categories_present %in% names(ylab_map)) {
                ylab_to_use <- ylab_map[[categories_present]]
        } else {
                ylab_to_use <- "Mean ± SE"
        }
        
        # Colors and shapes for analyzers
        analyzer_colors <- c(
                "FTIR.1" = "#1b9e77", "FTIR.2" = "#d95f02", "FTIR.3" = "#7570b3",
                "FTIR.4" = "#e7298a", "CRDS.1" = "#66a61e", "CRDS.2" = "#e6ab02",
                "CRDS.3" = "#a6761d"
        )
        analyzer_shapes <- c(
                "FTIR.1" = 0, "FTIR.2" = 1, "FTIR.3" = 2,
                "FTIR.4" = 5, "CRDS.1" = 15, "CRDS.2" = 19, "CRDS.3" = 17
        )
        
        # Generate plot
        p <- ggplot(summary_data, aes_string(x = x, y = "mean_val", color = "analyzer", shape = "analyzer", group = "analyzer")) +
                geom_line() +
                geom_point(size = 2) +
                geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), width = 0.2) +
                scale_color_manual(values = analyzer_colors) +
                scale_shape_manual(values = analyzer_shapes) +
                facet_grid(variable ~ location, scales = "free_y", switch = "y") +
                scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
                labs(
                        x = x,
                        y = ylab_to_use,
                        color = "Analyzer",
                        shape = "Analyzer"
                ) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
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
        select(DATE.TIME, day, hour, everything())


# Write csv
input_combined <- input_combined %>% select(DATE.TIME, hour, everything())
write.csv(input_combined, "20250408-15_ringversuche_input_combined_data.csv", row.names = FALSE)

######## Computation of ratios, ventilation rates and emissions #########
# Convert DATE.TIME format
input_combined <- input_combined %>% filter(DATE.TIME >= "2025-04-08 12:00:00" & DATE.TIME <= "2025-04-14 10:00:00")

# Calculate emissions using the function
emission_combined  <- indirect.CO2.balance(input_combined)

emission_combined <- emission_combined %>%
        mutate(
                # NH3 to CO2
                NHCO_in = (NH3_in / CO2_in)*100,
                NHCO_N  = (NH3_N  / CO2_N)*100,
                NHCO_S  = (NH3_S  / CO2_S)*100,
                
                # NH3 to CH4
                NHCH_in = (NH3_in / CH4_in)*100,
                NHCH_N  = (NH3_N  / CH4_N)*100,
                NHCH_S  = (NH3_S  / CH4_S)*100,
                
                # CH4 to CO2
                CHCO_in = (CH4_in / CO2_in)*100,
                CHCO_N  = (CH4_N  / CO2_N)*100,
                CHCO_S  = (CH4_S  / CO2_S)*100
        )


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
        "delta_CO2_N", "delta_CH4_N", "delta_NH3_N",
        "delta_CO2_S", "delta_CH4_S", "delta_NH3_S"
)

evars <- c(
        "e_CH4_N", "e_NH3_N", "e_CH4_S", "e_NH3_S"
)

qvars <- c(
        "Q_Vent_rate_N", "Q_Vent_rate_S"
)

rvars <- c(
        "NHCO_in", "NHCO_N", "NHCO_S",
        "CHCO_in", "CHCO_N", "CHCO_S",
        "NHCH_in", "NHCH_N", "NHCH_S"
)


vars <- c(
        "CO2_in", "CH4_in", "NH3_in",
        "CO2_N", "CH4_N", "NH3_N",
        "CO2_S", "CH4_S", "NH3_S",
        "delta_CO2_N", "delta_CH4_N", "delta_NH3_N",
        "delta_CO2_S", "delta_CH4_S", "delta_NH3_S",
        "e_CH4_N", "e_NH3_N", "e_CH4_S", "e_NH3_S",
        "Q_Vent_rate_N", "Q_Vent_rate_S",
        "NHCO_in", "NHCO_N", "NHCO_S",
        "CHCO_in", "CHCO_N", "CHCO_S",
        "NHCH_in", "NHCH_N", "NHCH_S"
        )

# Tukey HSD
# Remove rows where any response variable has NA/NaN/Inf
emission_clean <- emission_combined %>%
        filter(if_all(all_of(vars), ~ !is.na(.) & is.finite(.)))

# Then perform Tukey HSD
result_HSD_summary <- HSD_table(
        data = emission_clean,
        response_vars = vars,
        group_var = "analyzer"
)

# Save the result
write_excel_csv(result_HSD_summary, "20250408_20250414_HSD_table.csv")


######## Trend Visualization ########
# Reshape the data
emission_reshaped <-  reparam(emission_combined) %>%
        mutate(
        DATE.TIME = as.POSIXct(DATE.TIME),
        day = as.factor(as.Date(DATE.TIME)),
        hour = as.factor(format(DATE.TIME, "%H:00"))) 

write_excel_csv(emission_reshaped, "20250408-15_ringversuche_emission_reshaped.csv")

# Concentration plots (ppm)
c_erbr_plot <- emiconplot(data = emission_reshaped,
                          x = "hour",
                          y = c("c_CO2", "c_CH4", "c_NH3"),
                          var_type_filter = "concentration")

r_erbr_plot <- emiconplot(data = emission_reshaped,
                          x = "hour",
                          y = c("NH3/CO2", "NH3/CH4", "CH4/CO2"),
                          var_type_filter = "ratio")

q_erbr_plot <- emiconplot(data = emission_reshaped,
                          x = "hour",
                          y = c("Q_Vent_rate"),
                          var_type_filter = "ventilation")

e_erbr_plot <- emiconplot(data = emission_reshaped,
                          x = "hour",
                          y = c("e_CH4", "e_NH3"),
                          var_type_filter = "emission")


# Named list of plots
dailyplots <- list(
        c_erbr_plot = c_erbr_plot,
        r_erbr_plot = r_erbr_plot,
        q_erbr_plot = q_erbr_plot,
        e_erbr_plot = e_erbr_plot)

# Corresponding size settings (width, height)
plot_sizes <- list(
        c_erbr_plot = c(15.5, 8.75),
        r_erbr_plot = c(15.5, 8.75),
        q_erbr_plot = c(13.5, 5.12),
        e_erbr_plot = c(12.69, 8.44))

# Save each plot using its specific size
for (plot_name in names(dailyplots)) {
        size <- plot_sizes[[plot_name]]
        ggsave(
                filename = paste0(plot_name, ".pdf"),
                plot = dailyplots[[plot_name]],
                width = size[1], height = size[2], dpi = 300)
        }


######## Stats Visualization ########
# Calculate mean, sd and cv
c_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("c_CO2", "c_CH4", "c_NH3"),
        group_vars = c("analyzer", "location"))

q_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("Q_Vent_rate"),
        group_vars = c("analyzer", "location"))

e_stat_sum <- stat_table(
        data = emission_reshaped,
        response_vars = c("e_CH4", "e_NH3"),
        group_vars = c("analyzer", "location"))

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



result <- stat_table(emission_reshaped, c("CO2", "CH4", "NH3", "Q"), c("var_type", "location"))



################# Sankey diagram ##############
# Create nodes dataframe for all unique nodes (conc, ventilation, emission)
conc_nodes <- c_stat_sum %>%
        mutate(node = paste0("Conc ", gas, " (", location, ")")) %>%
        distinct(node)

vent_nodes <- q_stat_sum %>%
        mutate(node = paste0("Ventilation rate (", location, ")")) %>%
        distinct(node)

emiss_nodes <- e_stat_sum %>%
        mutate(node = paste0("Emission ", gas, " (", location, ")")) %>%
        distinct(node)

# Combine all nodes and assign zero-based IDs
nodes <- bind_rows(conc_nodes, vent_nodes, emiss_nodes) %>%
        distinct() %>%
        mutate(id = row_number() - 1)

# Helper function to get node ID
node_id <- function(x) nodes$id[match(x, nodes$node)]

# Create links from concentration CV to ventilation rate CV by location and gas

# Because ventilation rate is not gas-specific, average concentration CVs per location (across gases)
conc_cv_by_location <- c_stat_sum %>%
        group_by(location) %>%
        summarise(mean_cv_conc = mean(cv_concentration, na.rm = TRUE))

vent_cv_by_location <- q_stat_sum %>%
        select(location, cv_Q) %>%
        distinct()

# Link conc -> ventilation
conc_to_vent <- conc_cv_by_location %>%
        left_join(vent_cv_by_location, by = "location") %>%
        rowwise() %>%
        mutate(
                source = node_id(paste0("Conc CO2 (", location, ")")),  # pick a representative gas node for source
                target = node_id(paste0("Ventilation rate (", location, ")")),
                value = mean_cv_conc
        ) %>%
        filter(!is.na(source), !is.na(target)) %>%
        select(source, target, value)

# Create links from ventilation rate CV to emission CV by location and gas
# For emission, keep separate gas nodes

vent_to_emiss <- e_stat_sum %>%
        rowwise() %>%
        mutate(
                source = node_id(paste0("Ventilation rate (", location, ")")),
                target = node_id(paste0("Emission ", gas, " (", location, ")")),
                value = cv_emission
        ) %>%
        filter(!is.na(source), !is.na(target)) %>%
        select(source, target, value)

# Combine all links
links <- bind_rows(conc_to_vent, vent_to_emiss)

# Remove links with NA or zero values
links <- links %>% filter(!is.na(value) & value > 0)

# Create Sankey Diagram
networkD3::sankeyNetwork(
        Links = links,
        Nodes = nodes,
        Source = "source",
        Target = "target",
        Value = "value",
        NodeID = "node",
        fontSize = 12,
        nodeWidth = 30,
        sinksRight = FALSE
)
