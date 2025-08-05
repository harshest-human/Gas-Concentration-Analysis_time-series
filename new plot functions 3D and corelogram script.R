###### Development of 3Dplot function #######
emiheatmap_3d <- function(data, response_vars, group_var = "analyzer", 
                          location_filter = NULL, variable_filter = NULL) {
        library(dplyr)
        library(plotly)
        
        # Filter for variables requested and remove FTIR.4 if present
        data_sub <- data %>%
                filter(variable %in% response_vars,
                       .data[[group_var]] != "FTIR.4")
        
        # Convert grouping vars to factors with known order
        data_sub[[group_var]] <- factor(data_sub[[group_var]], levels = sort(unique(data_sub[[group_var]])))
        data_sub$location <- factor(data_sub$location, levels = sort(unique(data_sub$location)))
        
        # Numeric positions for 3D axes
        data_sub$analyzer_num <- as.numeric(data_sub[[group_var]])
        data_sub$location_num <- as.numeric(data_sub$location)
        data_sub$hour_num <- as.numeric(data_sub$hour)
        
        # Plotly 3D scatter plot
        p <- plot_ly(data_sub,
                     x = ~hour_num,
                     y = ~analyzer_num,
                     z = ~location_num,
                     color = ~cv,
                     colors = viridis::viridis(100),
                     type = "scatter3d",
                     mode = "markers",
                     marker = list(size = 6)) %>%
                layout(scene = list(
                        xaxis = list(title = "Hour of Day", nticks = 24),
                        yaxis = list(title = group_var,
                                     tickvals = unique(data_sub$analyzer_num),
                                     ticktext = levels(data_sub[[group_var]])),
                        zaxis = list(title = "Location",
                                     tickvals = unique(data_sub$location_num),
                                     ticktext = levels(data_sub$location))
                ))
        
        print(p)
        return(p)
}

# Development of function correlogram plot
emicorrgram <- function(data, target_variable) {
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(scales)
        
        # Step 1: Filter & assign side
        filtered <- data %>%
                filter(variable == target_variable) %>%
                mutate(
                        side = case_when(
                                grepl("North", location, ignore.case = TRUE) ~ "North",
                                grepl("South", location, ignore.case = TRUE) ~ "South",
                                TRUE ~ NA_character_
                        )
                ) %>%
                filter(!is.na(side)) %>%
                select(DATE.TIME, analyzer, side, var_type, value)
        
        # Get var_type label for legend/title
        vtype <- unique(filtered$var_type)
        vtype <- ifelse(length(vtype) == 1, vtype, "Value")
        title_label <- paste0("Correlation of ", target_variable, " (", vtype, ")")
        
        # Step 2: Compute correlation for each side and reshape long
        cor_long <- filtered %>%
                group_by(side) %>%
                group_modify(~ {
                        pivoted <- .x %>%
                                pivot_wider(names_from = analyzer, values_from = value) %>%
                                select(-DATE.TIME) %>%
                                mutate(across(everything(), as.numeric))
                        
                        cor_mat <- cor(pivoted, use = "pairwise.complete.obs")
                        cor_mat[upper.tri(cor_mat, diag = TRUE)] <- NA
                        
                        as.data.frame(as.table(cor_mat)) %>%
                                filter(!is.na(Freq)) %>%
                                rename(Var1 = Var1, Var2 = Var2, correlation = Freq)
                }) %>%
                ungroup()
        
        # Define breaks and matching colors (11 each)
        breaks <- seq(-1, 1, length.out = 21)
        
        colors <- c(
                "darkred",    # -1.0
                "darkred",    # -0.9
                "darkred",    # -0.8
                "orangered4", # -0.7
                "orangered4", # -0.6
                "orangered3", # -0.5
                "orangered3", # -0.4
                "orangered3", # -0.3
                "orangered2", # -0.2
                "orange",     # -0.1
                "white",      #  0.0
                "yellow",     #  0.05
                "yellow",     #  0.1
                "yellow",     #  0.2
                "gold",       #  0.3
                "gold",       #  0.4
                "gold",       #  0.5
                "gold3",      #  0.6
                "green3",     #  0.7
                "green4",     #  0.8
                "darkgreen"   #  1.0
        )
        values <- scales::rescale(breaks, to = c(0, 1))
        
        p <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = correlation)) +
                geom_tile(color = "white") +
                geom_text(aes(label = round(correlation, 2)), size = 3) +
                scale_fill_gradientn(
                        colors = colors,
                        values = values,
                        limits = c(-1, 1),
                        name = "Correlation"
                ) +
                labs(title = title_label, x = NULL, y = NULL) +
                theme_classic() +
                theme(
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        panel.border = element_rect(colour = "black", fill = NA)
                ) +
                facet_wrap(~ side)
        
        print(p)
}

