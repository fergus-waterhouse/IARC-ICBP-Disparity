########################################################################
# ICBP INEQUALITIES IN NET SURVIVAL PROJECT CODE
# Author: Fergus Waterhouse
# Supervisors: Hanna Fink, Isabelle Soerjomataram
########################################################################

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(glue)
    library(broom)
})

# == CONFIGURATION =====================================================
# Set file paths
DIR_BASE   <- ""
FILE_SURV  <- file.path(DIR_BASE, "results/icbp_5y_ns_ses_agr.csv")
OUT_DETAIL <- file.path(DIR_BASE, "results/ses_detail_fig.svg")
OUT_SES    <- file.path(DIR_BASE, "results/ses_fig.svg")
OUT_AGE    <- file.path(DIR_BASE, "results/age_fig.svg")

registry_cols <- rev(c(
    "England"      = "#E64B35",
    "Scotland"     = "#00468B",
    "N Ireland"    = "chartreuse4",
    "Wales"        = "#AD002A", 
    "Ireland"      = "#00A087",
    "Victoria"     = "#FEB24C",
    "Queensland"   = "sienna2", 
    "New Zealand"  = "#925E9F" 
))


cancer_order <- rev(c(
    "Breast (C50)", "Rectum (C19-20)", "Colon (C18)", 
    "Ovary (...)", "Lung (C33-34)", "Oesophagus (C15)", 
    "Stomach (C16)", "Liver (C22)", "Pancreas (C25)"
))
registry_order <- names(registry_cols)

# == STEP 1: LOAD & PREP DATA ==========================================
# Applying factors globally here ensures ALL downstream plots use the correct order!
surv <- read.csv(FILE_SURV) %>% 
    mutate(
        cancer = factor(cancer, levels = cancer_order),
        registry = factor(registry, levels = registry_order)
    )

# == STEP 2: REGRESSION & DETAIL PLOT (SES) ============================
surv_ses <- surv %>% 
    filter(ses %in% 1:5) %>% 
    select(-agr)

ses_regression_results <- surv_ses %>% 
    mutate(
        # Standard Quintiles (1-5): map to 0, 0.25, 0.5, 0.75, 1
        # Queensland Tertiles (1-3): map to 1/12 (0.083), 1/2 (0.5), 11/12 (0.917)
        x = if_else(
            registry == "Queensland", 
            (1 / 12) + (ses - 1) * (5 / 12),  
            0.25 * (ses - 1)
        ),
        w = 1 / (netsurv.se^2)
    ) %>%
    group_by(registry, cancer) %>%
    reframe({
        fit <- lm(netsurv ~ x, weights = w)
        tidy_fit <- broom::tidy(fit, conf.int = TRUE)
        slope_row <- tidy_fit[tidy_fit$term == "x", ]
        
        # Predict at x = 0 (Equivalent to Quintile 1) and x = 1 (Equivalent to Quintile 5)
        preds <- predict(fit, newdata = data.frame(x = c(0, 1)), interval = "confidence")
        
        tibble(
            slope      = slope_row$estimate,
            slope.low  = slope_row$conf.low,
            slope.upp  = slope_row$conf.high,
            est_1      = preds[1, "fit"],
            est_1.low  = preds[1, "lwr"],
            est_1.upp  = preds[1, "upr"],
            est_5      = preds[2, "fit"],
            est_5.low  = preds[2, "lwr"],
            est_5.upp  = preds[2, "upr"]
        )
    })

surv_ses <- surv_ses %>% mutate(
    x = if_else(
        registry == "Queensland", 
        ((1 / 12) + (ses - 1) * (5 / 12)) * 4 + 1,  
        (0.25 * (ses - 1)) * 4 + 1
    )
)

p_detail <- ggplot(surv_ses, aes(x = x, y = netsurv, color = registry, fill = registry)) +
    geom_smooth(aes(weight = 1 / (netsurv.se^2)), method = "lm", alpha = 0.2, linewidth = 1) +
    geom_errorbar(aes(ymin = netsurv.low, ymax = netsurv.upp), width = 0.05, linewidth = 0.6) +
    geom_point(size = 2.5) +
    facet_grid(cancer ~ registry, scales = "free_y") +
    scale_x_continuous(breaks = seq(1, 5, 1), labels = c("1", "2", "3", "4", "5")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(values = registry_cols) +
    scale_fill_manual(values = registry_cols) +
    labs(x = "Deprivation Index", y = "5-Year Net Survival") +
    theme_minimal() +
    theme(legend.position = "none")

ggsave(OUT_DETAIL, p_detail, units = "px", width = 3000, height = 4000)


# == HELPER FUNCTION FOR LOLLIPOP PLOTS (AXES SWAPPED) =================
create_lollipop_plot <- function(data, col_empty, col_solid, col_seg_start, col_seg_end, 
                                 out_file, shape_empty = 21, size_empty = 2.5, shape_solid = 16, size_solid = 3.5) {
    
    n_reg <- nlevels(data$registry)
    dodge_width <- 0.75
    
    # Pre-calculate dodging on the X-axis
    plot_data <- data %>%
        mutate(
            cancer_num = as.numeric(cancer),
            registry_num = as.numeric(registry),
            dodge_offset = (registry_num - (n_reg + 1) / 2) * (dodge_width / n_reg),
            x_dodged = cancer_num + dodge_offset
        )
    
    cancer_labels <- levels(plot_data$cancer)
    
    p <- ggplot(plot_data, aes(color = registry)) +
        # 1. Background Zebra striping (Vertical now)
        geom_rect(
            data = data.frame(cancer_num = 1:length(cancer_labels)),
            aes(xmin = cancer_num - 0.4, xmax = cancer_num + 0.4, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.3, color = NA, inherit.aes = FALSE 
        ) +
        # geom_segment(
        #     aes(
        #         x = x_dodged,
        #         xend = x_dodged,
        #         y = .data[[col_seg_start]],
        #         yend = .data[[col_seg_end]] + 0.01
        #     ),
        #     arrow = arrow(length = unit(0.175, "cm"), type = "open"),
        #     linewidth = 0.8, alpha = 0.6
        # ) +
        geom_segment(
            aes(
                x = x_dodged,
                xend = x_dodged,
                y = .data[[col_seg_start]],
                yend = .data[[col_seg_end]]
            ),
            linewidth = 0.8, alpha = 0.6
        ) +
        # 3. Solid Shape
        geom_point(aes(x = x_dodged, y = .data[[col_solid]]), size = size_solid, shape = shape_solid) +
        # 4. Empty Shape
        geom_point(aes(x = x_dodged, y = .data[[col_empty]]), size = size_empty, shape = shape_empty, fill = "grey95", stroke = 1.2) +
        # Styling (Scales and Labels swapped)
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
        scale_x_continuous(breaks = 1:length(cancer_labels), labels = cancer_labels) +
        scale_color_manual(values = registry_cols) +
        labs(y = NULL, x = NULL, color = "Registry") +
        theme_minimal() +
        theme(
            panel.grid.major.y = element_blank(), 
            panel.grid.minor.y = element_blank(),
            axis.line.y = element_line(linewidth = 0.75, color = "grey30"),
            axis.ticks.y = element_line(linewidth = 0.75, color = "grey30"),
            # angle = 90 rotates text vertically (reading bottom-to-top)
            # hjust = 0.5 centers the text perfectly on each tick mark
            axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8), 
            text = element_text(size = 12),
            legend.position = "bottom",           
            legend.direction = "horizontal",      
            legend.box.margin = margin(t = 10),   
            legend.title = element_blank()        
        ) + coord_flip()
    
    ggsave(out_file, p, units = "px", width = 2660, height = 3800)
    return(p)
}


# == STEP 3: SES LOLLIPOP PLOT =========================================
p_lolly <- create_lollipop_plot(
    data          = ses_regression_results, 
    col_empty     = "est_1", 
    col_solid     = "est_5", 
    col_seg_start = "est_1", 
    col_seg_end   = "est_5", 
    out_file      = OUT_SES,
    shape_empty   = 21,  # Circle
    shape_solid   = 16   # Circle
)
print(p_lolly)


# == STEP 4: AGE LOLLIPOP PLOT =========================================
age_data <- surv %>% 
    filter(agr %in% c(1, 2)) %>% 
    select(registry, cancer, agr, netsurv) %>%
    pivot_wider(names_from = agr, values_from = netsurv, names_prefix = "age_")

age_lolly <- create_lollipop_plot(
    data          = age_data, 
    col_empty     = "age_1",    
    col_solid     = "age_2",    
    col_seg_start = "age_1", 
    col_seg_end   = "age_2", 
    out_file      = OUT_AGE,
    shape_empty   = 23,  # Fillable diamond
    size_empty    = 2,
    shape_solid   = 18,   # Solid diamond
    size_solid    = 4
)
print(age_lolly)



# == STEP 5: DISPARITY TABLES (SES & AGE) ==============================

# Helper function to format estimates and 95% CIs as a clean string
format_ci <- function(est, low, upp) {
    if_else(is.na(est), NA_character_, sprintf("%.3f [%.3f, %.3f]", est, low, upp))
}

# 1. SES DISPARITY TABLE
# The regression 'slope' perfectly represents the absolute gap (Most Deprived [x=1] - Least Deprived [x=0])
ses_disparity_table <- ses_regression_results %>%
    mutate(
        Least_Deprived_Surv = format_ci(est_1, est_1.low, est_1.upp),
        Most_Deprived_Surv  = format_ci(est_5, est_5.low, est_5.upp),
        Disparity_Gap       = format_ci(slope, slope.low, slope.upp) # Most deprived minus least deprived
    ) %>%
    select(registry, cancer, Least_Deprived_Surv, Most_Deprived_Surv, Disparity_Gap)

# Save SES Table
write.csv(
    ses_disparity_table, 
    file = file.path(DIR_BASE, "results/ses_disparity_table.csv"), 
    row.names = FALSE
)


# 2. AGE DISPARITY TABLE
# We pivot the main `surv` dataframe to get estimates and standard errors for both age groups on one row
age_disparity_raw <- surv %>%
    filter(agr %in% c(1, 2)) %>%
    select(registry, cancer, agr, netsurv, netsurv.se, netsurv.low, netsurv.upp) %>%
    pivot_wider(
        names_from = agr, 
        values_from = c(netsurv, netsurv.se, netsurv.low, netsurv.upp),
        names_sep = "_"
    )

age_disparity_table <- age_disparity_raw %>%
    mutate(
        Age_Group_1_Surv = format_ci(netsurv_1, netsurv.low_1, netsurv.upp_1),
        Age_Group_2_Surv = format_ci(netsurv_2, netsurv.low_2, netsurv.upp_2),
        
        # Calculate the gap (Age 2 - Age 1)
        gap_est = netsurv_2 - netsurv_1,
        
        # Calculate pooled standard error for the gap: sqrt(SE1^2 + SE2^2)
        gap_se  = sqrt(netsurv.se_1^2 + netsurv.se_2^2),
        
        # Calculate 95% CI for the gap
        gap_low = gap_est - (1.96 * gap_se),
        gap_upp = gap_est + (1.96 * gap_se),
        
        Disparity_Gap = format_ci(gap_est, gap_low, gap_upp)
    ) %>%
    select(registry, cancer, Age_Group_1_Surv, Age_Group_2_Surv, Disparity_Gap)

# Save Age Table
write.csv(
    age_disparity_table, 
    file = file.path(DIR_BASE, "results/age_disparity_table.csv"), 
    row.names = FALSE
)


# == STEP 6: AVOIDABLE DEATHS (SES & AGE) ==============================
# 1. Age-Associated Avoidable Deaths

age_ref <- df %>%
    filter(ses == 0 & agr == 1) %>%
    select(registry, cancer, R_ref = netsurv)

age_target <- df %>%
    filter(ses == 0 & agr == 2) %>%
    select(registry, cancer, 
           N_x = cases_in_period, 
           R_x = netsurv, 
           S_x_star = exp_surv)

age_ad <- age_target %>%
    inner_join(age_ref, by = c("registry", "cancer")) %>%
    mutate(
        age_avoidable_deaths = pmax(0, N_x * (R_ref - R_x) * S_x_star)
    ) %>%
    select(registry, cancer, age_avoidable_deaths)


# 2. SES-Associated Avoidable Deaths
# Get reference group (Least deprived: ses = 1)
ses_ref <- df %>%
    filter(agr == 0 & ses == 1) %>%
    select(registry, cancer, R_ref = netsurv)

# Get target groups (Other SES groups: ses = 2, 3, 4, 5)
ses_target <- df %>%
    filter(agr == 0 & ses %in% c(2, 3, 4, 5)) %>%
    select(registry, cancer, ses, 
           N_x = cases_in_period, 
           R_x = netsurv, 
           S_x_star = exp_surv)

# Calculate SES-associated Avoidable Deaths for each group, cap at 0, then sum
ses_ad <- ses_target %>%
    inner_join(ses_ref, by = c("registry", "cancer")) %>%
    mutate(
        ad_per_group = pmax(0, N_x * (R_ref - R_x) * S_x_star)
    ) %>%
    group_by(registry, cancer) %>%
    summarise(
        ses_avoidable_deaths = sum(ad_per_group),
        .groups = "drop"
    )


# 3. Excess Deaths
excess_df <- df %>%
    filter(agr == 0 & ses == 0) %>%
    select(registry, cancer, 
           N_total = cases_in_period, 
           os_km, 
           exp_surv) %>%
    mutate(
        # Cap at 0 in the rare event that observed survival > expected survival
        excess_deaths = pmax(0, N_total * (exp_surv - os_km))
    ) %>%
    select(registry, cancer, excess_deaths)


# 4. Combine Results & Calculate Proportions
total_summary <- age_ad %>%
    full_join(ses_ad, by = c("registry", "cancer")) %>%
    full_join(excess_df, by = c("registry", "cancer")) %>%
    mutate(
        # Calculate proportions
        # Using if_else ensures that if excess_deaths == 0, we don't get an Inf/NaN error
        prop_age_avoidable = if_else(excess_deaths > 0, age_avoidable_deaths / excess_deaths, 0),
        prop_ses_avoidable = if_else(excess_deaths > 0, ses_avoidable_deaths / excess_deaths, 0)
    )

print("Final Summary Dataset:")
print(total_summary)

# PLOT

# 1. Define the requested color palette
registry_cols <- rev(c(
    "England"      = "#E64B35",
    "Scotland"     = "#00468B",
    "N Ireland"    = "chartreuse4",
    "Wales"        = "#AD002A", 
    "Ireland"      = "#00A087",
    "Victoria"     = "#FEB24C",
    "Queensland"   = "sienna2", 
    "New Zealand"  = "#925E9F" 
))

cancer_order <- rev(c(
    "Breast (C50)", "Rectum (C19-20)", "Colon (C18)", 
    "Ovary (...)", "Lung (C33-34)", "Oesophagus (C15)", 
    "Stomach (C16)", "Liver (C22)", "Pancreas (C25)"
))

# 2. Define the new plotting function
create_avoidable_lollipop <- function(data, col_prop, col_abs, out_file, 
                                      shape_solid = 16, size_solid = 3.25) {
    
    # 1. NEW: Lock the registry factor levels to the exact order of your color palette
    data$registry <- factor(data$registry, levels = names(registry_cols))
    data$cancer <- factor(data$cancer, levels = cancer_order)
    if(!is.factor(data$cancer)) data$cancer <- as.factor(data$cancer)
    
    n_reg <- nlevels(data$registry)
    dodge_width <- 0.75
    
    # Pre-calculate dodging on the X-axis
    plot_data <- data %>%
        mutate(
            cancer_num = as.numeric(cancer),
            registry_num = as.numeric(registry),
            dodge_offset = (registry_num - (n_reg + 1) / 2) * (dodge_width / n_reg),
            x_dodged = cancer_num + dodge_offset
        )
    
    cancer_labels <- levels(plot_data$cancer)
    
    p <- ggplot(plot_data, aes(color = registry)) +
        # Background Zebra striping 
        geom_rect(
            data = data.frame(cancer_num = 1:length(cancer_labels)),
            aes(xmin = cancer_num - 0.4, xmax = cancer_num + 0.4, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.3, color = NA, inherit.aes = FALSE 
        ) +
        # Segment from 0 to the avoidable deaths proportion
        geom_segment(
            aes(
                x = x_dodged,
                xend = x_dodged,
                y = 0,
                yend = .data[[col_prop]]
            ),
            linewidth = 0.8, alpha = 0.6
        ) +
        # Solid Shape at the proportion value
        geom_point(
            aes(x = x_dodged, y = .data[[col_prop]]), 
            size = size_solid, shape = shape_solid
        ) +
        # Text representing absolute avoidable deaths
        geom_text(
            aes(
                x = x_dodged, 
                y = .data[[col_prop]], 
                label = scales::comma(round(.data[[col_abs]]))
            ),
            hjust = 0,  
            nudge_y = 0.005,
            size = 3.5,
            show.legend = FALSE
        ) +
        scale_y_continuous(
            labels = scales::percent_format(accuracy = 1),
            expand = expansion(mult = c(0, 0.15)) 
        ) +
        scale_x_continuous(breaks = 1:length(cancer_labels), labels = cancer_labels) +
        
        # 2. NEW: Force the legend to order starting from "England" to "New Zealand"
        scale_color_manual(
            values = registry_cols, 
            breaks = rev(names(registry_cols))
        ) +
        
        labs(y = "Proporiton of Excess Deaths Avoidable", x = NULL, color = "Registry") +
        theme_minimal() +
        theme(
            panel.grid.major.y = element_blank(), 
            panel.grid.minor.y = element_blank(),
            axis.line.y = element_line(linewidth = 0.75, color = "grey30"),
            axis.ticks.y = element_line(linewidth = 0.75, color = "grey30"),
            axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8), 
            text = element_text(size = 12),
            legend.position = "bottom",           
            legend.direction = "horizontal",      
            legend.box.margin = margin(t = 10),   
            legend.title = element_blank()        
        ) + 
        coord_flip()
    
    ggsave(out_file, p, units = "px", width = 2660, height = 3800)
    return(p)
}

# 3. Create the two separate plots

# 1. Age-Associated Avoidable Deaths Plot
plot_age <- create_avoidable_lollipop(
    data      = total_summary,
    col_prop  = "prop_age_avoidable",
    col_abs   = "age_avoidable_deaths",
    out_file  = "age_avoidable_deaths_plot.png"
)

# 2. SES-Associated Avoidable Deaths Plot
plot_ses <- create_avoidable_lollipop(
    data      = total_summary,
    col_prop  = "prop_ses_avoidable",
    col_abs   = "ses_avoidable_deaths",
    out_file  = "ses_avoidable_deaths_plot.png"
)