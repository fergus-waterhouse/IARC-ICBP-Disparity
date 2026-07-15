########################################################################
# ICBP INEQUALITIES IN NET SURVIVAL PROJECT CODE
# Author: Fergus Waterhouse
# Supervisors: Hanna Fink, Isabelle Soerjomataram
########################################################################

suppressPackageStartupMessages({
    library("tidyverse")
    library("glue")
    library("survival")
})

# NOTES:
# 1. IRELAND DLC NOT IN DATE FORMAT
# 2. All are 1 affluent - 5 deprived apart from England


# == STEP 1 ================== DEFINITIONS =============================

# ---- Analysis Approach ----
# TRUE = Period Approach (Left-truncated at start_date, right-censored at end_date)
# FALSE = Cohort Approach (Diagnosed between start_date & end_date, followed up to dlc)
use_period_approach <- TRUE

# ---- Cancer Sites ----
icbp_sites <- c(
    "Liver (C22)", 
    "Oesophagus (C15)",
    "Pancreas (C25)",
    "Lung (C33-34)",
    "Breast (C50)",
    "Ovary (...)",
    "Rectum (C19-20)",
    "Colon (C18)",
    "Stomach (C16)"
)

women_spec_sites <- c("Breast (C50)", "Ovary (...)")
men_spec_sites <- c()

# ---- Cohort ----
start_date <- dmy("01-01-2017")
end_date <- dmy("31-12-2021")

period <- 5 # e.g. 5 year net survival

# ---- Age Groups ----
age_splits <- c(-Inf, 70, Inf) # x-1 <= x < x+1
age_groups <- c(1, 2) # Age groups of interest e.g. 55-64 & 75-84

# ---- Age Standardisation ----
wcpp_age_splits <- c(-Inf, 50, 65, Inf) # Age groups 1: 15-49 2: 50-64 3: 65-99

wcpp_standards <- c(
    "Liver (C22)"      = 2, 
    "Oesophagus (C15)" = 2,
    "Pancreas (C25)"   = 2,
    "Lung (C33-34)"    = 2,
    "Rectum (C19-20)"  = 2,
    "Colon (C18)"      = 2,
    "Stomach (C16)"    = 2,
    "Breast (C50)"     = 3,
    "Ovary (...)"      = 3
)

# Standard Weights for age groups 1: 15-49 2: 50-64 3: 65-99 for each wcpp standard
wcpp_weights <- data.frame(
    wcpp_grp = c(
        "std1_age1", "std1_age2", "std1_age3",
        "std2_age1", "std2_age2", "std2_age3",
        "std3_age1", "std3_age2", "std3_age3"
    ),
    weight = c(
        0.7230, 0.1455, 0.1315,
        0.1213, 0.3242, 0.5545,
        0.3469, 0.3308, 0.3223
    )
)

# ---- Constants ----
days_in_year <- 365.241


# == STEP 2 ================== DIRECTORIES =============================

reg_path <- ""
out_path <- ""

reg_paths <- list(
    "England" = list(
        "surv" = glue("{reg_path}/EUUKM12_cc.dta"),
        "pop"  = glue("{reg_path}/popmort_EUUKM12_seg.dta"),
        "flip" = TRUE
    ),
    "N Ireland" = list(
        "surv" = glue("{reg_path}/EUUKMNIR_cc.dta"),
        "pop"  = glue("{reg_path}/popmort_EUUKMNIR_seg.dta"),
        "flip" = FALSE
    ),
    "Ireland" = list( 
        "surv" = glue("{reg_path}/EUIRL10_cc.dta"), 
        "pop"  = glue("{reg_path}/popmort_EUIRL10_seg.dta"),
        "flip" = FALSE
    ),
    "Victoria" = list(
        "surv" = glue("{reg_path}/OCAUSVIC_cc.dta"), 
        "pop"  = glue("{reg_path}/popmort_OCAUSVIC_seg.dta"),
        "flip" = FALSE
    ),
    "Scotland" = list(
        "surv" = glue("{reg_path}/EUUKMSCO_cc.dta"),
        "pop"  = glue("{reg_path}/popmort_EUUKMSCO_seg.dta"),
        "flip" = FALSE
    ),
    "Wales" = list(
        "surv" = glue("{reg_path}/EUUKMWAL_cc.dta"),
        "pop"  = glue("{reg_path}/popmort_EUUKMWAL_seg.dta"),
        "flip" = FALSE
    ),
    "Queensland" = list(
        "surv" = glue("{reg_path}/OCAUSQUE_cc.dta"),
        "pop"  = glue("{reg_path}/popmort_OCAUSQUE_seg.dta"),
        "flip" = FALSE
    ),
    "New Zealand" = list(
        "surv" = glue("{reg_path}/OCNEZ81_cc.dta"),
        "pop"  = glue("{reg_path}/popmort_OCNEZ81_seg.dta"),
        "flip" = FALSE
    )
)

approach_text <- ifelse(use_period_approach, "Period approach", "Standard Cohort approach")

cat(glue("

    ================ RUNNING ICBP INEQUALITIES ANALYSIS ================
    Author: Fergus Waterhouse
    Supervisors: Hanna Fink, Isabelle Soerjomataram
    ====================================================================
    
    For ICBP Phase 3 Sites: {glue_collapse(icbp_sites, sep = ', ')}
    
    Period: {start_date} to {end_date}
    Outcome: {period} year net & overall survival
    Estimator: Pohar-Perme ({approach_text}) & Kaplan-Meier
    
    ---------------- RUNNING SCRIPT ------------------------------------
    "))


# == STEP 3 ================== LOAD & FORMAT DATA ===============================
cat(paste0("\n\n [1] LOADING DATA FROM ", reg_path, " ..."))

# LOAD DATA
reg_data <- reg_paths %>%
    map(function(reg) {
        list(
            "surv" = haven::read_dta(reg$surv),
            "pop"  = haven::read_dta(reg$pop),
            "flip" = reg$flip
        )
    })

# FORMAT DATA
reg_data <- reg_data %>%
    map(function(reg) {
        
        # FOLLOW-UP DATA
        reg$surv <- reg$surv %>% 
            mutate(cancer = haven::as_factor(cancer)) %>%
            mutate(cancer = forcats::fct_recode(cancer,
                                                "Liver (C22)" = "Liver and intrahepatic bile ducts (C22)",
                                                "Lung (C33-34)" = "Trachea, bronchus and lung (C33-34)",
                                                "Ovary (...)" = "Ovary (C48.1-2, C56, C57.0)"
            )) %>%
            filter(cancer %in% icbp_sites) %>% # ONLY ICBP SELECT CANCER SITES
            mutate(
                agr = cut(age, breaks = age_splits, labels = F, right = F), 
                wcpp_agr = cut(age, breaks = wcpp_age_splits, labels = F, right = F), 
                wcpp_std = wcpp_standards[as.character(cancer)],
                wcpp_grp = paste0("std", wcpp_std, "_age", wcpp_agr)
            ) %>%
            mutate(ses = as.character(ses))
        
        if (reg$flip) {
            reg$surv <- reg$surv %>%
                mutate(ses = ifelse(ses == "Least deprived", "5", ses)) %>%
                mutate(ses = ifelse(ses == "Most deprived", "1", ses))
        } else {
            reg$surv <- reg$surv %>%
                mutate(ses = ifelse(ses == "Least deprived", "1", ses)) %>%
                mutate(ses = ifelse(ses == "Most deprived", "5", ses))
        }
        
        reg$surv <- reg$surv %>% 
            mutate(
                ses = as.integer(stringr::str_extract(ses, "\\d+")),
                ses = ifelse(ses == 99, NA, ses)
            )
        
        if(max(reg$surv$ses, na.rm = TRUE) == 10) {
            reg$surv <- reg$surv %>% mutate(ses = ceiling(ses / 2))
        }
        
        # RATETABLE
        reg$pop <- reg$pop %>%
            mutate(ses = as.character(seg)) %>%
            mutate(ses = ifelse(ses == "Least deprived", "1", ses)) %>%
            mutate(ses = ifelse(ses == "Most deprived", "5", ses)) %>%
            mutate(ses = as.integer(stringr::str_extract(ses, "\\d+")))
        
        if(max(reg$pop$ses, na.rm = TRUE) == 10) {
            reg$pop <- reg$pop %>%
                mutate(ses = ceiling(ses / 2)) %>%
                group_by(`_age`, `_year`, sex, ses) %>%
                summarise(rate = mean(rate, na.rm = TRUE), .groups = "drop") %>%
                as.data.frame()
        }
        
        reg$pop <- reg$pop %>%
            mutate(ses = as.numeric(ses))
        
        return(reg)
    })

# SEX-SPECIFIC
reg_data <- reg_data %>%
    map(function(reg) {
        reg$surv <- reg$surv %>% 
            filter(!(sex == 1 & cancer %in% women_spec_sites)) %>%
            filter(!(sex == 2 & cancer %in% men_spec_sites))
        return(reg)
    })

# REGISTRY-SPECIFIC DATE FIX
reg_data$Ireland$surv <- reg_data$Ireland$surv %>% mutate(
    dlc = as.Date(dlc, origin = "1960-01-01")
)

# FLIP SES VARIABLE EARLY (So it aligns in the raw_data snapshot below)
reg_data <- reg_data %>%
    map(function(reg) {
        if (reg$flip) {
            reg$surv <- mutate(reg$surv, ses = 6 - ses)
        }
        return(reg)
    })

# DETERMINE PRE-CLEANING ELIGIBILITY ("Identified" cohort)
reg_data <- reg_data %>%
    map(function(reg) {
        if (use_period_approach) {
            reg$surv <- reg$surv %>% 
                mutate(
                    tmp_time_in = pmax(0, as.numeric(start_date - doi)),
                    tmp_time_out = pmin(as.numeric(dlc - doi), as.numeric(end_date - doi)),
                    is_eligible = !is.na(doi) & !is.na(dlc) & (doi < dlc) & (tmp_time_out > tmp_time_in)
                )
        } else {
            reg$surv <- reg$surv %>% 
                mutate(
                    tmp_time_in = 0,
                    tmp_time_out = as.numeric(dlc - doi),
                    is_eligible = !is.na(doi) & !is.na(dlc) & (doi < dlc) & (doi >= start_date) & (doi <= end_date) & (tmp_time_out > tmp_time_in)
                )
        }
        return(reg)
    })

raw_reg_data <- reg_data


# == STEP 4 ================== DATA CLEANING ===========================
cat("\n [2] CLEANING DATA...")

clean_and_track <- function(surv_df, registry_name) {
    n0 <- nrow(surv_df)
    
    df1 <- surv_df %>% filter(exclude == 0)
    n1  <- nrow(df1)
    
    df2 <- df1 %>% filter(!is.na(doi), !is.na(dlc), !is.na(age))
    df2 <- df2 %>% filter(doi < dlc)
    n2  <- nrow(df2)
    
    df3 <- df2 %>% 
        arrange(id, desc(dlc)) %>% 
        distinct(id, .keep_all = TRUE) 
    n3  <- nrow(df3)
    
    df4 <- df3 %>% filter(ses %in% 1:5)
    n4 <- nrow(df4)
    
    tracking_table <- tibble(
        registry = registry_name,
        Step = c("Initial", "1. Exclude == 0", "2. Valid DOI/DLC/Age", 
                 "3. Deduplicate ID", "4. Valid SES"),
        Remaining_Rows = c(n0, n1, n2, n3, n4),
        Rows_Lost = c(0, n0-n1, n1-n2, n2-n3, n3-n4)
    )
    return(list(data = df4, tracking = tracking_table))
}

reg_data <- imap(reg_data, function(reg_data_list, reg_name) {
    clean <- clean_and_track(reg_data_list$surv, reg_name)
    reg_data_list$surv <- clean$data
    reg_data_list$track <- clean$tracking
    
    cat(paste0("\n        ------ ", reg_name, " ------------------\n"))
    tbl_text <- capture.output(print(as.data.frame(clean$tracking) %>% select(-registry), row.names = FALSE))
    cat(paste0("        ", tbl_text, collapse = "\n"))
    
    return(reg_data_list)
})


# == STEP 5 ================== DEFINE SURVIVAL VARIABLES ===============
cat("\n [3] DEFINING SURVIVAL VARIABLES...")

reg_data <- imap(reg_data, function(reg_list, reg_name) {
    if (use_period_approach) {
        reg_list$surv <- reg_list$surv %>% 
            mutate(
                event = status == 2 & dlc <= end_date,
                time_in = pmax(0, as.numeric(start_date - doi)),
                time_out = pmin(as.numeric(dlc - doi), as.numeric(end_date - doi))
            )
    } else {
        reg_list$surv <- reg_list$surv %>% 
            mutate(
                event = status == 2,
                time_in = 0,
                time_out = as.numeric(dlc - doi)
            )
    }
    return(reg_list)
})


# == STEP 6 ================== COHORT SELECTION ========================
cat("\n [4] COHORT SELECTION...")

select_and_track <- function(surv_df, registry_name, start, end, period_app) {
    n0 <- nrow(surv_df)
    df1 <- surv_df %>% filter(age > 15, age < 99)
    n1  <- nrow(df1)
    
    if (period_app) {
        df2 <- df1 %>% filter(time_out > time_in)
        n2  <- nrow(df2)
        tracking_table <- tibble(
            registry = registry_name,
            Step = c("1. 15 < Age < 99", "2. During Period", "Final"),
            Remaining_Rows = c(n1, n2, n2),
            Rows_Lost = c(n0-n1, n1-n2, "")
        )
    } else {
        df1b <- df1 %>% filter(doi >= start & doi <= end)
        n1b <- nrow(df1b)
        df2 <- df1b %>% filter(time_out > time_in)
        n2 <- nrow(df2)
        tracking_table <- tibble(
            registry = registry_name,
            Step = c("1. 15 < Age < 99", "2. Diagnosed in Window", "3. Valid Follow-up", "Final"),
            Remaining_Rows = c(n1, n1b, n2, n2),
            Rows_Lost = c(n0-n1, n1-n1b, n1b-n2, "")
        )
    }
    return(list(data = df2, tracking = tracking_table))
}

reg_data <- imap(reg_data, function(reg_data_list, reg_name) {
    selection <- select_and_track(reg_data_list$surv, reg_name, start_date, end_date, use_period_approach)
    reg_data_list$surv <- selection$data
    reg_data_list$track <- rbind(reg_data_list$track, selection$tracking)
    
    cat(paste0("\n        ------ ", reg_name, " ------------------\n"))
    tbl_text <- capture.output(print(as.data.frame(selection$tracking) %>% select(-registry), row.names = FALSE))
    cat(paste0("        ", tbl_text, collapse = "\n"))
    
    return(reg_data_list)
})


# == STEP 7 ================== RATETABLE SETUP =========================
cat("\n [5] PREPARING RATETABLES...")

prep_ratetable <- function(pop_df) {
    pop_df <- as.data.frame(pop_df)
    pop_df$age  <- as.numeric(pop_df$`_age`)
    pop_df$year <- as.numeric(pop_df$`_year`)
    pop_df$sex  <- as.numeric(pop_df$sex)
    pop_df$ses  <- as.numeric(pop_df$ses)
    
    pop_df <- pop_df[!is.na(pop_df$age) & !is.na(pop_df$year) & !is.na(pop_df$sex) & !is.na(pop_df$ses), ]
    pop_df$daily_rate <- as.numeric(pop_df$rate) / days_in_year
    
    ages  <- sort(unique(pop_df$age))
    years <- sort(unique(pop_df$year))
    sexes <- sort(unique(pop_df$sex))
    ses_groups <- sort(unique(pop_df$ses))
    
    rate_array <- tapply(
        pop_df$daily_rate, 
        list(pop_df$age, pop_df$year, pop_df$sex, pop_df$ses), 
        FUN = mean
    )
    
    attr(rate_array, "dimnames") <- list(
        age = as.character(ages), 
        year = as.character(years), 
        sex = as.character(sexes),
        ses = as.character(ses_groups)
    )
    
    attr(rate_array, "dimid") <- c("age", "year", "sex", "ses")
    attr(rate_array, "factor") <- c(0, 0, 1, 1) 
    attr(rate_array, "type") <- c(2, 3, 1, 1)
    
    attr(rate_array, "cutpoints") <- list(
        age = ages * days_in_year, 
        year = as.Date(paste0(years, "-01-01")), 
        sex = NULL,
        ses = NULL
    )
    
    attr(rate_array, "summary") <- function(R) paste("Ratetable with dimensions:", paste(dim(R), collapse = " x "))
    class(rate_array) <- c("ratetable", "matrix", "array")
    
    return(rate_array)
}

reg_data <- imap(reg_data, function(reg_list, reg_name) {
    reg_list$pop <- prep_ratetable(reg_list$pop)
    return(reg_list)
})


# == STEP 8 ================== SURVIVAL CALCULATIONS ===================
cat("\n [6] RUNNING SURVIVAL ANALYSIS...")

cat(paste0("\n\t\t", "Jurisdiction", "\t", "Site", "\t\tStatus"))

results_list <- list()
summary_list <- list()
prog <- 1

# Helper function to extract a value safely
get_val <- function(x) if (length(x) > 0) x else NA_real_

# Helper function to format survival with 95% CI strings
format_surv <- function(est, low, up) {
    if (is.na(est)) return(NA_character_)
    sprintf("%.3f [%.3f, %.3f]", est, low, up)
}

# --- HELPER FUNCTION FOR EXPECTED SURVIVAL ---
calc_expected_survival <- function(df, pop_rt, period_days) {
    
    period_days <- 1825 # TEMPORARY FIX -- THIS SHOULD BE SOLVED SO THAT IT DEPENDS ON SPECIFIED "PERIOD"
    
    w_df <- wcpp_weights %>% 
        filter(wcpp_grp %in% unique(df$wcpp_grp)) %>%
        mutate(wcpp_grp = as.character(wcpp_grp)) %>% 
        as.data.frame()
    w_df$weight <- w_df$weight / sum(w_df$weight)
    
    exp_fit <- survexp(
        time_at_risk ~ wcpp_grp, 
        data = df, 
        rmap = list(age = age_in_days, year = year_in, sex = sex, ses = ses), 
        method = "ederer",       
        ratetable = pop_rt,
        times = period_days 
    )
    
    exp_surv_mat <- matrix(exp_fit$surv, nrow = length(exp_fit$time))
    strata_names <- sub("wcpp_grp=", "", as.character(unique(df$wcpp_grp)))
    colnames(exp_surv_mat) <- strata_names
    
    exp_surv_df <- as.data.frame(exp_surv_mat)
    exp_surv_df$time <- exp_fit$time
    
    std_exp <- exp_surv_df %>%
        pivot_longer(cols = -time, names_to = "wcpp_grp", values_to = "exp_surv_stratum") %>%
        left_join(w_df, by = "wcpp_grp") %>%
        summarize(std_exp_surv = sum(exp_surv_stratum * weight, na.rm = TRUE), .groups = 'drop') %>%
        pull(std_exp_surv)
    
    return(std_exp)
}
# ---------------------------------------------


for (reg_name in names(reg_data)) {
    
    surv_df <- reg_data[[reg_name]]$surv
    pop_rt  <- reg_data[[reg_name]]$pop
    raw_df  <- raw_reg_data[[reg_name]]$surv 
    
    for (c_site in unique(surv_df$cancer)) {
        
        subset_df <- surv_df %>% filter(cancer == c_site) %>% as.data.frame()
        if (nrow(subset_df) == 0) next
        
        # PREP EXPECTED SURVIVAL VARIABLES ONCE FOR THE ENTIRE CANCER SITE COHORT
        subset_df$time_at_risk <- subset_df$time_out - subset_df$time_in
        subset_df$age_in_days  <- (subset_df$age * days_in_year) + subset_df$time_in
        subset_df$year_in      <- subset_df$doi + subset_df$time_in
        
        # FILTER RAW DATA FOR AVOIDABLE DEATHS INCIDENCE CALCULATION
        raw_sub <- raw_df %>% filter(cancer == c_site)
        incid_sub <- raw_sub %>% filter(doi >= start_date, doi <= end_date, age > 15, age < 99)
        
        cat(paste0("\n        ", "[", prog, "]\t", reg_name, "     \t", c_site, "\t"))
        flush.console()
        
        site_results <- list()
        site_summary <- list()
        
        # ---- 1. OVERALL NET, KM, & EXPECTED SURVIVAL ----
        tryCatch({
            w_ov <- wcpp_weights %>% filter(wcpp_grp %in% unique(subset_df$wcpp_grp)) %>%
                mutate(wcpp_grp = as.character(wcpp_grp)) %>% as.data.frame()
            w_ov$weight <- w_ov$weight / sum(w_ov$weight)
            
            # Diagnoses explicitly in the 5 year period
            cases_period_ov <- nrow(incid_sub)
            
            # Net Survival (Pohar-Perme)
            fit_ov <- relsurv::rs.surv(
                survival::Surv(time_in, time_out, event) ~ 1,
                rmap = list(age = age * 365.241, ses = ses), data = subset_df, ratetable = pop_rt,
                method = "pohar-perme", weight.table = w_ov, weight.names = c("wcpp_grp", "weight")
            )
            summ <- summary(fit_ov, times = period * days_in_year)
            
            # Overall Survival (Kaplan-Meier)
            fit_km <- survival::survfit(survival::Surv(time_in, time_out, event) ~ 1, data = subset_df)
            summ_km <- summary(fit_km, times = period * days_in_year, extend = TRUE)
            
            # Age-Standardized Expected Survival
            exp_est <- calc_expected_survival(subset_df, pop_rt, period * days_in_year)
            
            ns_est <- get_val(summ$surv)
            ns_low <- get_val(summ$lower)
            ns_upp <- get_val(summ$upper)
            
            os_est <- get_val(summ_km$surv)
            os_low <- get_val(summ_km$lower)
            os_upp <- get_val(summ_km$upper)
            
            # --- Standard Table Append ---
            site_results[[length(site_results) + 1]] <- data.frame(
                registry = reg_name, cancer = c_site, ses = 0, agr = 0,
                cases_in_period = cases_period_ov, # <-- INCIDENCE ADDED
                netsurv = ns_est, netsurv.se = get_val(summ$std.err),
                netsurv.low = ns_low, netsurv.upp = ns_upp,
                os_km = os_est, os_km.low = os_low, os_km.upp = os_upp,
                exp_surv = exp_est  
            )
            
            # --- Summary Table Append ---
            id_cnt  <- sum(raw_sub$is_eligible, na.rm = TRUE)
            inc_cnt <- nrow(subset_df)
            
            site_summary[[length(site_summary) + 1]] <- data.frame(
                registry = reg_name, cancer = c_site, ses = 0, agr = 0,
                Cases_in_Period = cases_period_ov, # <-- INCIDENCE ADDED
                Identified = id_cnt, 
                Excluded   = id_cnt - inc_cnt, 
                Included   = inc_cnt, 
                Mean_Age   = mean(subset_df$age, na.rm = TRUE),
                Net_Survival = format_surv(ns_est, ns_low, ns_upp),
                OS_KM        = format_surv(os_est, os_low, os_upp),
                Exp_Survival = if(is.na(exp_est)) NA_character_ else sprintf("%.3f", exp_est)
            )
            
        }, error = function(e) cat(paste0("[Overall FAILED] ", e$message)))
        
        
        # ---- 2. SES SPECIFIC NET, KM, & EXPECTED SURVIVAL ----
        for(s in sort(unique(subset_df$ses))) {
            tryCatch({
                ses_df <- subset_df %>% filter(ses == s)
                
                w_ses <- wcpp_weights %>% filter(wcpp_grp %in% unique(ses_df$wcpp_grp)) %>%
                    mutate(wcpp_grp = as.character(wcpp_grp)) %>% as.data.frame()
                w_ses$weight <- w_ses$weight / sum(w_ses$weight)
                
                # Diagnoses explicitly in the 5 year period
                cases_period_ses <- sum(incid_sub$ses == s, na.rm = TRUE)
                
                fit_ses <- relsurv::rs.surv(
                    survival::Surv(time_in, time_out, event) ~ 1,
                    rmap = list(age = age * 365.241, ses = ses), data = ses_df, ratetable = pop_rt,
                    method = "pohar-perme", weight.table = w_ses, weight.names = c("wcpp_grp", "weight")
                )
                summ <- summary(fit_ses, times = period * days_in_year)
                
                fit_km <- survival::survfit(survival::Surv(time_in, time_out, event) ~ 1, data = ses_df)
                summ_km <- summary(fit_km, times = period * days_in_year, extend = TRUE)
                
                # Age-Standardized Expected Survival for SES
                exp_est <- calc_expected_survival(ses_df, pop_rt, period * days_in_year)
                
                ns_est <- get_val(summ$surv)
                ns_low <- get_val(summ$lower)
                ns_upp <- get_val(summ$upper)
                
                os_est <- get_val(summ_km$surv)
                os_low <- get_val(summ_km$lower)
                os_upp <- get_val(summ_km$upper)
                
                # --- Standard Table Append ---
                site_results[[length(site_results) + 1]] <- data.frame(
                    registry = reg_name, cancer = c_site, ses = as.integer(s), agr = 0,
                    cases_in_period = cases_period_ses, # <-- INCIDENCE ADDED
                    netsurv = ns_est, netsurv.se = get_val(summ$std.err),
                    netsurv.low = ns_low, netsurv.upp = ns_upp,
                    os_km = os_est, os_km.low = os_low, os_km.upp = os_upp,
                    exp_surv = exp_est
                )
                
                # --- Summary Table Append ---
                id_cnt  <- sum(raw_sub$is_eligible & raw_sub$ses == s, na.rm = TRUE)
                inc_cnt <- nrow(ses_df)
                
                site_summary[[length(site_summary) + 1]] <- data.frame(
                    registry = reg_name, cancer = c_site, ses = as.integer(s), agr = 0,
                    Cases_in_Period = cases_period_ses, # <-- INCIDENCE ADDED
                    Identified = id_cnt, 
                    Excluded   = id_cnt - inc_cnt, 
                    Included   = inc_cnt, 
                    Mean_Age   = mean(ses_df$age, na.rm = TRUE),
                    Net_Survival = format_surv(ns_est, ns_low, ns_upp),
                    OS_KM        = format_surv(os_est, os_low, os_upp),
                    Exp_Survival = if(is.na(exp_est)) NA_character_ else sprintf("%.3f", exp_est)
                )
            }, error = function(e) cat(paste0("[SES ", s, " FAILED] ", e$message)))
        }
        
        
        # ---- 3. AGE SPECIFIC NET, KM, & EXPECTED SURVIVAL ----
        for(a in age_groups) {
            tryCatch({
                agr_df <- subset_df %>% filter(agr == a)
                if(nrow(agr_df) == 0) stop("No data")
                
                w_agr <- wcpp_weights %>% filter(wcpp_grp %in% unique(agr_df$wcpp_grp)) %>%
                    mutate(wcpp_grp = as.character(wcpp_grp)) %>% as.data.frame()
                w_agr$weight <- w_agr$weight / sum(w_agr$weight)
                
                # Diagnoses explicitly in the 5 year period
                cases_period_agr <- sum(incid_sub$agr == a, na.rm = TRUE)
                
                fit_agr <- relsurv::rs.surv(
                    survival::Surv(time_in, time_out, event) ~ 1,
                    rmap = list(age = age * 365.241, ses = ses), data = agr_df, ratetable = pop_rt,
                    method = "pohar-perme",  weight.table = w_agr, weight.names = c("wcpp_grp", "weight")
                )
                summ <- summary(fit_agr, times = period * days_in_year)
                
                fit_km <- survival::survfit(survival::Surv(time_in, time_out, event) ~ 1, data = agr_df)
                summ_km <- summary(fit_km, times = period * days_in_year, extend = TRUE)
                
                # Age-Standardized Expected Survival for Age Groups
                exp_est <- calc_expected_survival(agr_df, pop_rt, period * days_in_year)
                
                ns_est <- get_val(summ$surv)
                ns_low <- get_val(summ$lower)
                ns_upp <- get_val(summ$upper)
                
                os_est <- get_val(summ_km$surv)
                os_low <- get_val(summ_km$lower)
                os_upp <- get_val(summ_km$upper)
                
                # --- Standard Table Append ---
                site_results[[length(site_results) + 1]] <- data.frame(
                    registry = reg_name, cancer = c_site, ses = 0, agr = a,
                    cases_in_period = cases_period_agr, # <-- INCIDENCE ADDED
                    netsurv = ns_est, netsurv.se = get_val(summ$std.err),
                    netsurv.low = ns_low, netsurv.upp = ns_upp,
                    os_km = os_est, os_km.low = os_low, os_km.upp = os_upp,
                    exp_surv = exp_est 
                )
                
                # --- Summary Table Append ---
                id_cnt  <- sum(raw_sub$is_eligible & raw_sub$agr == a, na.rm = TRUE)
                inc_cnt <- nrow(agr_df)
                
                site_summary[[length(site_summary) + 1]] <- data.frame(
                    registry = reg_name, cancer = c_site, ses = 0, agr = a,
                    Cases_in_Period = cases_period_agr, # <-- INCIDENCE ADDED
                    Identified = id_cnt, 
                    Excluded   = id_cnt - inc_cnt, 
                    Included   = inc_cnt, 
                    Mean_Age   = mean(agr_df$age, na.rm = TRUE),
                    Net_Survival = format_surv(ns_est, ns_low, ns_upp),
                    OS_KM        = format_surv(os_est, os_low, os_upp),
                    Exp_Survival = if(is.na(exp_est)) NA_character_ else sprintf("%.3f", exp_est)
                )
            }, error = function(e) cat(paste0("[Age ", a, " FAILED] ")))
        }
        
        # ---- BIND RESULTS ----
        if(length(site_results) > 0) {
            results_list[[length(results_list) + 1]] <- bind_rows(site_results)
            summary_list[[length(summary_list) + 1]] <- bind_rows(site_summary)
            cat("COMPLETE")
        } else {
            cat("FAILED ALL")
        }
        
        prog <- prog + 1
        flush.console()
    }
}

final_ns_results <- bind_rows(results_list) %>%
    arrange(registry, cancer, ses, agr)

final_summary_results <- bind_rows(summary_list) %>%
    arrange(registry, cancer, ses, agr)


# == STEP 9 ================== SAVING RESULTS ===================

cat(paste0("\n [7] SAVING RESULTS TO ", out_path, " ..."))
write.csv(final_ns_results, file = paste0(out_path, "/icbp_", period, "y_ns_ses_agr.csv"), row.names = FALSE)
write.csv(final_summary_results, file = paste0(out_path, "/icbp_", period, "y_descriptive_summary.csv"), row.names = FALSE)

cat("\n\n================ ICBP INEQUALITIES ANALYSIS COMPLETE ================\n\n")
