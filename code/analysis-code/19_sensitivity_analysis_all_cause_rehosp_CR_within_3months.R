


#Option to pull data or refresh from data cleaning
clean_cohort <- read.csv("clean_cohort.csv")



library(survival)
library(dplyr)
library(lubridate)
library(mice)
library(broom)
library(readr)



# Load imputation object and recreate completed datasets
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")
completed_datasets <- lapply(1:5, function(i) complete(imp, i))


##--------------------------------------------------
## 0. Add outcome + date vars from clean_cohort
##    (run once before modelling)
##    Assumes clean_cohort has:
##    - person_id
##    - all_cause_rehosp_censor_time
##    - out_readm_allcause_flag
##    - tavi_cips_disdate (Date)
##--------------------------------------------------

# Create a lookup table with person_id and all-cause outcome variables
allcause_lookup <- clean_cohort %>%
  select(
    person_id,
    out_readm_all_cause_date,   # date of all-cause rehosp
    out_readm_all_cause_flag,
    date_of_death   
  )

# Join these into each completed dataset  
completed_datasets <- lapply(
  completed_datasets,
  function(df) {
    df %>% left_join(allcause_lookup, by = "person_id") %>%
      mutate(
        year_tavi_discharge = factor(year(tavi_cips_disdate)),
        tavi_cips_disdate = as.Date(tavi_cips_disdate),
        out_readm_all_cause_date = as.Date(out_readm_all_cause_date),
        date_of_death = as.Date(date_of_death),
        rehab_indicator = ifelse(is.na(exp_cardiac_rehab_3m) | exp_cardiac_rehab_3m == 0, 0, 1),
        all_cause_rehosp_indicator = ifelse(
          is.na(out_readm_all_cause_flag) | out_readm_all_cause_flag == 0, 0, 1
        )
      )
  }
)



#--------------------------------------------------
# Helper function: run all-cause rehosp analysis
#   on one completed dataset
#--------------------------------------------------
run_allcause_90d_sensitivity <- function(df) {
  
  study_end <- as.Date("2024-03-31")
  
  
  df <- df %>%
    mutate(
      tavi_cips_disdate_plus90_days = tavi_cips_disdate + days(90),
      study_start_new = tavi_cips_disdate_plus90_days,
      
      # Time from 90 days post-discharge to event/censoring
      allcause_rehosp_censor_time = pmin(
        as.numeric(out_readm_all_cause_date - study_start_new),  # <-- change name if needed
        as.numeric(date_of_death         - study_start_new),
        as.numeric(study_end             - study_start_new),
        na.rm = TRUE
      ),
      
      # This is 0 for everyone (time origin), used for filtering early events
      days_since_tavi_dis_plus90 = as.numeric(tavi_cips_disdate_plus90_days - study_start_new),
      
      # Rehab exposure within 3 months
      rehab_indicator = ifelse(is.na(exp_cardiac_rehab_3m) | exp_cardiac_rehab_3m == 0, 0, 1),
      
      # All-cause rehosp indicator (1 = rehosp, 0 = censored)
      #   Change out_readm_allcause_flag to your actual flag variable if different
      allcause_rehosp_indicator = ifelse(is.na(out_readm_all_cause_flag) | 
                                           out_readm_all_cause_flag == 0, 0, 1)
    ) %>%
    # Drop people with event/censoring before or at 90 days post discharge
    filter(allcause_rehosp_censor_time > days_since_tavi_dis_plus90)
  
  # Year of TAVI discharge (for calendar time)
  df <- df %>%
    mutate(year_tavi_discharge = factor(year(tavi_cips_disdate)))
  
  #-------------------------------
  # 1) Unadjusted model
  #-------------------------------
  mod_unadj <- coxph(
    Surv(allcause_rehosp_censor_time, allcause_rehosp_indicator) ~ rehab_indicator,
    data = df
  )
  
 
  #-------------------------------
  # 2) Fully adjusted model with stratification (same covariates as HF)
  #-------------------------------
  mod_full <- coxph(
    Surv(allcause_rehosp_censor_time, allcause_rehosp_indicator) ~
      rehab_indicator + sex_code  + IMD_2019_QUINTILES + SMOKING_STATUS + region_name +
      BMI_category + KATZ_INDEX_CAT + binary_known_poor_mobility + TAVI_procedural_complications +
      previous_MI_time_grouped + binary_known_prev_cardiac_surgery + binary_Known_DIABETES + ethnicity_5_group + LV_FUNCTION + MITRAL_REGURGITATION +
      strata(age_groups, NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE, year_tavi_discharge, COMORBIDITY_pulmonary_liver_neuro_excardiacvascu)
    ,
    data = df
  )
  
  list(
    data      = df,
    unadj     = mod_unadj,
    full      = mod_full
  )
}

#--------------------------------------------------
# Apply to all 5 imputed datasets
#
#   imp <- readRDS(".../imp_mice_object.rds")
#   completed_datasets <- lapply(1:5, function(i) complete(imp, i))
#--------------------------------------------------

allcause_analyses <- lapply(completed_datasets, run_allcause_90d_sensitivity)






###################  Pool the Cox model results across the 5 imputations

# Helper: turn list of models into mira object
get_mira <- function(analyses_list, element) {
  models <- lapply(analyses_list, `[[`, element)
  as.mira(models)
}

mira_unadj   <- get_mira(allcause_analyses, "unadj")
mira_full    <- get_mira(allcause_analyses, "full")

pool_unadj   <- pool(mira_unadj)
pool_full    <- pool(mira_full)

#--------------------------------------------------
# Tidy & export pooled results (HR, 95% CI, p) 
#   + number of events (rounded to nearest 5)
#--------------------------------------------------

tidy_pooled <- function(pool_obj, analyses_list, model_name, out_path) {
  
  # Summary from mice::pool (Rubin's rules)
  s <- summary(pool_obj, conf.int = TRUE, exponentiate = TRUE)
  # Columns are: term, estimate, std.error, statistic, df, p.value, 2.5 %, 97.5 %
  
  res <- s %>%
    rename(
      HR       = estimate,
      lower_CI = `2.5 %`,
      upper_CI = `97.5 %`,
      p_value  = p.value
    ) %>%
    mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3)) %>%
    mutate(extra_info = NA_character_) %>%
    select(term, HR, std.error, statistic, p_value, lower_CI, upper_CI, extra_info)
  
  # Number of events: should be identical across imputations (outcome is observed)
  events_vec <- sapply(analyses_list, function(x) summary(x[[model_name]])$nevent)
  rounded_events <- 5 * round(mean(events_vec) / 5)
  
  events_row <- tibble(
    term      = "Number of events (rounded to nearest 5)",
    HR        = NA_real_,
    std.error = NA_real_,
    statistic = NA_real_,
    p_value   = NA_real_,
    lower_CI  = NA_real_,
    upper_CI  = NA_real_,
    extra_info = as.character(rounded_events)
  )
  
  out <- bind_rows(events_row, res)
  write.csv(out, out_path, row.names = FALSE)
  out
}

# Choose your output directory
out_dir <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Sensitivity analysis_3months_allcause"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pooled_unadj_tbl <- tidy_pooled(
  pool_unadj,
  allcause_analyses,
  "unadj",
  file.path(out_dir, "allcause_90d_unadjusted_pooled.csv")
)



pooled_full_tbl <- tidy_pooled(
  pool_full,
  allcause_analyses,
  "full",
  file.path(out_dir, "allcause_90d_fully_adjusted_pooled.csv")
)




