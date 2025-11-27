

#Option to pull data or refresh from data cleaning
clean_cohort <- read.csv("clean_cohort.csv")
View(clean_cohort)


install.packages("mice")      # if not already installed
install.packages("survival")  # if you're using Cox models
library(mice)
library(survival)
library(dplyr)
library(lubridate)



# Create survival object for time to cardiac rehab exposure while accounting for censoring due to death or the end of follow-up (same as per MICE - see MICE code)
clean_cohort <- clean_cohort %>%
  mutate(
    rehab_status = ifelse(is.na(exp_cardiac_rehab_6m), 0, ifelse(exp_cardiac_rehab_6m == 1, 1, 0)),  # 1 for exposed, 0 for censored
    cardiac_rehab_time = cardiac_rehab_follow_up_days        # Time in days (time to rehab or censoring)
  )



##convert all relevant date variables to Date format 
clean_cohort <- clean_cohort %>%
  mutate(
    tavi_cips_disdate_plus180days = as.Date(tavi_cips_disdate_plus180days),
    date_of_death = as.Date(date_of_death),
    out_readm_hf_date = as.Date(out_readm_hf_date),
    out_readm_all_cause_date = as.Date(out_readm_all_cause_date),
    out_readm_non_cvd_date = as.Date(out_readm_non_cvd_date)
  )

######Create survival object for time to HF rehosp while accounting for censoring due to death or the end of follow-up (same as per MICE - see MICE code)
##Explanation:
# Define study start and end dates
study_start_new <- clean_cohort$tavi_cips_disdate_plus180days
study_end <- as.Date("2024-03-31")



# Add necessary variables (same as per MICE - see MICE code)
clean_cohort <- clean_cohort %>%
  mutate(
    # Days since study start (calendar time)
    days_since_tavi_dis_plus180 = as.numeric(tavi_cips_disdate_plus180days - study_start_new),
    
    # Time of HF rehospitalization or censoring (death or end of study period)
    hf_rehosp_censor_time = pmin(as.numeric(out_readm_hf_date - study_start_new),         # Days to HF readmission
                                 as.numeric(date_of_death - study_start_new),                # Days to death
                                 as.numeric(study_end - study_start_new), na.rm = TRUE),     # Days to study end
    
    # Rehab exposure indicator
    rehab_indicator = ifelse(is.na(exp_cardiac_rehab_6m) | exp_cardiac_rehab_6m == 0, 0, 1),
    
    # Rehospitalization indicator (1 for HF rehospitalization, 0 for censoring)
    hf_rehosp_indicator = ifelse(is.na(out_readm_hf_flag) | out_readm_hf_flag == 0, 0, 1)
    
  )




#####Create survival object for time to ALL_CAUSE rehosp while accounting for censoring due to death or the end of follow-up (same as per MICE - see MICE code)


# Add necessary variables
clean_cohort <- clean_cohort %>%
  mutate(
    # Days since study start (calendar time)
    days_since_tavi_dis_plus180 = as.numeric(tavi_cips_disdate_plus180days - study_start_new),
    
    # Time of HF rehospitalization or censoring (death or end of study period)
    all_cause_rehosp_censor_time = pmin(as.numeric(out_readm_all_cause_date - study_start_new),         # Days to all cause readmission
                                        as.numeric(date_of_death - study_start_new),                # Days to death
                                        as.numeric(study_end - study_start_new), na.rm = TRUE),     # Days to study end
    
    # Rehab exposure indicator
    rehab_indicator = ifelse(is.na(exp_cardiac_rehab_6m) | exp_cardiac_rehab_6m == 0, 0, 1),
    
    # Rehospitalization indicator (1 for all cause rehospitalization, 0 for censoring)
    all_cause_rehosp_indicator = ifelse(is.na(out_readm_all_cause_flag) | out_readm_all_cause_flag == 0, 0, 1)
    
  )





#####Create survival object for time to NON-CARDIOVASCULAR rehosp while accounting for censoring due to death or the end of follow-up (same as per MICE - see MICE code)


# Add necessary variables
clean_cohort <- clean_cohort %>%
  mutate(
    # Days since study start (calendar time)
    days_since_tavi_dis_plus180 = as.numeric(tavi_cips_disdate_plus180days - study_start_new),
    
    # Time of non_CVD rehospitalization or censoring (death or end of study period)
    non_cvd_rehosp_censor_time = pmin(as.numeric(out_readm_non_cvd_date - study_start_new),         # Days to non_cvd readmission
                                      as.numeric(date_of_death - study_start_new),                # Days to death
                                      as.numeric(study_end - study_start_new), na.rm = TRUE),     # Days to study end
    
    # Rehab exposure indicator
    rehab_indicator = ifelse(is.na(exp_cardiac_rehab_6m) | exp_cardiac_rehab_6m == 0, 0, 1),
    
    # non_cvd Rehospitalization indicator (1 for non-CVD rehospitalization, 0 for censoring)
    non_cvd_rehosp_indicator = ifelse(is.na(out_readm_non_cvd_flag) | out_readm_non_cvd_flag == 0, 0, 1)
    
  )





#####Create survival object for time to MORTALITY while accounting for censoring due to death or the end of follow-up (same as per MICE - see MICE code)


# Add necessary variables
clean_cohort <- clean_cohort %>%
  mutate(
    # Days since study start (calendar time)
    days_since_tavi_dis_plus180 = as.numeric(tavi_cips_disdate_plus180days - study_start_new),
    
    # Time to death or censoring (end of study period)
    all_cause_mortality_censor_time = pmin (as.numeric(date_of_death - study_start_new),            # Days to death
                                            as.numeric(study_end - study_start_new), na.rm = TRUE),     # Days to study end
    
    # Rehab exposure indicator
    rehab_indicator = ifelse(is.na(exp_cardiac_rehab_6m) | exp_cardiac_rehab_6m == 0, 0, 1),
    
    # mortality indicator (1 for all cause mortality, 0 for censoring)
    all_cause_mortality_indicator = ifelse(!is.na(date_of_death), 1, 0)
  )








#########################################################################################################################################################################################################
##                                                                                START FINE-GREY FROM HERE DOWN - ALWAYS RESET THE DATA BEFORE EACH MODEL      
#########################################################################################################################################################################################################


############################################
# NEW: helper to pool Fine-Gray (crr) models   Helper for pooling Fine-Gray models (add once) - do this before every new model below after resetting the data
############################################


library(dplyr)

# Rubin pooling for crr models
library(cmprsk)

pool_crr <- function(fit_list) {
  m <- length(fit_list)
  if (m < 2L) stop("Need at least 2 imputations")
  
  # coefficient matrix: m x p
  coef_mat <- t(sapply(fit_list, function(f) as.numeric(f$coef)))
  
  # variance matrix: m x p (diagonal of var-cov)
  var_mat  <- t(sapply(fit_list, function(f) diag(f$var)))
  
  term_names <- names(fit_list[[1]]$coef)
  
  # Rubin's rules
  Q_bar <- colMeans(coef_mat)           # pooled log-HR
  U_bar <- colMeans(var_mat)           # within-imputation variance
  B     <- apply(coef_mat, 2, var)     # between-imputation variance
  
  T_var <- U_bar + (1 + 1/m) * B       # total variance
  se    <- sqrt(T_var)
  
  # degrees of freedom
  df <- (m - 1) * (1 + U_bar / ((1 + 1/m) * B))^2
  
  t_val <- Q_bar / se
  p_val <- 2 * pt(abs(t_val), df = df, lower.tail = FALSE)
  
  HR        <- exp(Q_bar)
  lower_CI  <- exp(Q_bar - 1.96 * se)
  upper_CI  <- exp(Q_bar + 1.96 * se)
  
  data.frame(
    term        = term_names,
    logHR       = Q_bar,
    se          = se,
    df          = df,
    HR          = HR,
    lower_95_CI = lower_CI,
    upper_95_CI = upper_CI,
    p.value     = p_val,
    row.names   = NULL
  )
}



##########################################################################################################################################################
##                                                       1. Heart Failure rehospitalisation - add competing-risk time & status + Fine-Gray
##########################################################################################################################################################




##FULLY ADJUSTED MODEL: HEART FAILURE REHOSP *************************************RESTART THE DATA from TOP of sheet****************************************
###############################
# NEW: covariate matrix builder for HF Fine-Gray
###############################
hf_covars <- function(data) {
  model.matrix(
    ~ rehab_indicator + age_groups + sex_code + IMD_2019_QUINTILES + SMOKING_STATUS + 
      year_tavi_discharge + TAVI_procedural_complications + ethnicity_5_group + 
      BMI_category + KATZ_INDEX_CAT +  COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
      previous_MI_time_grouped + binary_known_prev_cardiac_surgery +
      binary_known_poor_mobility + region_name + binary_Known_DIABETES + 
      NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + LV_FUNCTION + MITRAL_REGURGITATION,
    data = data
  )[ , -1, drop = FALSE]  # remove the intercept column
}





##########################
# NEW: competing-risk vars for HF rehosp  (a) add the competing-risk time/status variables, then (b) add a Fine-Gray fit after your Cox models.
##########################

##rebuild completed_datasets and add the date columns + competing-risk vars

library(dplyr)
library(survival)
library(mice)
library(cmprsk)


# Load imputation object and recreate completed datasets
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")
completed_datasets <- lapply(1:5, function(i) complete(imp, i))


study_end <- as.Date("2024-03-31")

# Add original dates + competing-risk time/status into EACH completed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  # bring original date variables back in by matching person_id
  df$tavi_cips_disdate_plus180days <- clean_cohort$tavi_cips_disdate_plus180days[
    match(df$person_id, clean_cohort$person_id)
  ]  
    df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[  # Add days_since_tavi_dis_plus180 back to each imputed dataset ## Ensure it has been added to clean_cohort (see code at the start of this script)
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_hf_date <- clean_cohort$out_readm_hf_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$date_of_death <- clean_cohort$date_of_death[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_all_cause_date <- clean_cohort$out_readm_all_cause_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_non_cvd_date <- clean_cohort$out_readm_non_cvd_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  
  
  df <- df %>%
    mutate(
      t0      = tavi_cips_disdate_plus180days,
      t_hf    = as.numeric(out_readm_hf_date   - t0),
      t_death = as.numeric(date_of_death       - t0),
      t_end   = as.numeric(study_end           - t0),
      # time to first of HF rehosp, death, or end of follow-up
      time_hf_cr = pmin(t_hf, t_death, t_end, na.rm = TRUE),
      # 0 = censored, 1 = HF rehosp, 2 = death before HF rehosp
      status_hf_cr = dplyr::case_when(
        !is.na(t_hf)    & t_hf <= t_death & t_hf <= t_end ~ 1L,
        !is.na(t_death) & (is.na(t_hf) | t_death < t_hf) & t_death <= t_end ~ 2L,
        TRUE ~ 0L
      )
    )
  return(df)
})


# Sanity check
str(completed_datasets[[1]][, c("time_hf_cr", "status_hf_cr")])

##Key point: now every completed dataset has time_hf_cr and status_hf_cr before we do any filtering.


### Need to filter out events within 180 days for each of the outcomes (do each outcome separately)

# Pre-filter all datasets first
filtered_datasets_hf_rehosp <- lapply(completed_datasets, function(df) {
  df[df$hf_rehosp_censor_time > df$days_since_tavi_dis_plus180, ]
})


# Check you still have rows and the new variables are present
##should see a positive row count and both columns as numeric/integer.
nrow(filtered_datasets_hf_rehosp[[1]])
str(filtered_datasets_hf_rehosp[[1]][, c("time_hf_cr", "status_hf_cr")])
# n=2190 (rounded to nearets 5) patients filtered - this aligns with cause specific cox models. Confirmed correct 

#run the Fine-Gray models and pool
#Then you can use time_hf_cr and status_hf_cr in the Fine-Gray model
# Fit Fine-Gray model in each imputed dataset
fg_models_hf <- lapply(filtered_datasets_hf_rehosp, function(data) {    
  crr(
    ftime   = data$time_hf_cr,
    fstatus = data$status_hf_cr,  # 0=censor, 1=HF rehosp, 2=death
    cov1    = hf_covars(data),
    failcode = 1,
    cencode  = 0
  )
})

lapply(fg_models_hf, summary) ##extract the summary from each model

#Extract term, HR, CI, p for each model
extract_fg <- function(fg_obj, imp_number) {
  
  s <- summary(fg_obj)
  coef_df <- as.data.frame(s$coef)
  
  # Add term names
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  # Compute HR and CI
  coef_df$HR          <- exp(coef_df$coef)
  coef_df$lower_95_CI <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
  coef_df$upper_95_CI <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
  
  # Keep requested columns + imputation number
  coef_df <- coef_df[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p-value")]
  coef_df$imputation <- imp_number
  
  return(coef_df)
}



#Apply to all 5 imputed models
all_fg_results <- do.call(
  rbind,
  lapply(seq_along(fg_models_hf), function(i) extract_fg(fg_models_hf[[i]], i))
)

#Save to CSV
write.csv(
  all_fg_results,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/hf_rehosp/hf_finegray_all_imputations_fully_adjusted.csv",
  row.names = FALSE
)



# Pool Fine-Gray models
fg_pooled_hf <- pool_crr(fg_models_hf)

# Keep term, HR, CI, p-value
fg_hf_df <- fg_pooled_hf[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
fg_hf_df[, 2:5] <- round(fg_hf_df[, 2:5], 3)

# Save to CSV
write.csv(
  fg_hf_df,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/hf_rehosp/hf_rehosp_finegray_pooled_fully_adjusted.csv",
  row.names = FALSE
)






##DEMOGRAPHICS ADJUSTED MODEL: HEART FAILURE REHOSP *************************************RESTART THE DATA from TOP of sheet****************************************
###############################
# NEW: covariate matrix builder for HF Fine-Gray
###############################
hf_covars <- function(data) {
  model.matrix(
    ~ rehab_indicator + age_groups + sex_code + IMD_2019_QUINTILES + SMOKING_STATUS + 
      ethnicity_5_group + region_name,
    data = data
  )[ , -1, drop = FALSE]  # remove the intercept column
}




##########################
# NEW: competing-risk vars for HF rehosp  (a) add the competing-risk time/status variables, then (b) add a Fine-Gray fit after your Cox models.
##########################

##rebuild completed_datasets and add the date columns + competing-risk vars

library(dplyr)
library(survival)
library(mice)
library(cmprsk)


# Load imputation object and recreate completed datasets
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")
completed_datasets <- lapply(1:5, function(i) complete(imp, i))


study_end <- as.Date("2024-03-31")

# Add original dates + competing-risk time/status into EACH completed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  # bring original date variables back in by matching person_id
  df$tavi_cips_disdate_plus180days <- clean_cohort$tavi_cips_disdate_plus180days[
    match(df$person_id, clean_cohort$person_id)
  ]  
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[  # Add days_since_tavi_dis_plus180 back to each imputed dataset ## Ensure it has been added to clean_cohort (see code at the start of this script)
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_hf_date <- clean_cohort$out_readm_hf_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$date_of_death <- clean_cohort$date_of_death[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_all_cause_date <- clean_cohort$out_readm_all_cause_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_non_cvd_date <- clean_cohort$out_readm_non_cvd_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  
  
  df <- df %>%
    mutate(
      t0      = tavi_cips_disdate_plus180days,
      t_hf    = as.numeric(out_readm_hf_date   - t0),
      t_death = as.numeric(date_of_death       - t0),
      t_end   = as.numeric(study_end           - t0),
      # time to first of HF rehosp, death, or end of follow-up
      time_hf_cr = pmin(t_hf, t_death, t_end, na.rm = TRUE),
      # 0 = censored, 1 = HF rehosp, 2 = death before HF rehosp
      status_hf_cr = dplyr::case_when(
        !is.na(t_hf)    & t_hf <= t_death & t_hf <= t_end ~ 1L,
        !is.na(t_death) & (is.na(t_hf) | t_death < t_hf) & t_death <= t_end ~ 2L,
        TRUE ~ 0L
      )
    )
  return(df)
})


# Sanity check
str(completed_datasets[[1]][, c("time_hf_cr", "status_hf_cr")])

##Key point: now every completed dataset has time_hf_cr and status_hf_cr before we do any filtering.


### Need to filter out events within 180 days for each of the outcomes (do each outcome separately)

# Pre-filter all datasets first
filtered_datasets_hf_rehosp <- lapply(completed_datasets, function(df) {
  df[df$hf_rehosp_censor_time > df$days_since_tavi_dis_plus180, ]
})


# Check you still have rows and the new variables are present
##should see a positive row count and both columns as numeric/integer.
nrow(filtered_datasets_hf_rehosp[[1]])
str(filtered_datasets_hf_rehosp[[1]][, c("time_hf_cr", "status_hf_cr")])
# n=2190 (rounded to nearets 5) patients filtered - this aligns with cause specific cox models. Confirmed correct 

#run the Fine-Gray models and pool
#Then you can use time_hf_cr and status_hf_cr in the Fine-Gray model
# Fit Fine-Gray model in each imputed dataset
fg_models_hf <- lapply(filtered_datasets_hf_rehosp, function(data) {    
  crr(
    ftime   = data$time_hf_cr,
    fstatus = data$status_hf_cr,  # 0=censor, 1=HF rehosp, 2=death
    cov1    = hf_covars(data),
    failcode = 1,
    cencode  = 0
  )
})

lapply(fg_models_hf, summary) ##extract the summary from each model

#Extract term, HR, CI, p for each model
extract_fg <- function(fg_obj, imp_number) {
  
  s <- summary(fg_obj)
  coef_df <- as.data.frame(s$coef)
  
  # Add term names
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  # Compute HR and CI
  coef_df$HR          <- exp(coef_df$coef)
  coef_df$lower_95_CI <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
  coef_df$upper_95_CI <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
  
  # Keep requested columns + imputation number
  coef_df <- coef_df[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p-value")]
  coef_df$imputation <- imp_number
  
  return(coef_df)
}



#Apply to all 5 imputed models
all_fg_results_demographics_adjusted <- do.call(
  rbind,
  lapply(seq_along(fg_models_hf), function(i) extract_fg(fg_models_hf[[i]], i))
)

#Save to CSV
write.csv(
  all_fg_results_demographics_adjusted,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/hf_rehosp/hf_finegray_all_imputations_demographics_adjusted.csv",
  row.names = FALSE
)





# Pool Fine-Gray models
fg_pooled_hf <- pool_crr(fg_models_hf)

# Keep term, HR, CI, p-value
fg_hf_df <- fg_pooled_hf[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
fg_hf_df[, 2:5] <- round(fg_hf_df[, 2:5], 3)

# Save to CSV
write.csv(
  fg_hf_df,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/hf_rehosp/hf_rehosp_finegray_pooled_demographics_adjusted.csv",
  row.names = FALSE
)











##Heart FAILURE UNADJUSTED MODEL  *****************************RESTART THE DATA from TOP of sheet******************
###############################
# NEW: covariate matrix builder for HF Fine-Gray
###############################
hf_covars_unadjusted <- function(data) {
  model.matrix(
    ~ rehab_indicator,
    data = data
  )[ , -1, drop = FALSE]  # remove the intercept column
}




##########################
# NEW: competing-risk vars for HF rehosp  (a) add the competing-risk time/status variables, then (b) add a Fine-Gray fit after your Cox models.
##########################

##rebuild completed_datasets and add the date columns + competing-risk vars

library(dplyr)
library(survival)
library(mice)
library(cmprsk)


# Load imputation object and recreate completed datasets
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")
completed_datasets <- lapply(1:5, function(i) complete(imp, i))


study_end <- as.Date("2024-03-31")

# Add original dates + competing-risk time/status into EACH completed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  # bring original date variables back in by matching person_id
  df$tavi_cips_disdate_plus180days <- clean_cohort$tavi_cips_disdate_plus180days[
    match(df$person_id, clean_cohort$person_id)
  ]  
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[  # Add days_since_tavi_dis_plus180 back to each imputed dataset ## Ensure it has been added to clean_cohort (see code at the start of this script)
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_hf_date <- clean_cohort$out_readm_hf_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$date_of_death <- clean_cohort$date_of_death[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_all_cause_date <- clean_cohort$out_readm_all_cause_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_non_cvd_date <- clean_cohort$out_readm_non_cvd_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  
  
  df <- df %>%
    mutate(
      t0      = tavi_cips_disdate_plus180days,
      t_hf    = as.numeric(out_readm_hf_date   - t0),
      t_death = as.numeric(date_of_death       - t0),
      t_end   = as.numeric(study_end           - t0),
      # time to first of HF rehosp, death, or end of follow-up
      time_hf_cr = pmin(t_hf, t_death, t_end, na.rm = TRUE),
      # 0 = censored, 1 = HF rehosp, 2 = death before HF rehosp
      status_hf_cr = dplyr::case_when(
        !is.na(t_hf)    & t_hf <= t_death & t_hf <= t_end ~ 1L,
        !is.na(t_death) & (is.na(t_hf) | t_death < t_hf) & t_death <= t_end ~ 2L,
        TRUE ~ 0L
      )
    )
  return(df)
})


# Sanity check
str(completed_datasets[[1]][, c("time_hf_cr", "status_hf_cr")])

##Key point: now every completed dataset has time_hf_cr and status_hf_cr before we do any filtering.


### Need to filter out events within 180 days for each of the outcomes (do each outcome separately)
# Pre-filter all datasets first
filtered_datasets_hf_rehosp <- lapply(completed_datasets, function(df) {
  df[df$hf_rehosp_censor_time > df$days_since_tavi_dis_plus180, ]
})


# Check you still have rows and the new variables are present
##should see a positive row count and both columns as numeric/integer.
nrow(filtered_datasets_hf_rehosp[[1]])
str(filtered_datasets_hf_rehosp[[1]][, c("time_hf_cr", "status_hf_cr")])
# n=2190 (rounded to nearets 5) patients filtered - this aligns with cause specific cox models. Confirmed correct 

#run the Fine-Gray models and pool
#Then you can use time_hf_cr and status_hf_cr in the Fine-Gray model
# Fit Fine-Gray model in each imputed dataset
fg_models_hf <- lapply(filtered_datasets_hf_rehosp, function(data) {    
  crr(
    ftime   = data$time_hf_cr,
    fstatus = data$status_hf_cr,  # 0=censor, 1=HF rehosp, 2=death
    cov1    = hf_covars_unadjusted(data),
    failcode = 1,
    cencode  = 0
  )
})

lapply(fg_models_hf, summary) ##extract the summary from each model

#Extract term, HR, CI, p for each model
extract_fg <- function(fg_obj, imp_number) {
  
  s <- summary(fg_obj)
  coef_df <- as.data.frame(s$coef)
  
  # Add term names
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  # Compute HR and CI
  coef_df$HR          <- exp(coef_df$coef)
  coef_df$lower_95_CI <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
  coef_df$upper_95_CI <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
  
  # Keep requested columns + imputation number
  coef_df <- coef_df[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p-value")]
  coef_df$imputation <- imp_number
  
  return(coef_df)
}



#Apply to all 5 imputed models
all_fg_results_unadjusted <- do.call(
  rbind,
  lapply(seq_along(fg_models_hf), function(i) extract_fg(fg_models_hf[[i]], i))
)

#Save to CSV
write.csv(
  all_fg_results_unadjusted,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/hf_rehosp/hf_finegray_all_imputations_unadjusted.csv",
  row.names = FALSE
)




# Pool Fine-Gray models
fg_pooled_hf <- pool_crr(fg_models_hf)

# Keep term, HR, CI, p-value
fg_hf_df <- fg_pooled_hf[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
fg_hf_df[, 2:5] <- round(fg_hf_df[, 2:5], 3)

# Save to CSV
write.csv(
  fg_hf_df,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/hf_rehosp/hf_rehosp_finegray_pooled_unadjusted.csv",
  row.names = FALSE
)






###RESET THE DATA BEFORE RUNNING THIS
##########################################################################################################################################################
##                                               2. All-cause rehospitalisation - add competing-risk time & status + Fine-Gray
##########################################################################################################################################################



##FULLY ADJUSTED MODEL ALL-CAUSE REHOSP  *************************************RESTART THE DATA from TOP of sheet****************************************
###############################
# NEW: covariate matrix builder for HF Fine-Gray
###############################
hf_covars <- function(data) {
  model.matrix(
    ~ rehab_indicator + age_groups + sex_code + IMD_2019_QUINTILES + SMOKING_STATUS + 
      year_tavi_discharge + TAVI_procedural_complications + ethnicity_5_group + 
      BMI_category + KATZ_INDEX_CAT +  COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
      previous_MI_time_grouped + binary_known_prev_cardiac_surgery +
      binary_known_poor_mobility + region_name + binary_Known_DIABETES + 
      NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + LV_FUNCTION + MITRAL_REGURGITATION,
    data = data
  )[ , -1, drop = FALSE]  # remove the intercept column
}





library(dplyr)
library(survival)
library(mice)
library(cmprsk)

# Load imputation object and recreate completed datasets
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")
completed_datasets <- lapply(1:5, function(i) complete(imp, i))

study_end <- as.Date("2024-03-31")

# Add original date variables + competing-risk time/status into EACH completed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  # bring original date variables back in by matching person_id
  df$tavi_cips_disdate_plus180days <- clean_cohort$tavi_cips_disdate_plus180days[
    match(df$person_id, clean_cohort$person_id)
  ]  
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[  # Add days_since_tavi_dis_plus180 back to each imputed dataset ## Ensure it has been added to clean_cohort (see code at the start of this script)
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_hf_date <- clean_cohort$out_readm_hf_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$date_of_death <- clean_cohort$date_of_death[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_all_cause_date <- clean_cohort$out_readm_all_cause_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_non_cvd_date <- clean_cohort$out_readm_non_cvd_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  
  df <- df %>%
    mutate(
      t0       = tavi_cips_disdate_plus180days,
      t_all    = as.numeric(out_readm_all_cause_date - t0),
      t_death  = as.numeric(date_of_death           - t0),
      t_end    = as.numeric(study_end               - t0),
      # time to first of all-cause rehosp, death, or end of follow-up
      time_all_cr = pmin(t_all, t_death, t_end, na.rm = TRUE),
      # 0 = censored, 1 = all-cause rehosp, 2 = death before rehosp
      status_all_cr = dplyr::case_when(
        !is.na(t_all)   & t_all <= t_death & t_all <= t_end ~ 1L,
        !is.na(t_death) & (is.na(t_all) | t_death < t_all) & t_death <= t_end ~ 2L,
        TRUE ~ 0L
      )
    )
  
  return(df)
})

# Quick sanity check
str(completed_datasets[[1]][, c("time_all_cr", "status_all_cr")])

# Filter to valid times > 0  
filtered_datasets_all_rehosp <- lapply(completed_datasets, function(df) {
  df[df$time_all_cr > 0 & !is.na(df$status_all_cr), ]
})



# Check still have rows and the new variables are present
##should see a positive row count and both columns as numeric/integer.
nrow(filtered_datasets_all_rehosp[[1]])
str(filtered_datasets_all_rehosp[[1]][, c("time_all_cr", "status_all_cr")])
##NOTES, n=8915 (rounded to nearest 5)  patients filtred out which exactly aligns with previous cox-model analyses for all cause rehosps. Confirmed correct


# Fit Fine-Gray model in each imputed dataset (same covariates as HF via hf_covars)
fg_models_all <- lapply(filtered_datasets_all_rehosp, function(data) {
  crr(
    ftime   = data$time_all_cr,
    fstatus = data$status_all_cr,   # 0=censor, 1=all-cause rehosp, 2=death
    cov1    = hf_covars(data),
    failcode = 1,
    cencode  = 0
  )
})






lapply(fg_models_all, summary) ##extract the summary from each model

#Extract term, HR, CI, p for each model
extract_fg <- function(fg_obj, imp_number) {
  
  s <- summary(fg_obj)
  coef_df <- as.data.frame(s$coef)
  
  # Add term names
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  # Compute HR and CI
  coef_df$HR          <- exp(coef_df$coef)
  coef_df$lower_95_CI <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
  coef_df$upper_95_CI <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
  
  # Keep requested columns + imputation number
  coef_df <- coef_df[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p-value")]
  coef_df$imputation <- imp_number
  
  return(coef_df)
}



#Apply to all 5 imputed models
all_fg_results_fully_adjusted <- do.call(
  rbind,
  lapply(seq_along(fg_models_all), function(i) extract_fg(fg_models_all[[i]], i))
)

#Save to CSV
write.csv(
  all_fg_results_fully_adjusted,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/all_cause_rehosp/all_cause_finegray_all_imputations_fully_adjusted.csv",
  row.names = FALSE
)



##POOL RESULTS AFTER 

# Pool Fine-Gray models
fg_pooled_all <- pool_crr(fg_models_all)

# Keep term, HR, CI, p-value
fg_all_df <- fg_pooled_all[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
fg_all_df[, 2:5] <- round(fg_all_df[, 2:5], 3)

# Save to CSV
write.csv(
  fg_all_df,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/all_cause_rehosp/all_cause_rehosp_finegray_pooled_fully_adjusted.csv",
  row.names = FALSE
)




## OPTIONAL: CIF + Gray's test for all-cause rehosp (first imputed dataset)
d1_all <- filtered_datasets_all_rehosp[[1]]

ci_all <- cuminc(
  ftime   = d1_all$time_all_cr,
  fstatus = d1_all$status_all_cr,      # 0=censor, 1=all-cause rehosp, 2=death
  group   = d1_all$rehab_indicator     # 0 = no rehab, 1 = rehab
)

ci_all     # includes Gray's test
# plot(ci_all)  # uncomment to see CIF curves










##################################################################
################      UNADJUSTED MODEL ALL-CAUSE REHOSP *************************************RESTART THE DATA from TOP of sheet****************************************
###############################
# NEW: covariate matrix builder for HF Fine-Gray
###############################
hf_covars_unadjusted <- function(data) {
  model.matrix(
    ~ rehab_indicator,
    data = data
  )[ , -1, drop = FALSE]  # remove the intercept column
}




library(dplyr)
library(survival)
library(mice)
library(cmprsk)

# Load imputation object and recreate completed datasets
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")
completed_datasets <- lapply(1:5, function(i) complete(imp, i))

study_end <- as.Date("2024-03-31")

# Add original date variables + competing-risk time/status into EACH completed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  # bring original date variables back in by matching person_id
  df$tavi_cips_disdate_plus180days <- clean_cohort$tavi_cips_disdate_plus180days[
    match(df$person_id, clean_cohort$person_id)
  ]  
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[  # Add days_since_tavi_dis_plus180 back to each imputed dataset ## Ensure it has been added to clean_cohort (see code at the start of this script)
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_hf_date <- clean_cohort$out_readm_hf_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$date_of_death <- clean_cohort$date_of_death[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_all_cause_date <- clean_cohort$out_readm_all_cause_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_non_cvd_date <- clean_cohort$out_readm_non_cvd_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  
  df <- df %>%
    mutate(
      t0       = tavi_cips_disdate_plus180days,
      t_all    = as.numeric(out_readm_all_cause_date - t0),
      t_death  = as.numeric(date_of_death           - t0),
      t_end    = as.numeric(study_end               - t0),
      # time to first of all-cause rehosp, death, or end of follow-up
      time_all_cr = pmin(t_all, t_death, t_end, na.rm = TRUE),
      # 0 = censored, 1 = all-cause rehosp, 2 = death before rehosp
      status_all_cr = dplyr::case_when(
        !is.na(t_all)   & t_all <= t_death & t_all <= t_end ~ 1L,
        !is.na(t_death) & (is.na(t_all) | t_death < t_all) & t_death <= t_end ~ 2L,
        TRUE ~ 0L
      )
    )
  
  return(df)
})

# Quick sanity check
str(completed_datasets[[1]][, c("time_all_cr", "status_all_cr")])

# Filter to valid times > 0  ##DOUBLE CHECK THIS IS ALL-CAUSE REHOSP AFTER 180 DAYS _ confirmed its the same: summary(completed_datasets[[1]]$time_all_cr) = summary(clean_cohort$all_cause_rehosp_censor_time)
filtered_datasets_all_rehosp <- lapply(completed_datasets, function(df) {
  df[df$time_all_cr > 0 & !is.na(df$status_all_cr), ]
})



# Check still have rows and the new variables are present
##should see a positive row count and both columns as numeric/integer.
nrow(filtered_datasets_all_rehosp[[1]])
str(filtered_datasets_all_rehosp[[1]][, c("time_all_cr", "status_all_cr")])
##NOTES, n=8915 (rounded to nearest 5)  patients filtred out which exactly aligns with previous cox-model analyses for all cause rehosps. Confirmed correct


# Fit Fine-Gray model in each imputed dataset (same covariates as HF via hf_covars)
fg_models_all <- lapply(filtered_datasets_all_rehosp, function(data) {
  crr(
    ftime   = data$time_all_cr,
    fstatus = data$status_all_cr,   # 0=censor, 1=all-cause rehosp, 2=death
    cov1    = hf_covars_unadjusted(data),
    failcode = 1,
    cencode  = 0
  )
})






lapply(fg_models_all, summary) ##extract the summary from each model

#Extract term, HR, CI, p for each model
extract_fg <- function(fg_obj, imp_number) {
  
  s <- summary(fg_obj)
  coef_df <- as.data.frame(s$coef)
  
  # Add term names
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  # Compute HR and CI
  coef_df$HR          <- exp(coef_df$coef)
  coef_df$lower_95_CI <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
  coef_df$upper_95_CI <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
  
  # Keep requested columns + imputation number
  coef_df <- coef_df[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p-value")]
  coef_df$imputation <- imp_number
  
  return(coef_df)
}



#Apply to all 5 imputed models
all_fg_results_unadjusted <- do.call(
  rbind,
  lapply(seq_along(fg_models_all), function(i) extract_fg(fg_models_all[[i]], i))
)

#Save to CSV
write.csv(
  all_fg_results_unadjusted,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/all_cause_rehosp/all_cause_finegray_all_imputations_unadjusted.csv",
  row.names = FALSE
)


##POOL RESULTS AFTER 

# Pool Fine-Gray models
fg_pooled_all <- pool_crr(fg_models_all)

# Keep term, HR, CI, p-value
fg_all_df <- fg_pooled_all[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
fg_all_df[, 2:5] <- round(fg_all_df[, 2:5], 3)

# Save to CSV
write.csv(
  fg_all_df,
  "D:/PhotonUser/My Files\Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/all_cause_rehosp/all_cause_rehosp_finegray_pooled_unadjusted.csv",
  row.names = FALSE
)









##DEMOGRAPHICS ADJUSTED MODEL - ALL CAUSE REHOSP *************************************RESTART THE DATA from TOP of sheet****************************************
###############################
# NEW: covariate matrix builder for HF Fine-Gray
###############################
hf_covars <- function(data) {
  model.matrix(
    ~ rehab_indicator + age_groups + sex_code + IMD_2019_QUINTILES + SMOKING_STATUS + 
      ethnicity_5_group + region_name,
    data = data
  )[ , -1, drop = FALSE]  # remove the intercept column
}





library(dplyr)
library(survival)
library(mice)
library(cmprsk)

# Load imputation object and recreate completed datasets
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")
completed_datasets <- lapply(1:5, function(i) complete(imp, i))

study_end <- as.Date("2024-03-31")

# Add original date variables + competing-risk time/status into EACH completed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  # bring original date variables back in by matching person_id
  df$tavi_cips_disdate_plus180days <- clean_cohort$tavi_cips_disdate_plus180days[
    match(df$person_id, clean_cohort$person_id)
  ]  
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[  # Add days_since_tavi_dis_plus180 back to each imputed dataset ## Ensure it has been added to clean_cohort (see code at the start of this script)
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_hf_date <- clean_cohort$out_readm_hf_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$date_of_death <- clean_cohort$date_of_death[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_all_cause_date <- clean_cohort$out_readm_all_cause_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_non_cvd_date <- clean_cohort$out_readm_non_cvd_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  
  df <- df %>%
    mutate(
      t0       = tavi_cips_disdate_plus180days,
      t_all    = as.numeric(out_readm_all_cause_date - t0),
      t_death  = as.numeric(date_of_death           - t0),
      t_end    = as.numeric(study_end               - t0),
      # time to first of all-cause rehosp, death, or end of follow-up
      time_all_cr = pmin(t_all, t_death, t_end, na.rm = TRUE),
      # 0 = censored, 1 = all-cause rehosp, 2 = death before rehosp
      status_all_cr = dplyr::case_when(
        !is.na(t_all)   & t_all <= t_death & t_all <= t_end ~ 1L,
        !is.na(t_death) & (is.na(t_all) | t_death < t_all) & t_death <= t_end ~ 2L,
        TRUE ~ 0L
      )
    )
  
  return(df)
})

# Quick sanity check
str(completed_datasets[[1]][, c("time_all_cr", "status_all_cr")])

# Filter to valid times > 0  
filtered_datasets_all_rehosp <- lapply(completed_datasets, function(df) {
  df[df$time_all_cr > 0 & !is.na(df$status_all_cr), ]
})



# Check still have rows and the new variables are present
##should see a positive row count and both columns as numeric/integer.
nrow(filtered_datasets_all_rehosp[[1]])
str(filtered_datasets_all_rehosp[[1]][, c("time_all_cr", "status_all_cr")])
##NOTES, n=8915 (rounded to nearest 5)  patients filtred out which exactly aligns with previous cox-model analyses for all cause rehosps. Confirmed correct


# Fit Fine-Gray model in each imputed dataset (same covariates as HF via hf_covars)
fg_models_all <- lapply(filtered_datasets_all_rehosp, function(data) {
  crr(
    ftime   = data$time_all_cr,
    fstatus = data$status_all_cr,   # 0=censor, 1=all-cause rehosp, 2=death
    cov1    = hf_covars(data),
    failcode = 1,
    cencode  = 0
  )
})



lapply(fg_models_all, summary) ##extract the summary from each model

#Extract term, HR, CI, p for each model
extract_fg <- function(fg_obj, imp_number) {
  
  s <- summary(fg_obj)
  coef_df <- as.data.frame(s$coef)
  
  # Add term names
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  # Compute HR and CI
  coef_df$HR          <- exp(coef_df$coef)
  coef_df$lower_95_CI <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
  coef_df$upper_95_CI <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
  
  # Keep requested columns + imputation number
  coef_df <- coef_df[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p-value")]
  coef_df$imputation <- imp_number
  
  return(coef_df)
}



#Apply to all 5 imputed models
all_fg_results_demographics_adjusted <- do.call(
  rbind,
  lapply(seq_along(fg_models_all), function(i) extract_fg(fg_models_all[[i]], i))
)

#Save to CSV
write.csv(
  all_fg_results_demographics_adjusted,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/all_cause_rehosp/all_cause_finegray_all_imputations_demographics_adjusted.csv",
  row.names = FALSE
)


##POOL RESULTS AFTER 

# Pool Fine-Gray models
fg_pooled_all <- pool_crr(fg_models_all)

# Keep term, HR, CI, p-value
fg_all_df <- fg_pooled_all[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
fg_all_df[, 2:5] <- round(fg_all_df[, 2:5], 3)

# Save to CSV
write.csv(
  fg_all_df,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/all_cause_rehosp/all_cause_rehosp_finegray_pooled_demographics_adjusted.csv",
  row.names = FALSE
)












###RESET THE DATA BEFORE RUNNING THIS
##########################################################################################################################################################
##                                               3. non-CVD rehospitalisation - add competing-risk time & status + Fine-Gray
##########################################################################################################################################################



##FULLY ADJUSTED MODEL Non-CVD REHOSP  *************************************RESTART THE DATA from TOP of sheet****************************************
###############################
# NEW: covariate matrix builder for HF Fine-Gray
###############################
hf_covars <- function(data) {
  model.matrix(
    ~ rehab_indicator + age_groups + sex_code + IMD_2019_QUINTILES + SMOKING_STATUS + 
      year_tavi_discharge + TAVI_procedural_complications + ethnicity_5_group + 
      BMI_category + KATZ_INDEX_CAT +  COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
      previous_MI_time_grouped + binary_known_prev_cardiac_surgery +
      binary_known_poor_mobility + region_name + binary_Known_DIABETES + 
      NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + LV_FUNCTION + MITRAL_REGURGITATION,
    data = data
  )[ , -1, drop = FALSE]  # remove the intercept column
}



################################################################################
## Non-CVD rehospitalisation - add competing-risk time & status + Fine-Gray
################################################################################

library(dplyr)
library(survival)
library(mice)
library(cmprsk)

# Load imputation object and recreate completed datasets
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")
completed_datasets <- lapply(1:5, function(i) complete(imp, i))

study_end <- as.Date("2024-03-31")

# Add original date variables + competing-risk time/status into EACH completed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  
  df$tavi_cips_disdate_plus180days <- clean_cohort$tavi_cips_disdate_plus180days[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_non_cvd_date <- clean_cohort$out_readm_non_cvd_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$date_of_death <- clean_cohort$date_of_death[
    match(df$person_id, clean_cohort$person_id)
  ]
  
  df <- df %>%
    mutate(
      t0        = tavi_cips_disdate_plus180days,
      t_noncvd  = as.numeric(out_readm_non_cvd_date - t0),
      t_death   = as.numeric(date_of_death         - t0),
      t_end     = as.numeric(study_end             - t0),
      # time to first of non-CVD rehosp, death, or end of follow-up
      time_noncvd_cr = pmin(t_noncvd, t_death, t_end, na.rm = TRUE),
      # 0 = censored, 1 = non-CVD rehosp, 2 = death before rehosp
      status_noncvd_cr = dplyr::case_when(
        !is.na(t_noncvd) & t_noncvd <= t_death & t_noncvd <= t_end ~ 1L,
        !is.na(t_death)  & (is.na(t_noncvd) | t_death < t_noncvd) & t_death <= t_end ~ 2L,
        TRUE ~ 0L
      )
    )
  
  return(df)
})

# Quick sanity check
nrow(completed_datasets[[1]])
str(completed_datasets[[1]][, c("time_noncvd_cr", "status_noncvd_cr")])

# Filter to valid times > 0
filtered_datasets_noncvd_rehosp <- lapply(completed_datasets, function(df) {
  df[df$time_noncvd_cr > 0 & !is.na(df$status_noncvd_cr), ]
})



# Check still have rows and the new variables are present
##should see a positive row count and both columns as numeric/integer.
nrow(filtered_datasets_noncvd_rehosp[[1]])
str(filtered_datasets_noncvd_rehosp[[1]][, c("time_noncvd_cr", "status_noncvd_cr")])
##NOTES, n=7790 (rounded to nearest 5)  patients filtered out which exactly aligns with previous cox-model analyses for non-cvd rehosps
##Filtering confirmed correct


# Fit Fine-Gray model in each imputed dataset (same covariates as HF via hf_covars)
fg_models_noncvd <- lapply(filtered_datasets_noncvd_rehosp, function(data) {
  crr(
    ftime   = data$time_noncvd_cr,
    fstatus = data$status_noncvd_cr,   # 0=censor, 1=non-CVD rehosp, 2=death
    cov1    = hf_covars(data),   ##Fully adjusted model with all covariates 
    failcode = 1,
    cencode  = 0
  )
})







lapply(fg_models_noncvd, summary) ##extract the summary from each model

#Extract term, HR, CI, p for each model
extract_fg <- function(fg_obj, imp_number) {
  
  s <- summary(fg_obj)
  coef_df <- as.data.frame(s$coef)
  
  # Add term names
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  # Compute HR and CI
  coef_df$HR          <- exp(coef_df$coef)
  coef_df$lower_95_CI <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
  coef_df$upper_95_CI <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
  
  # Keep requested columns + imputation number
  coef_df <- coef_df[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p-value")]
  coef_df$imputation <- imp_number
  
  return(coef_df)
}



#Apply to all 5 imputed models
non_cvd_fg_results_fully_adjusted <- do.call(
  rbind,
  lapply(seq_along(fg_models_noncvd), function(i) extract_fg(fg_models_noncvd[[i]], i))
)

#Save to CSV
write.csv(
  non_cvd_fg_results_fully_adjusted,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/non_cvd_rehosp/non_cvd_finegray_all_imputations_fully_adjusted.csv",
  row.names = FALSE
)




# Pool Fine-Gray models
fg_pooled_noncvd <- pool_crr(fg_models_noncvd)

# Keep term, HR, CI, p-value
fg_noncvd_df <- fg_pooled_noncvd[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
fg_noncvd_df[, 2:5] <- round(fg_noncvd_df[, 2:5], 3)

# Save to CSV
write.csv(
  fg_noncvd_df,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/non_cvd_rehosp/non_cvd_rehosp_finegray_pooled_fully_adjusted.csv",
  row.names = FALSE
)

## OPTIONAL: CIF + Gray's test for non-CVD rehosp (first imputed dataset)
d1_noncvd <- filtered_datasets_noncvd_rehosp[[1]]

ci_noncvd <- cuminc(
  ftime   = d1_noncvd$time_noncvd_cr,
  fstatus = d1_noncvd$status_noncvd_cr,   # 0=censor, 1=non-CVD rehosp, 2=death
  group   = d1_noncvd$rehab_indicator
)

ci_noncvd     # includes Gray's test
# plot(ci_noncvd)  # uncomment to see CIF curves










##################################################################
################      UNADJUSTED MODEL NON-CVD REHOSP *************************************RESTART THE DATA from TOP of sheet****************************************
###############################
# NEW: covariate matrix builder for  Fine-Gray
###############################
hf_covars_unadjusted <- function(data) {
  model.matrix(
    ~ rehab_indicator,
    data = data
  )[ , -1, drop = FALSE]  # remove the intercept column
}



################################################################################
## Non-CVD rehospitalisation - add competing-risk time & status + Fine-Gray
################################################################################

library(dplyr)
library(survival)
library(mice)
library(cmprsk)

# Load imputation object and recreate completed datasets
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")
completed_datasets <- lapply(1:5, function(i) complete(imp, i))

study_end <- as.Date("2024-03-31")

# Add original date variables + competing-risk time/status into EACH completed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  
  df$tavi_cips_disdate_plus180days <- clean_cohort$tavi_cips_disdate_plus180days[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_non_cvd_date <- clean_cohort$out_readm_non_cvd_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$date_of_death <- clean_cohort$date_of_death[
    match(df$person_id, clean_cohort$person_id)
  ]
  
  df <- df %>%
    mutate(
      t0        = tavi_cips_disdate_plus180days,
      t_noncvd  = as.numeric(out_readm_non_cvd_date - t0),
      t_death   = as.numeric(date_of_death         - t0),
      t_end     = as.numeric(study_end             - t0),
      # time to first of non-CVD rehosp, death, or end of follow-up
      time_noncvd_cr = pmin(t_noncvd, t_death, t_end, na.rm = TRUE),
      # 0 = censored, 1 = non-CVD rehosp, 2 = death before rehosp
      status_noncvd_cr = dplyr::case_when(
        !is.na(t_noncvd) & t_noncvd <= t_death & t_noncvd <= t_end ~ 1L,
        !is.na(t_death)  & (is.na(t_noncvd) | t_death < t_noncvd) & t_death <= t_end ~ 2L,
        TRUE ~ 0L
      )
    )
  
  return(df)
})

# Quick sanity check
nrow(completed_datasets[[1]])
str(completed_datasets[[1]][, c("time_noncvd_cr", "status_noncvd_cr")])

# Filter to valid times > 0
filtered_datasets_noncvd_rehosp <- lapply(completed_datasets, function(df) {
  df[df$time_noncvd_cr > 0 & !is.na(df$status_noncvd_cr), ]
})



# Check still have rows and the new variables are present
##should see a positive row count and both columns as numeric/integer.
nrow(filtered_datasets_noncvd_rehosp[[1]])
str(filtered_datasets_noncvd_rehosp[[1]][, c("time_noncvd_cr", "status_noncvd_cr")])
##NOTES, n=7790 (rounded to nearest 5)  patients filtered out which exactly aligns with previous cox-model analyses for non-cvd rehosps
##Filtering confirmed correct


# Fit Fine-Gray model in each imputed dataset (same covariates as HF via hf_covars)
fg_models_noncvd <- lapply(filtered_datasets_noncvd_rehosp, function(data) {
  crr(
    ftime   = data$time_noncvd_cr,
    fstatus = data$status_noncvd_cr,   # 0=censor, 1=non-CVD rehosp, 2=death
    cov1    = hf_covars_unadjusted(data),   ##unadjusted model
    failcode = 1,
    cencode  = 0
  )
})





lapply(fg_models_noncvd, summary) ##extract the summary from each model

#Extract term, HR, CI, p for each model
extract_fg <- function(fg_obj, imp_number) {
  
  s <- summary(fg_obj)
  coef_df <- as.data.frame(s$coef)
  
  # Add term names
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  # Compute HR and CI
  coef_df$HR          <- exp(coef_df$coef)
  coef_df$lower_95_CI <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
  coef_df$upper_95_CI <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
  
  # Keep requested columns + imputation number
  coef_df <- coef_df[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p-value")]
  coef_df$imputation <- imp_number
  
  return(coef_df)
}



#Apply to all 5 imputed models
non_cvd_fg_results_unadjusted <- do.call(
  rbind,
  lapply(seq_along(fg_models_noncvd), function(i) extract_fg(fg_models_noncvd[[i]], i))
)

#Save to CSV
write.csv(
  non_cvd_fg_results_unadjusted,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/non_cvd_rehosp/non_cvd_finegray_all_imputations_unadjusted.csv",
  row.names = FALSE
)




# Pool Fine-Gray models
fg_pooled_noncvd <- pool_crr(fg_models_noncvd)

# Keep term, HR, CI, p-value
fg_noncvd_df <- fg_pooled_noncvd[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
fg_noncvd_df[, 2:5] <- round(fg_noncvd_df[, 2:5], 3)

# Save to CSV
write.csv(
  fg_noncvd_df,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/non_cvd_rehosp/non_cvd_rehosp_finegray_pooled_unadjusted.csv",
  row.names = FALSE
)





## OPTIONAL: CIF + Gray's test for non-CVD rehosp (first imputed dataset)
d1_noncvd <- filtered_datasets_noncvd_rehosp[[1]]

ci_noncvd <- cuminc(
  ftime   = d1_noncvd$time_noncvd_cr,
  fstatus = d1_noncvd$status_noncvd_cr,   # 0=censor, 1=non-CVD rehosp, 2=death
  group   = d1_noncvd$rehab_indicator
)

ci_noncvd     # includes Gray's test
# plot(ci_noncvd)  # uncomment to see CIF curves












##DEMOGRAPHICS ADJUSTED MODEL - NON-CVD REHOSP *************************************RESTART THE DATA from TOP of sheet****************************************
###############################
# NEW: covariate matrix builder for HF Fine-Gray
###############################
hf_covars <- function(data) {
  model.matrix(
    ~ rehab_indicator + age_groups + sex_code + IMD_2019_QUINTILES + SMOKING_STATUS + 
      ethnicity_5_group + region_name,
    data = data
  )[ , -1, drop = FALSE]  # remove the intercept column
}





################################################################################
## Non-CVD rehospitalisation - add competing-risk time & status + Fine-Gray
################################################################################

library(dplyr)
library(survival)
library(mice)
library(cmprsk)

# Load imputation object and recreate completed datasets
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")
completed_datasets <- lapply(1:5, function(i) complete(imp, i))

study_end <- as.Date("2024-03-31")

# Add original date variables + competing-risk time/status into EACH completed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  
  df$tavi_cips_disdate_plus180days <- clean_cohort$tavi_cips_disdate_plus180days[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$out_readm_non_cvd_date <- clean_cohort$out_readm_non_cvd_date[
    match(df$person_id, clean_cohort$person_id)
  ]
  df$date_of_death <- clean_cohort$date_of_death[
    match(df$person_id, clean_cohort$person_id)
  ]
  
  df <- df %>%
    mutate(
      t0        = tavi_cips_disdate_plus180days,
      t_noncvd  = as.numeric(out_readm_non_cvd_date - t0),
      t_death   = as.numeric(date_of_death         - t0),
      t_end     = as.numeric(study_end             - t0),
      # time to first of non-CVD rehosp, death, or end of follow-up
      time_noncvd_cr = pmin(t_noncvd, t_death, t_end, na.rm = TRUE),
      # 0 = censored, 1 = non-CVD rehosp, 2 = death before rehosp
      status_noncvd_cr = dplyr::case_when(
        !is.na(t_noncvd) & t_noncvd <= t_death & t_noncvd <= t_end ~ 1L,
        !is.na(t_death)  & (is.na(t_noncvd) | t_death < t_noncvd) & t_death <= t_end ~ 2L,
        TRUE ~ 0L
      )
    )
  
  return(df)
})

# Quick sanity check
nrow(completed_datasets[[1]])
str(completed_datasets[[1]][, c("time_noncvd_cr", "status_noncvd_cr")])

# Filter to valid times > 0
filtered_datasets_noncvd_rehosp <- lapply(completed_datasets, function(df) {
  df[df$time_noncvd_cr > 0 & !is.na(df$status_noncvd_cr), ]
})



# Check still have rows and the new variables are present
##should see a positive row count and both columns as numeric/integer.
nrow(filtered_datasets_noncvd_rehosp[[1]])
str(filtered_datasets_noncvd_rehosp[[1]][, c("time_noncvd_cr", "status_noncvd_cr")])
##NOTES, n=7790 (rounded to nearest 5)  patients filtered out which exactly aligns with previous cox-model analyses for non-cvd rehosps
##Filtering confirmed correct


# Fit Fine-Gray model in each imputed dataset (same covariates as HF via hf_covars)
fg_models_noncvd <- lapply(filtered_datasets_noncvd_rehosp, function(data) {
  crr(
    ftime   = data$time_noncvd_cr,
    fstatus = data$status_noncvd_cr,   # 0=censor, 1=non-CVD rehosp, 2=death
    cov1    = hf_covars(data),   ##demographics adjusted model with demographic covariates 
    failcode = 1,
    cencode  = 0
  )
})






lapply(fg_models_noncvd, summary) ##extract the summary from each model

#Extract term, HR, CI, p for each model
extract_fg <- function(fg_obj, imp_number) {
  
  s <- summary(fg_obj)
  coef_df <- as.data.frame(s$coef)
  
  # Add term names
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  # Compute HR and CI
  coef_df$HR          <- exp(coef_df$coef)
  coef_df$lower_95_CI <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
  coef_df$upper_95_CI <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
  
  # Keep requested columns + imputation number
  coef_df <- coef_df[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p-value")]
  coef_df$imputation <- imp_number
  
  return(coef_df)
}



#Apply to all 5 imputed models
non_cvd_fg_results_demographics_adjusted <- do.call(
  rbind,
  lapply(seq_along(fg_models_noncvd), function(i) extract_fg(fg_models_noncvd[[i]], i))
)

#Save to CSV
write.csv(
  non_cvd_fg_results_demographics_adjusted,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/non_cvd_rehosp/non_cvd_finegray_all_imputations_demographics_adjusted.csv",
  row.names = FALSE
)





# Pool Fine-Gray models
fg_pooled_noncvd <- pool_crr(fg_models_noncvd)

# Keep term, HR, CI, p-value
fg_noncvd_df <- fg_pooled_noncvd[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
fg_noncvd_df[, 2:5] <- round(fg_noncvd_df[, 2:5], 3)

# Save to CSV
write.csv(
  fg_noncvd_df,
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/competing risk of death models/non_cvd_rehosp/non_cvd_rehosp_finegray_pooled_demographic_adjusted.csv",
  row.names = FALSE
)




## OPTIONAL: CIF + Gray's test for non-CVD rehosp (first imputed dataset)
d1_noncvd <- filtered_datasets_noncvd_rehosp[[1]]

ci_noncvd <- cuminc(
  ftime   = d1_noncvd$time_noncvd_cr,
  fstatus = d1_noncvd$status_noncvd_cr,   # 0=censor, 1=non-CVD rehosp, 2=death
  group   = d1_noncvd$rehab_indicator
)

ci_noncvd     # includes Gray's test
# plot(ci_noncvd)  # uncomment to see CIF curves



######DONE:
#HF rehosp -pooled fully adjusted, pooled demographics adjusted and unadjusted 
#all cause rehosp - pooled fully adjusted, pooled demographics adjusted and unadjusted 
#non-CVD rehosps - pooled fully adjusted, pooled demographics adjusted and unadjusted 










