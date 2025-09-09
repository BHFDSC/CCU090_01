#####Multiple Imputation by chained equations for dealing with missing data ##########

install.packages("mice")      # if not already installed
install.packages("survival")  # if you're using Cox models
library(mice)
library(survival)


library(dplyr)
library(lubridate)

###convert all "missing" values to actual NA values across relevant columns###

# Define the variables as a character vector
vars_to_check <- c("age_at_tavi", "sex_code", "ethnicity_5_group", "IMD_2019_QUINTILES", "region_name", "SMOKING_STATUS",
                   "BMI_category", "NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE",
                   "binary_known_poor_mobility", "binary_Known_DIABETES",
                   "KATZ_INDEX_CAT", "COMORBIDITY_pulmonary_liver_neuro_excardiacvascu", "LV_FUNCTION",
                   "previous_MI_time_grouped", "binary_known_prev_cardiac_surgery", 
                   "MITRAL_REGURGITATION", "TAVI_procedural_complications", "binary_known_on_dialysis", "BINARY_known_SEVERE_LIVER_DISEASE", "binary_known_pulmonary_disease",
                   "BINARY_known_HISTORY_OF_NEUROLOGICAL_DISEASE", "binary_known_extracardiac_arteriopathy",
                   "binary_known_prev_mi", "binary_known_prev_pci","binary_previous_known_CABG", "binary_previous_known_valve_operation",
                   "binary_previous_known_other_operation_pericardium",
                   "BINARY_known_PERMANENT_PACING", "BINARY_known_BLEEDING", "BINARY_known_USE_OF_CARDIOPULMONARY_BYPASS")

 



# Replace all instances of "missing" with NA across these variables
clean_cohort[vars_to_check] <- lapply(clean_cohort[vars_to_check], function(x) {
  x <- as.character(x)
  x[x %in% c("Missing", "missing")] <- NA  # cover both cases
  return(as.factor(x))
})



# Create survival object for time to cardiac rehab exposure while accounting for censoring due to death or the end of follow-up
clean_cohort <- clean_cohort %>%
  mutate(
    rehab_status = ifelse(is.na(exp_cardiac_rehab_6m), 0, ifelse(exp_cardiac_rehab_6m == 1, 1, 0)),  # 1 for exposed, 0 for censored
    cardiac_rehab_time = cardiac_rehab_follow_up_days        # Time in days (time to rehab or censoring)
  )



######Create survival object for time to HF rehosp while accounting for censoring due to death or the end of follow-up
##Explanation:
# Define study start and end dates
study_start_new <- clean_cohort$tavi_cips_disdate_plus180days
study_end <- as.Date("2024-03-31")

# Add necessary variables
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




#####Create survival object for time to ALL_CAUSE rehosp while accounting for censoring due to death or the end of follow-up

##Explanation:
# Define study start and end dates
study_start_new <- clean_cohort$tavi_cips_disdate_plus180days
study_end <- as.Date("2024-03-31")

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





#####Create survival object for time to NON-CARDIOVASCULAR rehosp while accounting for censoring due to death or the end of follow-up

##Explanation:
# Define study start and end dates
study_start_new <- clean_cohort$tavi_cips_disdate_plus180days
study_end <- as.Date("2024-03-31")

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





#####Create survival object for time to MORTALITY while accounting for censoring due to death or the end of follow-up

##Explanation:
# Define study start and end dates
study_start_new <- clean_cohort$tavi_cips_disdate_plus180days
study_end <- as.Date("2024-03-31")

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







#####ADD Nelson-Aalen Estimators for Each Outcome####

library(survival)

# ==== 1. Cardiac Rehab ====
clean_cohort$na_cumhaz_rehab <- nelsonaalen(clean_cohort,"cardiac_rehab_time","rehab_status")


# ==== 2. Heart Failure Rehospitalisation ====
clean_cohort$na_cumhaz_hf <- nelsonaalen(clean_cohort,"hf_rehosp_censor_time","hf_rehosp_indicator")

# ==== 3. All-Cause Rehospitalisation ====
clean_cohort$na_cumhaz_allcause <- nelsonaalen(clean_cohort,"all_cause_rehosp_censor_time","all_cause_rehosp_indicator")


# ==== 4. Non-CVD Rehospitalisation ====
clean_cohort$na_cumhaz_noncvd <- nelsonaalen(clean_cohort,"non_cvd_rehosp_censor_time","non_cvd_rehosp_indicator")


# ==== 5. Mortality ====
clean_cohort$na_cumhaz_mortality <- nelsonaalen(clean_cohort,"all_cause_mortality_censor_time","all_cause_mortality_indicator")







######  Prepare for Multiple Imputation with Predictor Matrix  ########

 


library(mice)



# Define the variables to impute
vars_to_impute <- c("ethnicity_5_group", "SMOKING_STATUS",
                    "BMI_category", "NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE", "KATZ_INDEX_CAT", "LV_FUNCTION",
                    "previous_MI_time_grouped", "MITRAL_REGURGITATION")


##Incorporate Nelson-Aalen Estimators and event indicators into the MICE setup 
surv_vars <- c(
  "na_cumhaz_rehab", "na_cumhaz_hf", "na_cumhaz_allcause",
  "na_cumhaz_noncvd", "na_cumhaz_mortality", 
  "rehab_status", "hf_rehosp_indicator", "all_cause_rehosp_indicator", "non_cvd_rehosp_indicator", "all_cause_mortality_indicator"
)


# Define variables to exclude from imputation and prediction
vars_exclude <- c(
  "rehab_indicator", "age_groups", "year_tavi_discharge",
  "exp_cardiac_rehab_attends_6m", "exp_cardiac_rehab_6m", "exp_cardiac_rehab_3m", "rehab_attendance_groups",
  "cardiac_rehab_time", "hf_rehosp_censor_time", "all_cause_rehosp_censor_time","non_cvd_rehosp_censor_time", "all_cause_mortality_censor_time",
  "tavi_cips_disdate", "tavi_cips_disdate_plus180days", "tavi_cips_disdate_plus90_days" 
)

# 1. Create a working copy of dataset
clean_cohort_mice <- clean_cohort [,c("person_id", vars_to_check, surv_vars, vars_exclude)]
                                      

# 2. Let mice guess default method and predictor matrix
meth <- make.method(clean_cohort_mice)
pred <- make.predictorMatrix(clean_cohort_mice)

#Specify criteria for age
meth["age_at_tavi"] <- ""         # do not impute
pred["age_at_tavi", ] <- 1        # predict other variables
pred[, "age_at_tavi"] <- 0        # Do not allow it to be predicted

# Exclude ID
pred[, "person_id"] <- 0
pred["person_id", ] <- 0
meth["person_id"] <- ""

# Exclude survival/event variables from being imputed, but allow them as predictors
meth[surv_vars] <- ""  #disables imputation for these variables
pred[surv_vars, ] <- 1  #Tells MICE that these variables can be used to predict other variables
pred[, surv_vars] <- 0  #Tells MICE that these variables should not be predicted by others.


# Exclude variables to exclude from being imputed and Prevent them from predicting others
meth[vars_exclude] <- ""   # Exclude variables to exclude from being imputed
pred[vars_exclude, ] <- 0  # don't predict others
pred[ , vars_exclude] <- 0  # don't be used as predictors


# 3. Run MICE
imp <- mice(
  data = clean_cohort_mice,
  m = 5, #should be 5
  maxit=30, #should be 30
  method = meth,
  predictorMatrix = pred,
  seed = 123
)


# 4. View summary
summary(imp)

#check varibles it excludes from imputation
imp$loggedEvents

#check the traces
plot(imp)

##traces look good

# Show imputed values for a specific variable
imp$imp$SMOKING_STATUS
imp$imp$ethnicity_5_group
imp$imp$BMI_category
imp$imp$NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE
imp$imp$KATZ_INDEX_CAT




# 5. Create a list of all 5 completed datasets
completed_datasets <- lapply(1:5, function(i) complete(imp, i))




###SAVING ####
# Define  save path
save_path <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data"


#  1. Save the full mids object

saveRDS(imp, file = file.path(save_path, "imp_mice_object.rds"))


# 2. Save the 5 completed datasets as .RDS files

# Save each imputed dataset as an RDS
for (i in 1:5) {
  completed <- complete(imp, i)
  saveRDS(completed, file = file.path(save_path, paste0("completed_dataset_", i, ".rds")))
}



# 3. Save as .csv

# Save as CSVs too (optional)
for (i in 1:5) {
  completed <- complete(imp, i)
  write.csv(completed, file = file.path(save_path, paste0("completed_dataset_", i, ".csv")), row.names = FALSE)
}



# 4. save Long format with all imputations stacked
# Convert to long format
long_data <- complete(imp, action = "long", include = TRUE)
#Save long format data
# Save as RDS
saveRDS(long_data, file = file.path(save_path, "long_format_imputed_data.rds"))

# Save as CSV
write.csv(long_data, file = file.path(save_path, "long_format_imputed_data.csv"), row.names = FALSE)




#########################################################################################################################################################################################################
#                                                                               LOAD DATA and fit cox models
#########################################################################################################################################################################################################

### fit Cox models ### 


library(survival)
library(mice)

#Load full mice object and recreate completed datasets
# Load full MICE mids object
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")

# Recreate list of completed datasets
completed_datasets <- lapply(1:5, function(i) complete(imp, i))



# Add days_since_tavi_dis_plus180 back to each imputed dataset ## Ensure it has been added to clean_cohort (see code at the start of this script)
# Add the variable back to each imputed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[match(df$person_id, clean_cohort$person_id)]
  return(df)
})

# Confirm it's added correctly in the first dataset
head(completed_datasets[[1]]$days_since_tavi_dis_plus180)
summary(completed_datasets[[1]]$days_since_tavi_dis_plus180)


###################################################################################################################################################################################

############################################################ ==== 1. Heart Failure Rehospitalisation ==== #########################################################################

### Need to filter events within 180 days for each of the outcomes (do each outcome separately)

# Pre-filter all datasets first
filtered_datasets_hf_rehosp <- lapply(completed_datasets, function(df) {
  df[df$hf_rehosp_censor_time > df$days_since_tavi_dis_plus180, ]
})

                 # ==== Adjust for demographics # ==== adjusting for demographics: age, sex, region, ethnicity, SES, smoking status

#  Fit Cox models on each imputed dataset with filter
#Fit the Cox proportional hazards model
# For outcome: hf_rehosp_indicator
# Exposure: rehab_indicator
# Time: hf_rehosp_censor_time

# ==== hf rehosp Adjust for demographics # ====
cox_models_hf_rehosp_demographics_adjusted <- lapply(filtered_datasets_hf_rehosp, function(data) {
  coxph(Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator + age_groups + sex_code + ethnicity_5_group + IMD_2019_QUINTILES +strata(region_name) + SMOKING_STATUS,
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_hf_rehosp_demographics_adjusted <- pool(cox_models_hf_rehosp_demographics_adjusted)

# View results
summary(pooled_hf_rehosp_demographics_adjusted)


# Extract the summary of the pooled model
summary_df_hf_rehosp_adjusted_demographics <- summary(pooled_hf_rehosp_demographics_adjusted)

# Add columns for HR and 95% CI
summary_df_hf_rehosp_adjusted_demographics$HR <- exp(summary_df_hf_rehosp_adjusted_demographics$estimate)
summary_df_hf_rehosp_adjusted_demographics$lower_95_CI <- exp(summary_df_hf_rehosp_adjusted_demographics$estimate - 1.96 * summary_df_hf_rehosp_adjusted_demographics$std.error)
summary_df_hf_rehosp_adjusted_demographics$upper_95_CI <- exp(summary_df_hf_rehosp_adjusted_demographics$estimate + 1.96 * summary_df_hf_rehosp_adjusted_demographics$std.error)

# Round results for clarity (optional)
summary_df_hf_rehosp_adjusted_demographics <- summary_df_hf_rehosp_adjusted_demographics[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
summary_df_hf_rehosp_adjusted_demographics <- round(summary_df_hf_rehosp_adjusted_demographics, 3)

# Print the final table
print(summary_df_hf_rehosp_adjusted_demographics)


## SAVE OUTPUT ###

# Define output path
output_path_hf_rehosp <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/hf_rehosp_from_imputed_data"

# Define filename
output_file_hf_rehosp <- "output_df_hf_rehosp_adjusted_demographics.csv"

# Combine full path
full_output_path_hf_rehosp <- file.path(output_path_hf_rehosp, output_file_hf_rehosp)

# Save as CSV
write.csv(summary_df_hf_rehosp_adjusted_demographics, 
          file = full_output_path_hf_rehosp, 
          row.names = FALSE)





            # ====hf rehosp Adjusted for all covarites # ====
#Fit Cox model with calendar time as the timescale, with all covariates 
### stratify binary_known_poor_mobility and region_name

cox_models_hf_rehosp_fully_adjusted <- lapply(filtered_datasets_hf_rehosp, function(data) {
  coxph(Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator + age_groups + sex_code + IMD_2019_QUINTILES + SMOKING_STATUS + 
          BMI_category + KATZ_INDEX_CAT +  COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
          previous_MI_time_grouped + binary_known_prev_cardiac_surgery + year_tavi_discharge +  # Include year
          strata(binary_known_poor_mobility, region_name) +           # stratify hazard assumption violating variable 
          binary_Known_DIABETES + 
          TAVI_procedural_complications + ethnicity_5_group + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + LV_FUNCTION + MITRAL_REGURGITATION,
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_hf_rehosp_fully_adjusted <- pool(cox_models_hf_rehosp_fully_adjusted)

# View results
summary(pooled_hf_rehosp_fully_adjusted)




# === Extract summary and calculate HRs and 95% CIs ===
summary_df_hf_rehosp_adjusted_full <- summary(pooled_hf_rehosp_fully_adjusted)

# Add exponentiated HR and confidence intervals
summary_df_hf_rehosp_adjusted_full$HR <- exp(summary_df_hf_rehosp_adjusted_full$estimate)
summary_df_hf_rehosp_adjusted_full$lower_95_CI <- exp(summary_df_hf_rehosp_adjusted_full$estimate - 1.96 * summary_df_hf_rehosp_adjusted_full$std.error)
summary_df_hf_rehosp_adjusted_full$upper_95_CI <- exp(summary_df_hf_rehosp_adjusted_full$estimate + 1.96 * summary_df_hf_rehosp_adjusted_full$std.error)

# Select and round output 
summary_df_hf_rehosp_adjusted_full <- summary_df_hf_rehosp_adjusted_full[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
summary_df_hf_rehosp_adjusted_full <- round(summary_df_hf_rehosp_adjusted_full, 3)

# === Save to CSV ===
write.csv(
  summary_df_hf_rehosp_adjusted_full,
  file = "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/hf_rehosp_from_imputed_data/output_df_cox_model_hf_rehosp_adjusted_fully.csv",
  row.names = FALSE
)





############################################################################################################################################################################

############################################################==== 2. All-Cause Rehospitalisation ==== #################################################################

#re-start the data - need new filter

library(survival)
library(mice)

#Load full mice object and recreate completed datasets
# Load full MICE mids object
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")

# Recreate list of completed datasets
completed_datasets <- lapply(1:5, function(i) complete(imp, i))



# Add days_since_tavi_dis_plus180 back to each imputed dataset # Ensure it has been added to clean_cohort (see code at the start of this script)
# Add the variable back to each imputed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[match(df$person_id, clean_cohort$person_id)]
  return(df)
})

# Confirm it's added correctly in the first dataset
head(completed_datasets[[1]]$days_since_tavi_dis_plus180)
summary(completed_datasets[[1]]$days_since_tavi_dis_plus180)



### Need to filter events within 180 days for each of the outcomes (do each outcome separately)

# Pre-filter all datasets first
filtered_datasets_all_cause_rehosp <- lapply(completed_datasets, function(df) {
  df[df$all_cause_rehosp_censor_time > df$days_since_tavi_dis_plus180, ]
})

#Check filter is correct
# Compare number of rows before and after filtering
for (i in 1:length(completed_datasets)) {
  before_n <- nrow(completed_datasets[[i]])
  after_n <- nrow(filtered_datasets_all_cause_rehosp[[i]])
  cat("Dataset", i, ": Before =", before_n, "| After =", after_n, "| Filtered =", before_n - after_n, "\n")
}

#Filtering is confirmed correct 


# ==== Adjust for demographics # ==== adjusting for demographics: age, sex, region, ethnicity, SES, smoking status

#  Fit Cox models on each imputed dataset with filter
#Fit the Cox proportional hazards model
# For outcome: all_cause_rehosp_indicator
# Exposure: rehab_indicator
# Time: all_cause_rehosp_censor_time

# ==== all_cause rehosp Adjust for demographics # ====
cox_models_all_cause_rehosp_demographics_adjusted <- lapply(filtered_datasets_all_cause_rehosp, function(data) {
  coxph(Surv(all_cause_rehosp_censor_time, all_cause_rehosp_indicator) ~ rehab_indicator + strata(age_groups) +   #Stratify by age groups to account for proportional hazards 
        sex_code + ethnicity_5_group + IMD_2019_QUINTILES + region_name + SMOKING_STATUS,
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_all_cause_rehosp_demographics_adjusted <- pool(cox_models_all_cause_rehosp_demographics_adjusted)

# View results
summary(pooled_all_cause_rehosp_demographics_adjusted)



# === Extract summary, calculate HR and 95% CI, and save ===

# Extract the summary
summary_df_all_cause_rehosp_adjusted_demographics <- summary(pooled_all_cause_rehosp_demographics_adjusted)

# Add columns for HR and 95% CI
summary_df_all_cause_rehosp_adjusted_demographics$HR <- exp(summary_df_all_cause_rehosp_adjusted_demographics$estimate)
summary_df_all_cause_rehosp_adjusted_demographics$lower_95_CI <- exp(summary_df_all_cause_rehosp_adjusted_demographics$estimate - 1.96 * summary_df_all_cause_rehosp_adjusted_demographics$std.error)
summary_df_all_cause_rehosp_adjusted_demographics$upper_95_CI <- exp(summary_df_all_cause_rehosp_adjusted_demographics$estimate + 1.96 * summary_df_all_cause_rehosp_adjusted_demographics$std.error)

# Round the results
summary_df_all_cause_rehosp_adjusted_demographics <- round(summary_df_all_cause_rehosp_adjusted_demographics, 3)
View(summary_df_all_cause_rehosp_adjusted_demographics)

# Save to CSV file
write.csv(
  summary_df_all_cause_rehosp_adjusted_demographics[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")],
  file = "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/all_cause_rehosp_from_imputed_data/output_df_all_cause_rehosp_adjusted_demographics.csv",
  row.names = FALSE
)






# ====    all cause rehosp Adjusted for all covarites # ====
#Fit Cox model with calendar time as the timescale, with all covariates 
### ### stratify age_groups and COMORBIDITY

cox_models_all_cause_rehosp_fully_adjusted <- lapply(filtered_datasets_all_cause_rehosp, function(data) {
  coxph(Surv(all_cause_rehosp_censor_time, all_cause_rehosp_indicator) ~ rehab_indicator +
          sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + region_name +
          BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + 
          strata(COMORBIDITY_pulmonary_liver_neuro_excardiacvascu, age_groups) + #stratify hazard violating variables age groups and comorbidity
          previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications + year_tavi_discharge,  # Include year
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_all_cause_rehosp_fully_adjusted <- pool(cox_models_all_cause_rehosp_fully_adjusted)

# View results
summary(pooled_all_cause_rehosp_fully_adjusted)



# === Extract summary, calculate HR and 95% CI, and save ===

# Extract the summary
summary_df_all_cause_rehosp_fully_adjusted <- summary(pooled_all_cause_rehosp_fully_adjusted)

# Add columns for HR and 95% CI
summary_df_all_cause_rehosp_fully_adjusted$HR <- exp(summary_df_all_cause_rehosp_fully_adjusted$estimate)
summary_df_all_cause_rehosp_fully_adjusted$lower_95_CI <- exp(summary_df_all_cause_rehosp_fully_adjusted$estimate - 1.96 * summary_df_all_cause_rehosp_fully_adjusted$std.error)
summary_df_all_cause_rehosp_fully_adjusted$upper_95_CI <- exp(summary_df_all_cause_rehosp_fully_adjusted$estimate + 1.96 * summary_df_all_cause_rehosp_fully_adjusted$std.error)

# Round the results
summary_df_all_cause_rehosp_fully_adjusted <- round(summary_df_all_cause_rehosp_fully_adjusted, 3)

# Save to CSV file
write.csv(
  summary_df_all_cause_rehosp_fully_adjusted[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")],
  file = "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/all_cause_rehosp_from_imputed_data/output_df_all_cause_rehosp_fully_adjusted.csv",
  row.names = FALSE
)







############################################################################################################################################################################

##########################################         ==== 3. Non-CVD Rehospitalisation ====     #####################################################################################

#restart the data - need new filter

library(survival)
library(mice)

#Load full mice object and recreate completed datasets
# Load full MICE mids object
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")

# Recreate list of completed datasets
completed_datasets <- lapply(1:5, function(i) complete(imp, i))



# Add days_since_tavi_dis_plus180 back to each imputed dataset # Ensure it has been added to clean_cohort (see code at the start of this script)
# Add the variable back to each imputed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[match(df$person_id, clean_cohort$person_id)]
  return(df)
})

# Confirm it's added correctly in the first dataset
head(completed_datasets[[1]]$days_since_tavi_dis_plus180)
summary(completed_datasets[[1]]$days_since_tavi_dis_plus180)



### Need to filter events within 180 days for each of the outcomes (do each outcome separately)

# Pre-filter all datasets first
filtered_datasets_non_cvd_rehosp <- lapply(completed_datasets, function(df) {
  df[df$non_cvd_rehosp_censor_time > df$days_since_tavi_dis_plus180, ]
})

#Check filter is correct
# Compare number of rows before and after filtering
for (i in 1:length(completed_datasets)) {
  before_n <- nrow(completed_datasets[[i]])
  after_n <- nrow(filtered_datasets_non_cvd_rehosp[[i]])
  cat("Dataset", i, ": Before =", before_n, "| After =", after_n, "| Filtered =", before_n - after_n, "\n")
}

#Filtering is confirmed correct 





# ==== Adjust for demographics # ==== adjusting for demographics: age, sex, region, ethnicity, SES, smoking status

#  Fit Cox models on each imputed dataset with filter
#Fit the Cox proportional hazards model
# For outcome: non_cvd_rehosp_indicator
# Exposure: rehab_indicator
# Time: non_cvd_rehosp_censor_time

###startify by age groups b/c proportional hazards doesnt hold true for age

# ==== Non-CVD rehosp Adjust for demographics # ====
cox_models_non_cvd_rehosp_demographics_adjusted <- lapply(filtered_datasets_non_cvd_rehosp, function(data) {
  coxph(Surv(non_cvd_rehosp_censor_time, non_cvd_rehosp_indicator) ~ rehab_indicator + strata(age_groups) + sex_code + ethnicity_5_group + IMD_2019_QUINTILES +region_name + SMOKING_STATUS,
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_non_cvd_rehosp_demographics_adjusted <- pool(cox_models_non_cvd_rehosp_demographics_adjusted)

# View results
summary(pooled_non_cvd_rehosp_demographics_adjusted)



# === Extract summary of pooled model ===
summary_df_non_cvd_rehosp_demographics_adjusted <- summary(pooled_non_cvd_rehosp_demographics_adjusted)

# === Add HR and 95% CI columns ===
summary_df_non_cvd_rehosp_demographics_adjusted$HR <- exp(summary_df_non_cvd_rehosp_demographics_adjusted$estimate)
summary_df_non_cvd_rehosp_demographics_adjusted$lower_95_CI <- exp(summary_df_non_cvd_rehosp_demographics_adjusted$estimate - 1.96 * summary_df_non_cvd_rehosp_demographics_adjusted$std.error)
summary_df_non_cvd_rehosp_demographics_adjusted$upper_95_CI <- exp(summary_df_non_cvd_rehosp_demographics_adjusted$estimate + 1.96 * summary_df_non_cvd_rehosp_demographics_adjusted$std.error)

# === Select and round final output ===
summary_df_non_cvd_rehosp_demographics_adjusted <- summary_df_non_cvd_rehosp_demographics_adjusted[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
summary_df_non_cvd_rehosp_demographics_adjusted <- round(summary_df_non_cvd_rehosp_demographics_adjusted, 3)

# === Save as CSV ===
write.csv(
  summary_df_non_cvd_rehosp_demographics_adjusted,
  file = "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/non_cvd_rehosp_from_imputed_data/non_cvd_rehosp_demographics_adjusted.csv",
  row.names = FALSE
)





# ====    Non_CVD rehosp Adjusted for all covarites # ====
#Fit Cox model with calendar time as the timescale, with all covariates 
######Stratify by age_groups; NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE; and COMORBIDITY_pulmonary_liver_neuro_excardiacvascu - these did not meet proportional hazards assumption

cox_models_non_cvd_rehosp_fully_adjusted <- lapply(filtered_datasets_non_cvd_rehosp, function(data) {
  coxph(Surv(non_cvd_rehosp_censor_time, non_cvd_rehosp_indicator) ~ rehab_indicator + 
          sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + region_name +
          BMI_category + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES +
          previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications + year_tavi_discharge +  # Include year
          strata(age_groups, COMORBIDITY_pulmonary_liver_neuro_excardiacvascu, NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE), #stratify hazard violoting variable and #use age_groups categorical (reference group is <70years old)
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_non_cvd_rehosp_fully_adjusted <- pool(cox_models_non_cvd_rehosp_fully_adjusted)

# View results
summary(pooled_non_cvd_rehosp_fully_adjusted)




# === Extract summary of pooled model ===
summary_df_non_cvd_rehosp_fully_adjusted <- summary(pooled_non_cvd_rehosp_fully_adjusted)

# === Add HR and 95% CI columns ===
summary_df_non_cvd_rehosp_fully_adjusted$HR <- exp(summary_df_non_cvd_rehosp_fully_adjusted$estimate)
summary_df_non_cvd_rehosp_fully_adjusted$lower_95_CI <- exp(summary_df_non_cvd_rehosp_fully_adjusted$estimate - 1.96 * summary_df_non_cvd_rehosp_fully_adjusted$std.error)
summary_df_non_cvd_rehosp_fully_adjusted$upper_95_CI <- exp(summary_df_non_cvd_rehosp_fully_adjusted$estimate + 1.96 * summary_df_non_cvd_rehosp_fully_adjusted$std.error)

# === Select and round final output ===
summary_df_non_cvd_rehosp_fully_adjusted <- summary_df_non_cvd_rehosp_fully_adjusted[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
summary_df_non_cvd_rehosp_fully_adjusted <- round(summary_df_non_cvd_rehosp_fully_adjusted, 3)

# === Save as CSV ===
write.csv(
  summary_df_non_cvd_rehosp_fully_adjusted,
  file = "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/non_cvd_rehosp_from_imputed_data/non_cvd_rehosp_fully_adjusted.csv",
  row.names = FALSE
)








############################################################################################################################################################################

#####################################################   ==== 4. MORTALITY     ==== #####################################################################################

#restart the data - need new filter

library(survival)
library(mice)

#Load full mice object and recreate completed datasets
# Load full MICE mids object
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")

# Recreate list of completed datasets
completed_datasets <- lapply(1:5, function(i) complete(imp, i))



# Add days_since_tavi_dis_plus180 back to each imputed dataset # Ensure it has been added to clean_cohort (see code at the start of this script)
# Add the variable back to each imputed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[match(df$person_id, clean_cohort$person_id)]
  return(df)
})

# Confirm it's added correctly in the first dataset
head(completed_datasets[[1]]$days_since_tavi_dis_plus180)
summary(completed_datasets[[1]]$days_since_tavi_dis_plus180)





### Need to filter events within 180 days for each of the outcomes (do each outcome separately)

# Pre-filter all datasets first
filtered_datasets_mortality <- lapply(completed_datasets, function(df) {
  df[df$all_cause_mortality_censor_time > df$days_since_tavi_dis_plus180, ]
})

#Check filter is correct
# Compare number of rows before and after filtering
for (i in 1:length(completed_datasets)) {
  before_n <- nrow(completed_datasets[[i]])
  after_n <- nrow(filtered_datasets_mortality[[i]])
  cat("Dataset", i, ": Before =", before_n, "| After =", after_n, "| Filtered =", before_n - after_n, "\n")
}

#Filtering is confirmed correct 




# ==== Adjust for demographics # ==== adjusting for demographics: age, sex, region, ethnicity, SES, smoking status

#  Fit Cox models on each imputed dataset with filter
#Fit the Cox proportional hazards model
# For outcome: all_cause_mortality_indicator
# Exposure: rehab_indicator
# Time: all_cause_mortality_censor_time

###startify by age_groups, sex_code, region_name, ethnicity_5_group - these did not meet hazards assumption test

# ==== Mortality Adjust for demographics # ====
cox_models_mortality_demographics_adjusted <- lapply(filtered_datasets_mortality, function(data) {
  coxph(Surv(all_cause_mortality_censor_time, all_cause_mortality_indicator) ~ rehab_indicator + strata(age_groups, sex_code, region_name, ethnicity_5_group) + IMD_2019_QUINTILES + SMOKING_STATUS,
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_mortality_demographics_adjusted <- pool(cox_models_mortality_demographics_adjusted)

# View results
summary(pooled_mortality_demographics_adjusted)



# Extract the summary of the pooled model
summary_df_mortality_demographics_adjusted <- summary(pooled_mortality_demographics_adjusted)

# Add columns for HR and 95% CI
summary_df_mortality_demographics_adjusted$HR <- exp(summary_df_mortality_demographics_adjusted$estimate)
summary_df_mortality_demographics_adjusted$lower_95_CI <- exp(summary_df_mortality_demographics_adjusted$estimate - 1.96 * summary_df_mortality_demographics_adjusted$std.error)
summary_df_mortality_demographics_adjusted$upper_95_CI <- exp(summary_df_mortality_demographics_adjusted$estimate + 1.96 * summary_df_mortality_demographics_adjusted$std.error)

# Round results
summary_df_mortality_demographics_adjusted <- summary_df_mortality_demographics_adjusted[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
summary_df_mortality_demographics_adjusted <- round(summary_df_mortality_demographics_adjusted, 3)

# Save to CSV
write.csv(summary_df_mortality_demographics_adjusted, 
          "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/mortality_from_imputed_data/mortality_demographics_adjusted.csv",
          row.names = FALSE)








# ====    Mortality Adjusted for all covarites # ====
#Fit Cox model with calendar time as the timescale, with all covariates 

### stratify smoking as well as age_groups, sex_code, ethnicity_5_group, NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE, binary_Known_DIABETES, TAVI_procedural_complications, year_tavi_discharge - these did not meet proportional hazards assumption

cox_models_mortality_fully_adjusted <- lapply(filtered_datasets_mortality, function(data) {
  coxph(Surv(all_cause_mortality_censor_time, all_cause_mortality_indicator) ~ rehab_indicator + 
          IMD_2019_QUINTILES + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + region_name +
          BMI_category + KATZ_INDEX_CAT +  binary_known_poor_mobility +
          previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications +
          strata(age_groups, sex_code, ethnicity_5_group, NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE, binary_Known_DIABETES, SMOKING_STATUS, year_tavi_discharge), #stratify hazard violating variable and #use age_groups categorical (reference group is <70years old)
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_mortality_fully_adjusted <- pool(cox_models_mortality_fully_adjusted)

# View results
summary(pooled_mortality_fully_adjusted)




# Extract the summary of the pooled fully adjusted mortality model
summary_df_mortality_fully_adjusted <- summary(pooled_mortality_fully_adjusted)

# Add columns for HR and 95% CI
summary_df_mortality_fully_adjusted$HR <- exp(summary_df_mortality_fully_adjusted$estimate)
summary_df_mortality_fully_adjusted$lower_95_CI <- exp(summary_df_mortality_fully_adjusted$estimate - 1.96 * summary_df_mortality_fully_adjusted$std.error)
summary_df_mortality_fully_adjusted$upper_95_CI <- exp(summary_df_mortality_fully_adjusted$estimate + 1.96 * summary_df_mortality_fully_adjusted$std.error)

# Round results for clarity
summary_df_mortality_fully_adjusted <- summary_df_mortality_fully_adjusted[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
summary_df_mortality_fully_adjusted <- round(summary_df_mortality_fully_adjusted, 3)

# Save to CSV
write.csv(summary_df_mortality_fully_adjusted, 
          "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/mortality_from_imputed_data/mortality_fully_adjusted.csv",
          row.names = FALSE)



#########################################################################################################################################################################################################



############################################################################################################################################################################

##########################################         ==== 5. PREDICTORS of CARDIAC REHAB ====     #####################################################################################
#do 1 cox model predicting time to rehab 
#do 2 Multivariable logistic regression modelpredicting rehab vs no rehab binary


#restart the data - need new filter

library(survival)
library(mice)

#Load full mice object and recreate completed datasets
# Load full MICE mids object
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")

# Recreate list of completed datasets
completed_datasets <- lapply(1:5, function(i) complete(imp, i))


# Relevel 'ethnicity_5_group' to have 'White' as the reference level
completed_datasets <- lapply(completed_datasets, function(df) {
  df$ethnicity_5_group <- relevel(factor(df$ethnicity_5_group), ref = "White")
  return(df)
})


# Add days_since_tavi_dis_plus180 back to each imputed dataset # Ensure it has been added to clean_cohort (see code at the start of this script)
# Add the variable back to each imputed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[match(df$person_id, clean_cohort$person_id)]
  return(df)
})

# Confirm it's added correctly in the first dataset
head(completed_datasets[[1]]$days_since_tavi_dis_plus180)
summary(completed_datasets[[1]]$days_since_tavi_dis_plus180)




###########################   Cox proportional hazards model - time to rehab   ####################################################
###Cox proportional hazards model for Factors associated with cardiac  rehabilitation exposure stratification by region and ethnicity_5_group AND
## stratify by TAVI_procedural_complications given they dont meet assumptions for PH assumption test
##strata(region_name) and strata(ethnicity_5_group) and strata(tavi_complications) allows different baseline hazards , avoiding direct estimation of a regional and ethnicity_5_group effect and tavi_complications.



cox_models_time_to_rehab <- lapply(completed_datasets, function(data) {
  coxph(Surv(cardiac_rehab_time, rehab_status) ~
          age_groups + sex_code + IMD_2019_QUINTILES + 
          BMI_category + binary_known_poor_mobility + KATZ_INDEX_CAT + 
          COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
          NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + SMOKING_STATUS + binary_Known_DIABETES + 
          binary_known_prev_cardiac_surgery + 
          strata(TAVI_procedural_complications) + strata(region_name) + strata (ethnicity_5_group),  ##stratify violating PH variables
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_cox_models_time_to_rehab <- pool(cox_models_time_to_rehab)

# View results
summary(pooled_cox_models_time_to_rehab)


# Extract the summary of the pooled model
summary_df_pooled_cox_models_time_to_rehab <- summary(pooled_cox_models_time_to_rehab)

# Add columns for HR and 95% CI
summary_df_pooled_cox_models_time_to_rehab$HR <- exp(summary_df_pooled_cox_models_time_to_rehab$estimate)
summary_df_pooled_cox_models_time_to_rehab$lower_95_CI <- exp(summary_df_pooled_cox_models_time_to_rehab$estimate - 1.96 * summary_df_pooled_cox_models_time_to_rehab$std.error)
summary_df_pooled_cox_models_time_to_rehab$upper_95_CI <- exp(summary_df_pooled_cox_models_time_to_rehab$estimate + 1.96 * summary_df_pooled_cox_models_time_to_rehab$std.error)

# Round results for clarity
summary_df_pooled_cox_models_time_to_rehab <- summary_df_pooled_cox_models_time_to_rehab[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
summary_df_pooled_cox_models_time_to_rehab <- round(summary_df_pooled_cox_models_time_to_rehab, 3)

# Save to CSV
write.csv(summary_df_pooled_cox_models_time_to_rehab,
          "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/predictors_of_CR_from_imputed_data/summary_df_pooled_cox_models_time_to_rehab.csv",
          row.names = FALSE)






###########################   Multivariable logistic regression model   ####################################################


#Factors associated with cardiac  rehabilitation exposure - a binary analysis (exposed or not exposed to rehab within the follow-up period)

# Fit the logistic regression model, Use region_name as a Random Effect (Mixed-Effects Model)
#Regional variability is important, so use a Generalized Linear Mixed Model (GLMM) instead
#controls for region without estimating a coefficient for each region
#region affects baseline rehab likelihood
#Use if regional differences are important but coefficients are unstable



#Load full mice object and recreate completed datasets
# Load full MICE mids object
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")

# Recreate list of completed datasets
completed_datasets <- lapply(1:5, function(i) complete(imp, i))

# Relevel 'ethnicity_5_group' to have 'White' as the reference level
completed_datasets <- lapply(completed_datasets, function(df) {
  df$ethnicity_5_group <- relevel(factor(df$ethnicity_5_group), ref = "White")
  return(df)
})


install.packages("Matrix", type = "binary")
install.packages("lme4")
library(lme4)
library(Matrix)


# Fit logistic regression models on each imputed dataset (random effects for region)
logistic_model_CR_exposure <- lapply(completed_datasets, function(data) {
  glm(
    rehab_status ~ age_groups + sex_code + ethnicity_5_group + IMD_2019_QUINTILES +
      BMI_category + binary_known_poor_mobility + KATZ_INDEX_CAT +
      COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + SMOKING_STATUS +
      binary_Known_DIABETES + binary_known_prev_cardiac_surgery +
      TAVI_procedural_complications + strata(region_name),
    family = binomial(link = "logit"),
    data = data
  )
})

# Pool using Rubin's rules
pooled_logistic_model_CR_exposure <- pool(logistic_model_CR_exposure)

# View results
summary(pooled_logistic_model_CR_exposure)

# Extract and save results
summary_df_logistic_regression_CR <- summary(pooled_logistic_model_CR_exposure)
summary_df_logistic_regression_CR$OR <- exp(summary_df_logistic_regression_CR$estimate)
summary_df_logistic_regression_CR$lower_CI <- exp(summary_df_logistic_regression_CR$estimate - 1.96 * summary_df_logistic_regression_CR$std.error)
summary_df_logistic_regression_CR$upper_CI <- exp(summary_df_logistic_regression_CR$estimate + 1.96 * summary_df_logistic_regression_CR$std.error)

# Select relevant columns
output_df_summary_df_logistic_regression_CR <- summary_df_logistic_regression_CR[, c("term", "OR", "lower_CI", "upper_CI", "p.value")]

# Round numeric columns (excluding 'term')
output_df_summary_df_logistic_regression_CR[, 2:5] <- round(output_df_summary_df_logistic_regression_CR[, 2:5], 3)

# Save to CSV
write.csv(output_df_summary_df_logistic_regression_CR,
          "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/predictors_of_CR_from_imputed_data/logistic_model_CR_exposure.csv",
          row.names = FALSE
)






############################################################################################################################################################################

##########################################         ==== 6. DOSE RESPONSE ====     #####################################################################################



library(survival)
library(mice)

#Load full mice object and recreate completed datasets
# Load full MICE mids object
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")

# Recreate list of completed datasets
completed_datasets <- lapply(1:5, function(i) complete(imp, i))



# Add days_since_tavi_dis_plus180 back to each imputed dataset # Ensure it has been added to clean_cohort (see code at the start of this script)
# Add the variable back to each imputed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[match(df$person_id, clean_cohort$person_id)]
  return(df)
})

# Confirm it's added correctly in the first dataset
head(completed_datasets[[1]]$days_since_tavi_dis_plus180)
summary(completed_datasets[[1]]$days_since_tavi_dis_plus180)


###################################################################################################################################################################################

############################################################ ==== CR DOSE RESPONSE and RISK of Heart Failure Rehospitalisation ==== ############################################

### Need to filter events within 180 days for each of the outcomes (do each outcome separately)

# Pre-filter all datasets first
filtered_datasets_hf_rehosp <- lapply(completed_datasets, function(df) {
  df[df$hf_rehosp_censor_time > df$days_since_tavi_dis_plus180, ]
})


#Check filter is correct
# Compare number of rows before and after filtering
for (i in 1:length(completed_datasets)) {
  before_n <- nrow(completed_datasets[[i]])
  after_n <- nrow(filtered_datasets_hf_rehosp[[i]])
  cat("Dataset", i, ": Before =", before_n, "| After =", after_n, "| Filtered =", before_n - after_n, "\n")
}

#Filtering is confirmed correct 


#### Fit the adjusted Cox model with continuous rehab attendance - adjusting for covariates - continuous dose response relationship ###
###stratifying ethnicity_5_group and SMOKING_STATUS, these were the 2 problematic variables and once stratified the model runs well. And included REGION


cox_models_hf_rehosp_dose_response_continous_adjusted <- lapply(filtered_datasets_hf_rehosp, function(data) {
  coxph(Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ 
          exp_cardiac_rehab_attends_6m + age_groups + sex_code + region_name + strata(ethnicity_5_group, SMOKING_STATUS) + KATZ_INDEX_CAT + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + 
          IMD_2019_QUINTILES + BMI_category + 
          binary_known_poor_mobility + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
          previous_MI_time_grouped + binary_known_prev_cardiac_surgery + 
          MITRAL_REGURGITATION + TAVI_procedural_complications + 
          binary_Known_DIABETES + LV_FUNCTION + year_tavi_discharge,
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_hf_rehosp_hf_rehosp_dose_response_continous_adjusted <- pool(cox_models_hf_rehosp_dose_response_continous_adjusted)

# View results
summary(pooled_hf_rehosp_hf_rehosp_dose_response_continous_adjusted)


# Extract the summary of the pooled model
summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted <- summary(pooled_hf_rehosp_hf_rehosp_dose_response_continous_adjusted)

# Add columns for HR and 95% CI
summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted$HR <- exp(summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted$estimate)
summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted$lower_95_CI <- exp(summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted$estimate - 1.96 * summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted$std.error)
summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted$upper_95_CI <- exp(summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted$estimate + 1.96 * summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted$std.error)

# Round results for clarity (optional)
summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted <- summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted <- round(summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted, 3)

# Print the final table
print(summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted)


## SAVE OUTPUT ###
# Save to CSV
write.csv(summary_df_hf_rehosp_hf_rehosp_dose_response_continous_adjusted, 
          "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/dose_response_from_imputed_data/hf_rehosp_hf_rehosp_dose_response_continous_adjusted.csv",
          row.names = FALSE)








############################################################################################################################################################################

############################################################====  Dose response CR and Risk of All-Cause Rehospitalisation ==== #################################################################

#re-start the data - need new filter

library(survival)
library(mice)

#Load full mice object and recreate completed datasets
# Load full MICE mids object
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")

# Recreate list of completed datasets
completed_datasets <- lapply(1:5, function(i) complete(imp, i))



# Add days_since_tavi_dis_plus180 back to each imputed dataset # Ensure it has been added to clean_cohort (see code at the start of this script)
# Add the variable back to each imputed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[match(df$person_id, clean_cohort$person_id)]
  return(df)
})



### Need to filter events within 180 days for each of the outcomes (do each outcome separately)

# Pre-filter all datasets first
filtered_datasets_all_cause_rehosp <- lapply(completed_datasets, function(df) {
  df[df$all_cause_rehosp_censor_time > df$days_since_tavi_dis_plus180, ]
})

#Check filter is correct
# Compare number of rows before and after filtering
for (i in 1:length(completed_datasets)) {
  before_n <- nrow(completed_datasets[[i]])
  after_n <- nrow(filtered_datasets_all_cause_rehosp[[i]])
  cat("Dataset", i, ": Before =", before_n, "| After =", after_n, "| Filtered =", before_n - after_n, "\n")
}

#Filtering is confirmed correct 




#### Fit the adjusted Cox model with continuous rehab attendance - adjusting for covariates - continuous dose response relationship ###

# NOTES:  binary_known_poor_mobility and MITRAL_REGURGITATION and IMD do not meet hazard assumption test, therefore stratify these


# ====   dose respnse  all cause rehosp Adjusted for all covarites # ====
#Fit Cox model with calendar time as the timescale, with all covariates 
### ### stratifying binary_known_poor_mobility and MITRAL_REGURGITATION and inlcuded REGION and stratify IMD

cox_models_all_cause_rehosp_dose_response_continuous_fully_adjusted <- lapply(filtered_datasets_all_cause_rehosp, function(data) {
  coxph(Surv(all_cause_rehosp_censor_time, all_cause_rehosp_indicator) ~ 
          exp_cardiac_rehab_attends_6m + age_groups + sex_code + ethnicity_5_group + region_name + SMOKING_STATUS + KATZ_INDEX_CAT + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + 
          BMI_category + 
          strata(binary_known_poor_mobility, MITRAL_REGURGITATION, IMD_2019_QUINTILES) + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
          previous_MI_time_grouped + binary_known_prev_cardiac_surgery + 
          TAVI_procedural_complications + 
          binary_Known_DIABETES + LV_FUNCTION + year_tavi_discharge,
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_all_cause_rehosp_dose_response_continuous_fully_adjusted <- pool(cox_models_all_cause_rehosp_dose_response_continuous_fully_adjusted)

# View results
summary(pooled_all_cause_rehosp_dose_response_continuous_fully_adjusted)


# === Extract summary, calculate HR and 95% CI, and save ===

# Extract the summary
summary_df_all_cause_rehosp_dose_response_continuous_fully_adjusted <- summary(pooled_all_cause_rehosp_dose_response_continuous_fully_adjusted)

# Add columns for HR and 95% CI
summary_df_all_cause_rehosp_dose_response_continuous_fully_adjusted$HR <- exp(summary_df_all_cause_rehosp_dose_response_continuous_fully_adjusted$estimate)
summary_df_all_cause_rehosp_dose_response_continuous_fully_adjusted$lower_95_CI <- exp(summary_df_all_cause_rehosp_dose_response_continuous_fully_adjusted$estimate - 1.96 * summary_df_all_cause_rehosp_dose_response_continuous_fully_adjusted$std.error)
summary_df_all_cause_rehosp_dose_response_continuous_fully_adjusted$upper_95_CI <- exp(summary_df_all_cause_rehosp_dose_response_continuous_fully_adjusted$estimate + 1.96 * summary_df_all_cause_rehosp_dose_response_continuous_fully_adjusted$std.error)

# Round the results
summary_df_all_cause_rehosp_dose_response_continuous_fully_adjusted <- round(summary_df_all_cause_rehosp_dose_response_continuous_fully_adjusted, 3)


## SAVE OUTPUT ###
# Save to CSV
write.csv(summary_df_all_cause_rehosp_dose_response_continuous_fully_adjusted, 
          "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/dose_response_from_imputed_data/all_cause_rehosp_dose_response_continuous_fully_adjusted.csv",
          row.names = FALSE)


###########################################################################################################################################################################


##########################################         ==== 3. Dose response CR and risk of  Non-CVD Rehospitalisation ====     #####################################################################################

#restart the data - need new filter

library(survival)
library(mice)

#Load full mice object and recreate completed datasets
# Load full MICE mids object
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")

# Recreate list of completed datasets
completed_datasets <- lapply(1:5, function(i) complete(imp, i))



# Add days_since_tavi_dis_plus180 back to each imputed dataset # Ensure it has been added to clean_cohort (see code at the start of this script)
# Add the variable back to each imputed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[match(df$person_id, clean_cohort$person_id)]
  return(df)
})



### Need to filter events within 180 days for each of the outcomes (do each outcome separately)

# Pre-filter all datasets first
filtered_datasets_non_cvd_rehosp <- lapply(completed_datasets, function(df) {
  df[df$non_cvd_rehosp_censor_time > df$days_since_tavi_dis_plus180, ]
})

#Check filter is correct
# Compare number of rows before and after filtering
for (i in 1:length(completed_datasets)) {
  before_n <- nrow(completed_datasets[[i]])
  after_n <- nrow(filtered_datasets_non_cvd_rehosp[[i]])
  cat("Dataset", i, ": Before =", before_n, "| After =", after_n, "| Filtered =", before_n - after_n, "\n")
}

#Filtering is confirmed correct 



#### Fit the adjusted Cox model with continuous rehab attendance - adjusting for covariates - continuous dose response relationship ###


# NOTES:  binary_known_poor_mobility and MITRAL_REGURGITATION do not meet hazard assumption test, therefore stratify these

non_cvd_rehosp_cox_model_dose_response_continuous_fully_adjusted <- lapply(filtered_datasets_non_cvd_rehosp, function(data) {
  coxph(Surv(non_cvd_rehosp_censor_time, non_cvd_rehosp_indicator) ~ 
          exp_cardiac_rehab_attends_6m + age_groups + sex_code + ethnicity_5_group + region_name + SMOKING_STATUS + KATZ_INDEX_CAT + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + 
          IMD_2019_QUINTILES + BMI_category + 
          strata(binary_known_poor_mobility, MITRAL_REGURGITATION) + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
          previous_MI_time_grouped + binary_known_prev_cardiac_surgery + 
          TAVI_procedural_complications + 
          binary_Known_DIABETES + LV_FUNCTION + year_tavi_discharge, 
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_non_cvd_rehosp_cox_model_dose_response_continuous_fully_adjusted <- pool(non_cvd_rehosp_cox_model_dose_response_continuous_fully_adjusted)

# View results
summary(pooled_non_cvd_rehosp_cox_model_dose_response_continuous_fully_adjusted)




# === Extract summary of pooled model ===
summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted <- summary(pooled_non_cvd_rehosp_cox_model_dose_response_continuous_fully_adjusted)

# === Add HR and 95% CI columns ===
summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted$HR <- exp(summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted$estimate)
summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted$lower_95_CI <- exp(summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted$estimate - 1.96 * summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted$std.error)
summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted$upper_95_CI <- exp(summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted$estimate + 1.96 * summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted$std.error)

# === Select and round final output ===
summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted <- summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted <- round(summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted, 3)

# === Save as CSV ===
write.csv(
  summary_df_non_cvd_rehosp_cox_model_dose_response_continuous_full_adjusted,
  file = "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/dose_response_from_imputed_data/non_cvd_rehosp_dose_response_continuous_fully_adjusted.csv",
  row.names = FALSE
)


########################################################################################################################################################################################################


#####################################################   ==== Dose response CR appointments attended and risk of MORTALITY     ==== #####################################################################

#restart the data - need new filter

library(survival)
library(mice)

#Load full mice object and recreate completed datasets
# Load full MICE mids object
imp <- readRDS("D:/PhotonUser/My Files/Home Folder/Justin Braver/Data/Imputed_data/imp_mice_object.rds")

# Recreate list of completed datasets
completed_datasets <- lapply(1:5, function(i) complete(imp, i))



# Add days_since_tavi_dis_plus180 back to each imputed dataset # Ensure it has been added to clean_cohort (see code at the start of this script)
# Add the variable back to each imputed dataset
completed_datasets <- lapply(completed_datasets, function(df) {
  df$days_since_tavi_dis_plus180 <- clean_cohort$days_since_tavi_dis_plus180[match(df$person_id, clean_cohort$person_id)]
  return(df)
})




### Need to filter events within 180 days for each of the outcomes (do each outcome separately)

# Pre-filter all datasets first
filtered_datasets_mortality <- lapply(completed_datasets, function(df) {
  df[df$all_cause_mortality_censor_time > df$days_since_tavi_dis_plus180, ]
})

#Check filter is correct
# Compare number of rows before and after filtering
for (i in 1:length(completed_datasets)) {
  before_n <- nrow(completed_datasets[[i]])
  after_n <- nrow(filtered_datasets_mortality[[i]])
  cat("Dataset", i, ": Before =", before_n, "| After =", after_n, "| Filtered =", before_n - after_n, "\n")
}

#Filtering is confirmed correct 




#  Fit the adjusted Cox model with continuous rehab attendance - adjusting for covariates - continuous dose response relationship ###
# NOTES: tavi complications and year of discharge do not meet hazard assumption test, therefore stratify these


mortality_cox_model_dose_response_continuous_fully_adjusted <- lapply(filtered_datasets_mortality, function(data) {
  coxph(Surv(all_cause_mortality_censor_time, all_cause_mortality_indicator) ~ 
          exp_cardiac_rehab_attends_6m + age_groups + sex_code + ethnicity_5_group + region_name + SMOKING_STATUS + KATZ_INDEX_CAT + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + 
          IMD_2019_QUINTILES + BMI_category + binary_known_poor_mobility + MITRAL_REGURGITATION +
          strata(TAVI_procedural_complications, year_tavi_discharge) + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
          previous_MI_time_grouped + binary_known_prev_cardiac_surgery + 
          binary_Known_DIABETES + LV_FUNCTION,
        data = data)
})

# Pool using Rubin's rules (mice)
pooled_mortality_cox_model_dose_response_continuous_fully_adjusted <- pool(mortality_cox_model_dose_response_continuous_fully_adjusted)

# View results
summary(pooled_mortality_cox_model_dose_response_continuous_fully_adjusted)



# Extract the summary of the pooled model
summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted <- summary(pooled_mortality_cox_model_dose_response_continuous_fully_adjusted)

# Add columns for HR and 95% CI
summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted$HR <- exp(summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted$estimate)
summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted$lower_95_CI <- exp(summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted$estimate - 1.96 * summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted$std.error)
summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted$upper_95_CI <- exp(summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted$estimate + 1.96 * summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted$std.error)

# Round results
summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted <- summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted[, c("term", "HR", "lower_95_CI", "upper_95_CI", "p.value")]
summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted <- round(summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted, 3)

# Save to CSV
write.csv(summary_df_mortality_cox_model_dose_response_continuous_fully_adjusted, 
          "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Imputed data summary/dose_response_from_imputed_data/mortality_cox_model_dose_response_continuous_fully_adjusted.csv",
          row.names = FALSE)




##############################################################################################################################################################################################################








