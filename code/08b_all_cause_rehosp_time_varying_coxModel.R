#########################################################################################################################################
#####                            OBJECTIVE 4. Effectiveness of cardiac rehabilitation  post TAVI                          ###############
#####      SECONDARY outcome - all cause hospital readmission, cox model with CR time varying exposure   ###############################
#########################################################################################################################################


##Cox model with time-varying exposure to cardiac rehabilitation (CR) using calendar time as the timescale, 
##The observation period is from 1 January 2018 to 31 March 2024
#The model will evaluate how the patient's hazard for ALL CAUSE rehospitalization based on their participation in rehab during different time intervals

#####################################################
##Step 1: Define and create Calendar Time Variables##
#####################################################

# Load necessary libraries
library(survival)
library(dplyr)
library(lubridate)

#Create Relevant Time Variables for:
# tavi_time: days since 1 January 2018 to the date of TAVI discharge
# rehab_time: days since 1 January 2018 to the start of cardiac rehab (if applicable).
# rehosp_censor_time: days since 1 January 2018 to all cause rehospitalization or censoring (death or end of study period)

# Define study start and end dates
study_start <- as.Date("2018-01-01")
study_end <- as.Date("2024-03-31")

# Add necessary time variables
clean_cohort_all_cause_rehosp <- clean_cohort %>%
  mutate(
    # Days since 1 January 2018 (calendar time)
    days_since_2018 = as.numeric(tavi_cips_disdate - study_start),  #this will be used in step2 (as time_1 in the split data with time intervals)
    
    # TAVI time (days from 1 January 2018 to TAVI procedure discharge date)
    tavi_time = as.numeric(tavi_cips_disdate - study_start),
    
    # Rehab time (days from 1 January 2018 to the start of the 1st rehab appointment)
    rehab_time = ifelse(!is.na(exp_cardiac_rehab_appt_date_6m), 
                        as.numeric(exp_cardiac_rehab_appt_date_6m - study_start), NA),
    
    # Time of all cause rehospitalization or censoring (death or end of study period)
    rehosp_censor_time_all_cause = pmin(as.numeric(out_readm_all_cause_date - study_start),            # Days to all cause readmission  
                              as.numeric(date_of_death - study_start),                # Days to death
                              as.numeric(study_end - study_start), na.rm = TRUE),     # Days to study end
    
    # Rehab exposure indicator
    rehab_indicator = ifelse(is.na(exp_cardiac_rehab_6m) | exp_cardiac_rehab_6m == 0, 0, 1),
    
    # Rehospitalization indicator (1 for all cause rehospitalization, 0 for censoring)
    rehosp_indicator_all_cause = ifelse(is.na(out_readm_all_cause_flag) | out_readm_all_cause_flag == 0, 0, 1) 
    
  )
    
  

##check if any patients have rehosp or censoring time less than or equal to TAVI time and rehosps or death within 180 days of TAVI discharge
count_invalid_all_cause_rehosp <- sum(clean_cohort_all_cause_rehosp$rehosp_censor_time_all_cause <= clean_cohort_all_cause_rehosp$tavi_time + 180, na.rm = TRUE
)
print(count_invalid_all_cause_rehosp)


##NOTES, n=8915 (rounded to nearest 5)  patients have an all cause rehosp or censoring time less than or equal to TAVI time and/or less than tavi_time+180 days, only interested in post discharge from TAVI and post TAVI + 180 days  
##only interested in post discharge from TAVI plus 180 days to account for enough time post discharge to enrol in rehab,  so filter out these patients within 180 days post discharge



# Grouping by rehab_indicator and summarizing count
count_invalid_all_cause_rehosp <- clean_cohort_all_cause_rehosp %>%
  group_by(rehab_indicator) %>%  # Group by rehab exposure
  summarise(
    count_invalid = sum(
      rehosp_censor_time_all_cause <= tavi_time | 
        rehosp_censor_time_all_cause <= tavi_time + 180, na.rm = TRUE
    ) 
  )

# Print results
print(count_invalid_all_cause_rehosp)


### 8541 patients in the non CR group and 376 patients in the CR group have all cause rehosp within 180 days of TAVI discharge - therefore filter them out




# Filter out patients based on rehosp_censor_time conditions, filter any all-cause rehosps within 180 days post discharge
clean_cohort_all_cause_rehosp <- clean_cohort_all_cause_rehosp %>%
  filter(!(rehosp_censor_time_all_cause < tavi_time + 180))

# Count the number of invalid cases 
count_invalid_all_cause_rehosp <- sum(clean_cohort_all_cause_rehosp$rehosp_censor_time_all_cause <= clean_cohort_all_cause_rehosp$tavi_time| 
                                        clean_cohort_all_cause_rehosp$out_readm_all_cause_date < clean_cohort_all_cause_rehosp$tavi_time + 180, na.rm = TRUE)
print(count_invalid_all_cause_rehosp)


###NOTES: FILTERED OUT 8915 (rounded to nearest 5) PATIENTS who had a rehosp or censoring time less than or equal to TAVI time and or patients with a rehosp within 180 days of tavi discharge
###NOTE: 250 (rounded to nearest 5) patients had a rehosp censoring time less than TAVI time and the rest had a rehosp censor time less than TAVI_time+180 days


#####################################################################
##Step 2: Split Data into Time Intervals for Time-Varying Cox Model##
#####################################################################

#To model time-varying exposure, split data based on whether the patient was in CR (rehab_time):
#Pre-Rehab: from TAVI to the start of CR (if applicable).
#Post-Rehab: from start of CR to rehospitalization or censoring.


split_data_all_cause_rehosp <- clean_cohort_all_cause_rehosp |>                   #Part 1: Splitting the Dataset into Periods Based on Rehab Time
  mutate(
    period_start_1 = days_since_2018,           # The start of the first period is the start of the follow-up time.
    period_end_1 = case_when(
      is.na(rehab_time) ~ rehosp_censor_time_all_cause,   # If the rehab time is NA, the period ends at rehospitalization or censoring
      .default = rehab_time                     # Otherwise, the period ends at rehab start time.
    ),
    period_start_2 = case_when(
      !is.na(rehab_time) ~ rehab_time,          # The second period starts at rehab time if rehab time is not NA.
      .default = NA                             # Otherwise, no second period.
    ),
    period_end_2 = case_when(
      period_end_1 == rehosp_censor_time_all_cause ~ NA,  # If the first period already ends at censoring, no second period.
      .default = rehosp_censor_time_all_cause             # Otherwise, the second period ends at rehosp or censoring.
    )
  ) |>                                          #Part 2: Reshape the Split Data for Modeling
  dplyr::select(matches("person_id|^period")) |>       # Select person_id and the period columns.
  pivot_longer(cols = !person_id) |>            # Reshape to a long format, gathering columns into key-value pairs.
  drop_na(value) |>                             # Remove missing values.
  mutate(
    period = stringr::str_extract(name, "[1|2]"),     # Extract the period number (1 or 2).
    name = gsub("_[1|2]", "", name)                   # Remove the period number from the variable name.
  ) |>
  pivot_wider(id_cols = c(person_id, period), names_from = name, values_from = value)      # Reshape back to wide format.


#Part 3: Joining and Adding Covariates

split_data2_all_cause_rehosp <- split_data_all_cause_rehosp |> left_join(clean_cohort_all_cause_rehosp |> dplyr::select(person_id, days_since_2018, tavi_time, rehab_time, rehosp_censor_time_all_cause, rehab_indicator,
                                                                                                                        rehosp_indicator_all_cause, age_at_tavi, age_groups, sex_code, ethnicity_5_group, IMD_2019_QUINTILES, region_name, 
                                                                                                                        SMOKING_STATUS, KATZ_INDEX_CAT, KATZ_INDEX_CLEAN, COMORBIDITY_pulmonary_liver_neuro_excardiacvascu,
                                                                                                                        TAVI_procedural_complications,
                                                                                                                        BMI, BMI_category,obesity_binary, NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE,
                                                                                                                        CCS_ANGINA_STATUS_PRE_PRO_STABLE,
                                                                                                                        CSHA_CLINICAL_FRAILTY_SCALE_SCORE, binary_known_extracardiac_arteriopathy, 
                                                                                                                        binary_known_poor_mobility, binary_Known_DIABETES, binary_known_on_dialysis, 
                                                                                                                        BINARY_known_SEVERE_LIVER_DISEASE,
                                                                                                                        binary_known_pulmonary_disease,
                                                                                                                        BINARY_known_HISTORY_OF_NEUROLOGICAL_DISEASE, EXTENT_OF_CORONARY_VESSEL_DISEASE,LV_FUNCTION,
                                                                                                                        binary_known_prev_mi, previous_MI_time_grouped, binary_known_prev_cardiac_surgery,
                                                                                                                        binary_known_prev_pci,
                                                                                                                        binary_previous_known_CABG, binary_previous_known_valve_operation,
                                                                                                                        binary_previous_known_other_operation_pericardium, MITRAL_REGURGITATION,
                                                                                                                        PA_SYSTOLIC_PRESSURE_MMHG, 
                                                                                                                        BINARY_known_PERMANENT_PACING, BINARY_known_BLEEDING,
                                                                                                                        BINARY_known_VASCULAR_ACCESS_SITE_COMPLICATIONS,
                                                                                                                        BINARY_known_USE_OF_CARDIOPULMONARY_BYPASS,
                                                                                                                        BINARY_known_ACUTE_KIDNEY_INJURY,
                                                                                                                        exp_cardiac_rehab_appts_6m, exp_cardiac_rehab_attends_6m, rehab_attendance_groups)) |>
  mutate(
    rehab_indicator = if_else(period == "1", 0, rehab_indicator),     # Rehab indicator is 0 before rehab.
    dupes = if_else(period == 2, 1, NA)                               # Mark duplicates
  ) |> 
  group_by(person_id) |>                                              # Group by person
  fill(dupes, .direction = "up") |>                                   # Propagate the duplicate marker upward.
  ungroup() |>                                                        # Remove grouping.
  mutate(rehosp_indicator_all_cause = case_when(period == "1" & dupes == 1 ~ 0, .default = rehosp_indicator_all_cause)) |>         # Correct the rehosp indicator
  dplyr::select(-dupes)



#####################################################################
##Step 3: Fit Time-Varying Cox Proportional Hazards Model############
#####################################################################
#Fit the Cox proportional hazards model with the coxph function

# Fit time-varying Cox model with calendar time as the timescale
cox_model_timevaryingCR_all_cause_rehosp_unadjusted <- coxph(
  Surv(period_start, period_end, rehosp_indicator_all_cause) ~ rehab_indicator,
  data = split_data2_all_cause_rehosp
)

# Summary of the Cox model
summary(cox_model_timevaryingCR_all_cause_rehosp_unadjusted)


#count number of patients who had rehab after the rehosp (not after TAVI), therefore filter out these patients
count_invalid_all_cause_rehosp2 <- sum(split_data2_all_cause_rehosp$period_end <= split_data2_all_cause_rehosp$period_start, na.rm = TRUE)
print(count_invalid_all_cause_rehosp2)


### NOTES: n= 1 patient had rehab after the rehosp (not after TAVI), therefore filter out these patients
split_data2_all_cause_rehosp <- split_data2_all_cause_rehosp %>%
  filter(!(period_end <= period_start))

sum(split_data2_all_cause_rehosp$period_end <= split_data2_all_cause_rehosp$period_start, na.rm = TRUE)


# Fit time-varying Cox model with calendar time as the timescale, with 1 patient filtered out
cox_model_timevaryingCR_all_cause_rehosp_unadjusted_2 <- coxph(
  Surv(period_start, period_end, rehosp_indicator_all_cause) ~ rehab_indicator,
  data = split_data2_all_cause_rehosp
)

# Summary of the Cox model
summary(cox_model_timevaryingCR_all_cause_rehosp_unadjusted_2)


##Note p-value is 0.05, CR exposure is sig associated with rate of all cause rehosp compared to not exposed 


#####################################################################
##Step 4: Check Assumptions  ################################
#####################################################################

#Check proportional hazards assumption
# Test proportional hazards assumption
cox.zph(cox_model_timevaryingCR_all_cause_rehosp_unadjusted_2)

# rehab_indicator p-values are higher than 0.05 therefore hazard assumptions hold true. Need to add all other covarites




##extract unadjusted time varying CR exposure cox modle for all cause rehosp results 

library(broom)
library(dplyr)
library(readr)

# 1. Tidy model output
results_cox_model_timevaryingCR_all_cause_rehosp_unadjusted_2 <- tidy(
  cox_model_timevaryingCR_all_cause_rehosp_unadjusted_2,
  exponentiate = TRUE,   # gives HR instead of log(HR)
  conf.int = TRUE        # includes 95% CI
)

# 2. Clean column names and round
results_cox_model_timevaryingCR_all_cause_rehosp_unadjusted_2_clean <- results_cox_model_timevaryingCR_all_cause_rehosp_unadjusted_2 %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))


# Save the rounded table as a CSV file for extraction
output_file_results_cox_model_timevaryingCR_all_cause_rehosp_unadjusted_2_clean <- 
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_rehosp/all_cause_rehosp_CR_time_varying_cox_model/results_cox_model_timevaryingCR_all_cause_rehosp_unadjusted_2_clean.csv"
write.csv(results_cox_model_timevaryingCR_all_cause_rehosp_unadjusted_2_clean, file = output_file_results_cox_model_timevaryingCR_all_cause_rehosp_unadjusted_2_clean, row.names = FALSE)







##########################################################################################################################
##Adjusted model: All Covariates  ####################################################
##########################################################################################################################

#Covariate selection by study cardiologist Professor Tom Marwick
#   Complications after TAVI likely dependent on: age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
#  BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
#  previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications  


#####################################################################
## Model 3: Fit Time-Varying Cox Proportional Hazards Model############
#####################################################################
#Fit the Cox proportional hazards model with the coxph function

#Fit time-varying Cox model with calendar time as the timescale, with all covariates
cox_model_timevaryingCR_all_cause_rehosp_3a <- coxph(
  Surv(period_start, period_end, rehosp_indicator_all_cause) ~ rehab_indicator + age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
    BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications,
  data = split_data2_all_cause_rehosp
)

# Summary of the Cox model
summary(cox_model_timevaryingCR_all_cause_rehosp_3a)


#Check proportional hazards assumption
# Test proportional hazards assumption
cox.zph(cox_model_timevaryingCR_all_cause_rehosp_3a)


# rehab_indicator and Katz_index for the hazards assumptions test are less than 0.05
#therefore hazard assumptions do not hold true. Need to adjust these variables 
# Stratify Katz_index Categorical Variables That Violate PH:  stratification removes the need to estimate their hazard ratio over time.



##########################################################################################################################
##Adjusted model 3b: 
### ###Stratify Katz index to correct for hazrad assumptions in model above
##########################################################################################################################


cox_model_timevaryingCR_all_cause_rehosp_fully_adjusted <- coxph(
  Surv(period_start, period_end, rehosp_indicator_all_cause) ~ 
    rehab_indicator  + 
    age_at_tavi + sex_code + ethnicity_5_group +  
    binary_Known_DIABETES + 
    NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + strata(KATZ_INDEX_CAT, LV_FUNCTION) +  # Stratifying PH-violating factors, included LV function which improves the model 
    binary_known_prev_cardiac_surgery +  
    MITRAL_REGURGITATION + TAVI_procedural_complications + 
    SMOKING_STATUS + BMI_category + previous_MI_time_grouped +
  IMD_2019_QUINTILES + binary_known_poor_mobility + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu,  
  data = split_data2_all_cause_rehosp
)

# Summary of the Cox model
summary(cox_model_timevaryingCR_all_cause_rehosp_fully_adjusted)


#Check proportional hazards assumption
# Test proportional hazards assumption
cox.zph(cox_model_timevaryingCR_all_cause_rehosp_fully_adjusted)



#NOTES:Hazard assumptions hold true. all varaibles have p-value greater than 0.05 and global p-value = 0.331  
## # NOTES INTERPRETATION OF MAIN FINDING: after Adjusting for all key covariates, CR exposure is not associated with slower rate of all-cause rehosps compared with no CR exposure




##extract fully adjusted time varying CR exposure cox model for all cause rehosp results 

library(broom)
library(dplyr)
library(readr)

# 1. Tidy model output
results_cox_model_timevaryingCR_all_cause_rehosp_fully_adjusted <- tidy(
  cox_model_timevaryingCR_all_cause_rehosp_fully_adjusted,
  exponentiate = TRUE,   # gives HR instead of log(HR)
  conf.int = TRUE        # includes 95% CI
)

# 2. Clean column names and round
results_cox_model_timevaryingCR_all_cause_rehosp_fully_adjusted_clean <- results_cox_model_timevaryingCR_all_cause_rehosp_fully_adjusted %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))


# Save the rounded table as a CSV file for extraction
output_file_results_cox_model_timevaryingCR_all_cause_rehosp_fully_adjusted_clean <- 
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_rehosp/all_cause_rehosp_CR_time_varying_cox_model/results_cox_model_timevaryingCR_all_cause_rehosp_fully_adjusted_clean.csv"
write.csv(results_cox_model_timevaryingCR_all_cause_rehosp_fully_adjusted_clean, file = output_file_results_cox_model_timevaryingCR_all_cause_rehosp_fully_adjusted_clean, row.names = FALSE)







#########################################################################################################################################
#####                            SUMMARY OF REHOSP VOLUMES - all cause rehosps                                     ###############
#########################################################################################################################################


##NOTE: this is based on the filtered out patients from above (n=8915 - rounded to nearest 5), these patients had rehop within 180 days of TAVI - reset the data if you do not want the filtered baseline cohort. 


##raw volumes of all cause hospitalizations of entire cohort since TAVI + 180 days until study end
##hospitilisations are only rehosps post 180days of TAVI discharge

all_cause_rehosp_raw_vols <- clean_cohort_all_cause_rehosp |>
  mutate(tavi_cips_disdate_plus180_days = tavi_cips_disdate + days(180),  
         rehosp_indicator_all_cause = ifelse((out_readm_all_cause_date >= tavi_cips_disdate_plus180_days), 1,  0) 
         )
all_cause_rehosp_raw_vols$rehosp_indicator_all_cause[is.na(all_cause_rehosp_raw_vols$rehosp_indicator_all_cause)] <- 0

all_cause_rehosp_raw_vols_summary <- all_cause_rehosp_raw_vols |>
  group_by(rehab_indicator) |>    # Group by rehab exposure
  summarise(
    total_patients = round(n_distinct(person_id) / 5) * 5,      # Unique patients in each group and Round to nearest 5
    total_all_cause_rehosp = round(sum(rehosp_indicator_all_cause) / 5) * 5,  # Total all-cause rehospitalizations and # Round to nearest 5
    avg_all_cause_rehosp_per_patient = total_all_cause_rehosp / total_patients, # Avg all-cause rehosp per patient
  )
print(all_cause_rehosp_raw_vols_summary)  



# Save the rounded table as a CSV file for extraction
output_file_all_cause_rehosp_raw_vols_summary <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_rehosp/all_cause_rehosp_raw_vols_summary.csv"
write.csv(all_cause_rehosp_raw_vols_summary, file = output_file_all_cause_rehosp_raw_vols_summary, row.names = FALSE)

##note this excludes the volume of patients (n=8915 - rounded to nearest 5) who had an all cause rehosp within 180 days post TAVI, that's why total patients coloumn is not the full baseline cohort








#########################################################################################################################################












