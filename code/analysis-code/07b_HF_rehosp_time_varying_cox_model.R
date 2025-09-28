#####################################################################################################################################################################################################
##  HF rehospitalization (rehosp) - filerting  HF rehosps within 180 days of TAVi Discharge # 
## Therefore model will be Time since discharge from TAVI plus 180 days                                                                                                                           ###
#################################################################################################################################################################################################

##Filter out patinets rehospitalized within 180 days after the TAVI discharge during which patients should not be penalized for not attending rehab.


# Load necessary libraries
library(survival)
library(dplyr)
library(lubridate)

#Create Relevant Time Variables for:
# tavi_time: days since 1 January 2018 to the date of TAVI discharge
# rehab_time: days since 1 January 2018 to the start of cardiac rehab (if applicable).
# rehosp_censor_time: days since 1 January 2018 to rehospitalization for HF or censoring (death or end of study period)

# Define study start and end dates
study_start <- as.Date("2018-01-01")
study_end <- as.Date("2024-03-31")

# Add necessary time variables
clean_cohort <- clean_cohort %>%
  mutate(
    # Days since 1 January 2018 (calendar time)
    days_since_2018 = as.numeric(tavi_cips_disdate - study_start),  #this will be used in step2 (as time_1 in the split data with time intervals)
    
    # TAVI time (days from 1 January 2018 to TAVI procedure discharge date)
    tavi_time = as.numeric(tavi_cips_disdate - study_start),
    
    # Rehab time (days from 1 January 2018 to the start of the 1st rehab appointment)
    rehab_time = ifelse(!is.na(exp_cardiac_rehab_appt_date_6m), 
                        as.numeric(exp_cardiac_rehab_appt_date_6m - study_start), NA),
    
    # Time of HF rehospitalization or censoring (death or end of study period)
    rehosp_censor_time = pmin(as.numeric(out_readm_hf_date - study_start),            # Days to HF readmission
                              as.numeric(date_of_death - study_start),                # Days to death
                              as.numeric(study_end - study_start), na.rm = TRUE),     # Days to study end
    
    # Rehab exposure indicator
    rehab_indicator = ifelse(is.na(exp_cardiac_rehab_6m) | exp_cardiac_rehab_6m == 0, 0, 1),
    
    # Rehospitalization indicator (1 for HF rehospitalization, 0 for censoring)
    rehosp_indicator = ifelse(is.na(out_readm_hf_flag) | out_readm_hf_flag == 0, 0, 1)
    
  )




##check if any patients have rehosp or censoring time less than or equal to TAVI time and rehosps within 180 days of TAVI discharge
count_invalid_HF_rehosp <- sum(clean_cohort$rehosp_censor_time < clean_cohort$tavi_time + 180, na.rm = TRUE)
print(count_invalid_HF_rehosp)



# Grouping by rehab_indicator and summarizing count
count_invalid_HF_rehosp <- clean_cohort %>%
  group_by(rehab_indicator) %>%  # Group by rehab exposure
  summarise(
    count_invalid = sum(
      rehosp_censor_time <= tavi_time | 
        rehosp_censor_time < tavi_time + 180, na.rm = TRUE
    ) 
  )

# Print results
print(count_invalid_HF_rehosp)



##NOTES, just over 2K patients have a HF rehosp or censoring time less than or equal to TAVI time and/or have a rehosp within 180 days of TAVI discharge
## 2105(rounded to nearest 5) patients not exposed to rehab vs 80 (rounded to nearest 5) patients exposed to rehab
##only interested in post discharge from TAVI so filter out these patients 
##only interested in post discharge from TAVI plus 180 days to account for enough time post discharge to enrol in rehab,  so filter out these patients within 180 days post discharge



# Filter out these patients 
clean_cohort <- clean_cohort %>%
  filter(!(rehosp_censor_time <= tavi_time | rehosp_censor_time <= days_since_2018 | rehosp_censor_time < tavi_time + 180))

# Count the number of invalid cases 
count_invalid_HF_rehosp <- sum(clean_cohort$rehosp_censor_time <= clean_cohort$tavi_time| 
                                 clean_cohort$rehosp_censor_time < clean_cohort$tavi_time + 180, na.rm = TRUE)
print(count_invalid_HF_rehosp)


###NOTES: FILTERED OUT n=2185 (rounded) PATIENTS (2105 in no rehab and 80 in the rehab group - rounded to nearest 5) who had a rehosp or censoring time less than or equal to TAVI time and or patients with a rehosp censoring time within 180 days of tavi discharge






#####################################################################
##Step 2: Split Data into Time Intervals for Time-Varying Cox Model##
#####################################################################

#To model time-varying exposure, split data based on whether the patient was in CR (rehab_time):
#Pre-Rehab: from TAVI to the start of CR (if applicable).
#Post-Rehab: from start of CR to rehospitalization or censoring.


split_data <- clean_cohort |>                   #Part 1: Splitting the Dataset into Periods Based on Rehab Time
  mutate(
    period_start_1 = days_since_2018,           # The start of the first period is the start of the follow-up time.
    period_end_1 = case_when(
      is.na(rehab_time) ~ rehosp_censor_time,   # If the rehab time is NA, the period ends at rehospitalization or censoring
      .default = rehab_time                     # Otherwise, the period ends at rehab start time.
    ),
    period_start_2 = case_when(
      !is.na(rehab_time) ~ rehab_time,          # The second period starts at rehab time if rehab time is not NA.
      .default = NA                             # Otherwise, no second period.
    ),
    period_end_2 = case_when(
      period_end_1 == rehosp_censor_time ~ NA,  # If the first period already ends at censoring, no second period.
      .default = rehosp_censor_time             # Otherwise, the second period ends at rehosp or censoring.
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

split_data2 <- split_data |> left_join(clean_cohort |> dplyr::select(person_id, days_since_2018, tavi_time, rehab_time, rehosp_censor_time, rehab_indicator, rehosp_indicator, age_at_tavi, age_groups, sex_code, ethnicity_5_group, IMD_2019_QUINTILES, region_name, 
                                                                     SMOKING_STATUS, KATZ_INDEX_CAT, KATZ_INDEX_CLEAN, COMORBIDITY_pulmonary_liver_neuro_excardiacvascu, TAVI_procedural_complications,
                                                                     BMI, BMI_category,obesity_binary, NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE, CCS_ANGINA_STATUS_PRE_PRO_STABLE,
                                                                     CSHA_CLINICAL_FRAILTY_SCALE_SCORE, binary_known_extracardiac_arteriopathy, 
                                                                     binary_known_poor_mobility, binary_Known_DIABETES, binary_known_on_dialysis, BINARY_known_SEVERE_LIVER_DISEASE,
                                                                     binary_known_pulmonary_disease,
                                                                     BINARY_known_HISTORY_OF_NEUROLOGICAL_DISEASE, EXTENT_OF_CORONARY_VESSEL_DISEASE,LV_FUNCTION,
                                                                     binary_known_prev_mi, previous_MI_time_grouped, binary_known_prev_cardiac_surgery, binary_known_prev_pci,
                                                                     binary_previous_known_CABG, binary_previous_known_valve_operation,
                                                                     binary_previous_known_other_operation_pericardium, MITRAL_REGURGITATION, PA_SYSTOLIC_PRESSURE_MMHG, 
                                                                     BINARY_known_PERMANENT_PACING, BINARY_known_BLEEDING, BINARY_known_VASCULAR_ACCESS_SITE_COMPLICATIONS,
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
  mutate(rehosp_indicator = case_when(period == "1" & dupes == 1 ~ 0, .default = rehosp_indicator)) |>         # Correct the rehosp indicator
  dplyr::select(-dupes)



####################################################################
##Step 3: Fit Time-Varying Cox Proportional Hazards Model############
#####################################################################
#Fit the Cox proportional hazards model with the coxph function

# Fit time-varying Cox model with calendar time as the timescale
cox_model_timevaryingCR_HF_primary_outcome <- coxph(
  Surv(period_start, period_end, rehosp_indicator) ~ rehab_indicator,
  data = split_data2
)

# Summary of the Cox model
summary(cox_model_timevaryingCR_HF_primary_outcome)


##Notes: warning appears that there are patients who had rehab after the rehosp (not after TAVI)

# Count the number of invalid cases where patients had rehab after the rehosp (not after TAVI) 
count_invalid_cox_HF_model_1_rehab_post_rehosp <- sum(split_data2$period_end <= split_data2$period_start, na.rm = TRUE)

print(count_invalid_cox_HF_model_1_rehab_post_rehosp)



### NOTES: n= 1 patient had rehab after the rehosp (not after TAVI), therefore filter out this patient
split_data2 <- split_data2 %>%
  filter(!(period_end <= period_start))



# Fit time-varying Cox model with calendar time as the timescale - starting fom 180 days from TAVI discharge, with 1 patients filtered out who did rehab before rehosp
cox_model_timevaryingCR_HF_primary_outcome_2 <- coxph(
  Surv(period_start, period_end, rehosp_indicator) ~ rehab_indicator,
  data = split_data2
)

# Summary of the Cox model
summary(cox_model_timevaryingCR_HF_primary_outcome_2)



#####################################################################
##Step 4: Check Assumptions 
#####################################################################

# Test proportional hazards assumption
cox.zph(cox_model_timevaryingCR_HF_primary_outcome_2)


##Notes: model meets hazard assumption test for all varaibles. Now add all covariates to the model



#extract results and save as CSV
library(broom)
library(dplyr)
library(readr)

# Extract tidy results with HRs and confidence intervals
cox_model_timevaryingCR_HF_unadjusted_results <- tidy(
  cox_model_timevaryingCR_HF_primary_outcome_2,
  exponentiate = TRUE,
  conf.int = TRUE
)

# Clean and round output
cox_model_timevaryingCR_HF_unadjusted_results_clean <- cox_model_timevaryingCR_HF_unadjusted_results %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))


# Save the rounded table as a CSV file for extraction
output_file_cox_model_timevaryingCR_HF_unadjusted_results <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/HF Rehosp/Time varying cox model/cox_model_timevaryingCR_HF_unadjusted_results.csv"
write.csv(cox_model_timevaryingCR_HF_unadjusted_results_clean, file = output_file_cox_model_timevaryingCR_HF_unadjusted_results, row.names = FALSE)







##########################################################################################################################
##Adjusted model 2a: All Covariates  ####################################################
##########################################################################################################################

#Covariate selection by study cardiologist Professor Tom Marwick
#   Complications after TAVI likely dependent on: age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
#  BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
#  previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications  

#Note: did not include PA_SYSTOLIC_PRESSURE_MMHG b/c many missing observations,and was not a key varaible in the clinical selection process 

#############################################################
# Fit Time-Varying Cox Proportional Hazards Model 2a ########
#############################################################
#Fit the Cox proportional hazards model with the coxph function


# Fit time-varying Cox model with calendar time as the timescale, with all covariates
cox_model_timevaryingCR_HF_primary_outcome_model_2a <- coxph(
  Surv(period_start, period_end, rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
    BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications  
  ,
  data = split_data2
)

# Summary of the Cox model
summary(cox_model_timevaryingCR_HF_primary_outcome_model_2a)


# NOTES INTERPRETATION OF MAIN FINDING: Adjusting for all key covariates, CR exposure is not associated with a lower rate of HF rehosp compared to not being exposed to CR


# Check Assumptions
#Check proportional hazards assumption
# Test proportional hazards assumption
cox.zph(cox_model_timevaryingCR_HF_primary_outcome_model_2a)


###NOTES: Model runs with all covariates. 
## binary_known_poor_mobility and binary_known_prev_cardiac_surgery are less than 0.05 therefore proportional hazards assumptions do not hold for these. 
##Stratify by binary_known_poor_mobility and binary_known_prev_cardiac_surgery, to correct for hazrad assumptions 







##############################################################################################################################################################################
##Adjusted model 2B: All Covariates and stratify binary_known_poor_mobility and binary_known_prev_cardiac_surgery, to correct for hazrad assumptions #########################
##############################################################################################################################################################################



# Fit time-varying Cox model with calendar time as the timescale, with all covariates and stratify binary_known_poor_mobility and binary_known_prev_cardiac_surgery
cox_model_timevaryingCR_HF_primary_outcome_model_2b <- coxph(
  Surv(period_start, period_end, rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
    BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + binary_Known_DIABETES +
    previous_MI_time_grouped + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications +
    strata(binary_known_poor_mobility, binary_known_prev_cardiac_surgery)   #stratify PH violating variables
  ,
  data = split_data2
)

# Summary of the Cox model
summary(cox_model_timevaryingCR_HF_primary_outcome_model_2b)


# NOTES INTERPRETATION OF MAIN FINDING: Adjusting for all key covariates, CR exposure is not sig associated with the rate of HF rehosp 


# Check Assumptions
#Check proportional hazards assumption
# Test proportional hazards assumption
cox.zph(cox_model_timevaryingCR_HF_primary_outcome_model_2b)

##NOTES: all proportional hazards assumptions hold true - model functions well




##############################################################################################################################################################################

# Extract coefficients and confidence intervals

library(broom)
library(dplyr)
library(readr)

# 1. Extract tidy results with exponentiated HRs and confidence intervals
cox_model_timevaryingCR_HF_fully_adjusted <- tidy(
  cox_model_timevaryingCR_HF_primary_outcome_model_2b,
  exponentiate = TRUE,
  conf.int = TRUE
)

# 2. Clean and format the output
cox_model_timevaryingCR_HF_fully_adjusted_clean <- cox_model_timevaryingCR_HF_fully_adjusted %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))


# Save the rounded table as a CSV file for extraction
output_file_cox_model_timevaryingCR_HF_fully_adjusted_results <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/HF Rehosp/Time varying cox model/cox_model_timevaryingCR_HF_fully_adjusted_clean.csv"
write.csv(cox_model_timevaryingCR_HF_fully_adjusted_clean, file = output_file_cox_model_timevaryingCR_HF_fully_adjusted_results, row.names = FALSE)




# NOTES INTERPRETATION OF MAIN FINDING: Adjusting for all key covariates, CR exposure is not associated with fewer all-cause rehosps compared with no CR exposure







#########################################################################################################################################
#####                            SUMMARY OF HF REHOSP VOLUMES - HF rehosps                                                ###############
#########################################################################################################################################



##NOTE: this includes the filtered out patients from above - reset the data if you do not want the filtered baseline cohort. 


##raw volumes of HF hospitalizations of entire cohort since TAVI + 180 days until study end
## hospitilisations are only HF rehosps post 180days of TAVI discharge

HF_rehosp_raw_vols <- clean_cohort |>
  mutate(tavi_cips_disdate_plus180_days = tavi_cips_disdate + days(180),  
         rehosp_indicator_HF = ifelse((out_readm_hf_date >= tavi_cips_disdate_plus180_days), 1,  0) 
  )
HF_rehosp_raw_vols$rehosp_indicator_HF[is.na(HF_rehosp_raw_vols$rehosp_indicator_HF)] <- 0

HF_rehosp_raw_vols_summary <- HF_rehosp_raw_vols |>
  group_by(rehab_indicator) |>    # Group by rehab exposure
  summarise(
    total_patients = round(n_distinct(person_id) / 5) * 5,      # Unique patients in each group and Round to nearest 5
    total_HF_rehosp = round(sum(rehosp_indicator_HF) / 5) * 5,  # Total HF rehospitalizations and # Round to nearest 5
    avg_HF_rehosp_per_patient = total_HF_rehosp / total_patients, # Avg HF rehosp per patient
  )
print(HF_rehosp_raw_vols_summary)  





# Save the rounded table as a CSV file for extraction
output_file_HF_rehosp_raw_vols_summary <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/HF Rehosp/HF_rehosp_raw_vols_summary.csv"
write.csv(HF_rehosp_raw_vols_summary, file = output_file_HF_rehosp_raw_vols_summary, row.names = FALSE)

##note this excludes the volume of patients (n=2183) who had a rehosp within 180 days post TAVI, that's why total patients coloumn is not the full cohort






#########################################################################################################################################



