#########################################################################################################################################
#####                            OBJECTIVE 4. Effectiveness of cardiac rehabilitation  post TAVI                          ###############
#####      SECONDARY outcome - all cause Mortality, cox model with CR time varying exposure                                ###############
#########################################################################################################################################

#a.	all-cause mortality -  cox model with time varying exposure - take out death from censoring time and change outcome to all-cause mortality


##Cox model with time-varying exposure to cardiac rehabilitation (CR) using calendar time as the timescale, 
##The observation period is from 1 January 2018 to 31 March 2024
#The model will evaluate how the patient's hazard for ALL CAUSE mortality changes based on their participation in rehab during different time intervals

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
# censor_time: days since 1 January 2018 to censoring (end of study period)

# Define study start and end dates
study_start <- as.Date("2018-01-01")
study_end <- as.Date("2024-03-31")

# Add necessary time variables
clean_cohort_all_cause_mortality <- clean_cohort %>%
  mutate(
    # Days since 1 January 2018 (calendar time)
    days_since_2018 = as.numeric(tavi_cips_disdate - study_start),  #this will be used in step2 (as time_1 in the split data with time intervals)
    
    # TAVI time (days from 1 January 2018 to TAVI procedure discharge date)
    tavi_time = as.numeric(tavi_cips_disdate - study_start),
    
    # Rehab time (days from 1 January 2018 to the start of the 1st rehab appointment)
    rehab_time = ifelse(!is.na(exp_cardiac_rehab_appt_date_6m), 
                        as.numeric(exp_cardiac_rehab_appt_date_6m - study_start), NA),
    
    # # Define time to event: either death or end of study. Time of all cause mortality or censoring (death or end of study period)
    mortality_censor_time_all_cause = pmin(as.numeric(date_of_death - study_start),                # Days to death
                                        as.numeric(study_end - study_start), na.rm = TRUE),     # Days to study end
    
    # Rehab exposure indicator
    rehab_indicator = ifelse(is.na(exp_cardiac_rehab_6m) | exp_cardiac_rehab_6m == 0, 0, 1),
    
    # Define mortality indicator (1 for death, 0 for censoring)
    mortality_indicator_all_cause = ifelse(!is.na(date_of_death), 1, 0)
    
  )


##check if any patients die or censoring time less than or equal to TAVI time and die within 180 days of TAVI discharge
count_invalid_all_cause_mortality <- sum(clean_cohort_all_cause_mortality$mortality_censor_time_all_cause <= clean_cohort_all_cause_mortality$tavi_time + 180, na.rm = TRUE)
                                        
                                      
print(count_invalid_all_cause_mortality)


##NOTES, n=1295 (rounded to nearest 5)  patients die within 180 days of tavi discharge.
##only interested in post discharge from TAVI plus 180 days to account for enough time post discharge to enrol in rehab,  so filter out these patients within 180 days post discharge




# Grouping by rehab_indicator and summarizing count
count_invalid_all_cause_mortality <- clean_cohort_all_cause_mortality %>%
  group_by(rehab_indicator) %>%  # Group by rehab exposure
  summarise(
    count_invalid = sum(
      mortality_censor_time_all_cause <= tavi_time | 
        mortality_censor_time_all_cause <= tavi_time + 180, na.rm = TRUE
    ) 
  )

# Print results
print(count_invalid_all_cause_mortality)


### 1250 (rounded to nearest 5) patients in the non CR group and 40 (rounded to nearest 5) patients in the CR group have all cause mortality within 180 days of TAVI discharge - therefore filter them out






# Filter out patients based on mortalty_censor_time conditions, filter any deaths within 180 days post discharge
clean_cohort_all_cause_mortality <- clean_cohort_all_cause_mortality %>%
  filter(!(mortality_censor_time_all_cause <= tavi_time | mortality_censor_time_all_cause <= days_since_2018 | mortality_censor_time_all_cause <= tavi_time + 180))

# Count the number of invalid cases 
count_invalid_all_cause_mortality <- sum(clean_cohort_all_cause_mortality$mortality_censor_time_all_cause <= clean_cohort_all_cause_mortality$tavi_time,
                                         clean_cohort_all_cause_mortality$mortality_censor_time_all_cause <=  clean_cohort_all_cause_mortality$ tavi_time + 180, na.rm = TRUE)
print(count_invalid_all_cause_mortality)


###NOTES: FILTERED OUT 1295 (rounded to nearest 5) PATIENTS who died at a time less than or equal to TAVI time and or patients who died within 180 days of tavi discharge






#####################################################################
##Step 2: Split Data into Time Intervals for Time-Varying Cox Model##
#####################################################################

#To model time-varying exposure, split data based on whether the patient was in CR (rehab_time):
#Pre-Rehab: from TAVI to the start of CR (if applicable).
#Post-Rehab: from start of CR to all cause mortalty or end of study.


split_data_all_cause_mortality <- clean_cohort_all_cause_mortality |>                   #Part 1: Splitting the Dataset into Periods Based on Rehab Time
  mutate(
    period_start_1 = days_since_2018,           # The start of the first period is the start of the follow-up time.
    period_end_1 = case_when(
      is.na(rehab_time) ~ mortality_censor_time_all_cause,   # If the rehab time is NA, the period ends at death or study end
      .default = rehab_time                     # Otherwise, the period ends at rehab start time.
    ),
    period_start_2 = case_when(
      !is.na(rehab_time) ~ rehab_time,          # The second period starts at rehab time if rehab time is not NA.
      .default = NA                             # Otherwise, no second period.
    ),
    period_end_2 = case_when(
      period_end_1 == mortality_censor_time_all_cause ~ NA,  # If the first period already ends at censoring, no second period.
      .default = mortality_censor_time_all_cause             # Otherwise, the second period ends at death or censoring.
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

split_data2_all_cause_mortality <- split_data_all_cause_mortality |> left_join(clean_cohort_all_cause_mortality |> dplyr::select(person_id, days_since_2018, tavi_time, rehab_time, mortality_censor_time_all_cause, rehab_indicator,
                                                                                                                        mortality_indicator_all_cause, age_at_tavi, age_groups, sex_code, ethnicity_5_group, IMD_2019_QUINTILES, region_name, 
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
  mutate(mortality_indicator_all_cause = case_when(period == "1" & dupes == 1 ~ 0, .default = mortality_indicator_all_cause)) |>         # Correct the mortality indicator
  dplyr::select(-dupes)








#####################################################################
##Step 3: Fit Time-Varying Cox Proportional Hazards Model############
#####################################################################
#Fit the Cox proportional hazards model with the coxph function

# Fit time-varying Cox model with calendar time as the timescale
cox_model_timevaryingCR_all_cause_mortality_unadjusted <- coxph(
  Surv(period_start, period_end, mortality_indicator_all_cause) ~ rehab_indicator,
  data = split_data2_all_cause_mortality
)

# Summary of the Cox model
summary(cox_model_timevaryingCR_all_cause_mortality_unadjusted)



#count number of patients who had rehab after the death, therefore filter out these patients - given not posisble and likley error in data
count_invalid_all_cause_mortality2 <- sum(split_data2_all_cause_mortality$period_end <= split_data2_all_cause_mortality$period_start, na.rm = TRUE)
print(count_invalid_all_cause_mortality2)


### NOTES: n= 1 patient had rehab after death, therefore filter out this patients
split_data2_all_cause_mortality <- split_data2_all_cause_mortality %>%
  filter(!(period_end <= period_start))

sum(split_data2_all_cause_mortality$period_end <= split_data2_all_cause_mortality$period_start, na.rm = TRUE)


# Fit time-varying Cox model with calendar time as the timescale, with 1 patients filtered out
cox_model_timevaryingCR_all_cause_mortality_unadjusted <- coxph(
  Surv(period_start, period_end, mortality_indicator_all_cause) ~ rehab_indicator,
  data = split_data2_all_cause_mortality
)

# Summary of the Cox model
summary(cox_model_timevaryingCR_all_cause_mortality_unadjusted)

#this worked with no warnings




#####################################################################
##Step 4: Check Assumptions  ################################
#####################################################################

# Test proportional hazards assumption
cox.zph(cox_model_timevaryingCR_all_cause_mortality_unadjusted)

# Age  p-value is greater than 0.05 therefore hazard assumptions holds true for unadjusted model. Need to run model adjusting for all other covarites)





#extract results
library(broom)
library(dplyr)
library(readr)

# 1. Tidy model results for extraction
tidy_results_cox_model_timevaryingCR_all_cause_mortality_unadjusted <- tidy(cox_model_timevaryingCR_all_cause_mortality_unadjusted, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# 2. Extract number of events from the model summary and round
num_events_cox_model_timevaryingCR_all_cause_mortality_unadjusted <- summary(cox_model_timevaryingCR_all_cause_mortality_unadjusted)$nevent
rounded_events_cox_model_timevaryingCR_all_cause_mortality_unadjusted <- 5 * round(num_events_cox_model_timevaryingCR_all_cause_mortality_unadjusted / 5)

# 3. Add as a separate row
events_row_cox_model_timevaryingCR_all_cause_mortality_unadjusted <- tibble(
  term = "Number of events (rounded to nearest 5)",
  HR = NA, std.error = NA, statistic = NA,
  p_value = NA, lower_CI = NA, upper_CI = NA,
  extra_info = as.character(rounded_events_cox_model_timevaryingCR_all_cause_mortality_unadjusted)
)

# 4. Add column for extra info and combine
tidy_results_cox_model_timevaryingCR_all_cause_mortality_unadjusted <- tidy_results_cox_model_timevaryingCR_all_cause_mortality_unadjusted %>%
  mutate(extra_info = NA_character_) %>%
  select(term, HR, std.error, statistic, p_value, lower_CI, upper_CI, extra_info)

# 5. Combine and export
output_HR_cox_model_timevaryingCR_all_cause_mortality_unadjusted <- bind_rows(events_row_cox_model_timevaryingCR_all_cause_mortality_unadjusted,
                                                                                tidy_results_cox_model_timevaryingCR_all_cause_mortality_unadjusted)


# Save the rounded table as a CSV file for extraction
output_file_HR_cox_model_timevaryingCR_all_cause_mortality_unadjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_mortality/mortality_time_varying_CR_cox_model/output_HR_cox_model_timevaryingCR_all_cause_mortality_unadjusted.csv"
write.csv(output_HR_cox_model_timevaryingCR_all_cause_mortality_unadjusted, file = output_file_HR_cox_model_timevaryingCR_all_cause_mortality_unadjusted, row.names = FALSE)









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
cox_model_timevaryingCR_all_cause_mortality_3a <- coxph(
  Surv(period_start, period_end, mortality_indicator_all_cause) ~ rehab_indicator + age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
    BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications,
  data = split_data2_all_cause_mortality
)

# Summary of the Cox model
summary(cox_model_timevaryingCR_all_cause_mortality_3a)


#Check proportional hazards assumption
# Test proportional hazards assumption
cox.zph(cox_model_timevaryingCR_all_cause_mortality_3a)


#NOTES: model runs well. Proportional hazards assumption test shows that age_at_tavi; NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE; binary_Known_DIABETES ; LV_FUNCTION have p-values less than 0.05, 
##therefore hazard assumption test does not hold. 
##Need to stratify these variables That Violate PH For age_at_tavi; NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE; binary_Known_DIABETES ; LV_FUNCTION 
##stratification removes the need to estimate their hazard ratio over time.




########################################################################################################################################################################
##Adjusted model 3b: Stratify age_groups; NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE; binary_Known_DIABETES ; LV_FUNCTION, to correct for hazrad assumptions in model above
#######################################################################################################################################################################

#Fit the Cox proportional hazards model with the coxph function

#Fit time-varying Cox model with calendar time as the timescale, with all covariates
cox_model_timevaryingCR_all_cause_mortality_3b <- coxph(
  Surv(period_start, period_end, mortality_indicator_all_cause) ~ rehab_indicator + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
    BMI_category + KATZ_INDEX_CAT +  binary_known_poor_mobility + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + MITRAL_REGURGITATION + TAVI_procedural_complications +
  strata(age_groups, NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE, binary_Known_DIABETES, LV_FUNCTION),  # Stratifying PH-violating factors, 
  data = split_data2_all_cause_mortality
)

# Summary of the Cox model
summary(cox_model_timevaryingCR_all_cause_mortality_3b)



# Test proportional hazards assumption
cox.zph(cox_model_timevaryingCR_all_cause_mortality_3b)

## hazard assumption test holds true for all variables and global
##Interpretation: Stratifying PH-violating factors, including age groups, and after adjusting for all other covariates, CR exposure is not associated with a higher rate of all-cause mortality. 

##Stratification removes the need to estimate an explicit HR for age while accounting for different baseline hazards across age groups
##Removes the need to model age explicitly, allowing different baseline hazards for different ages.
## Useful when age has a strong effect on survival that changes over time.






#extract results and save in home folder as csv
library(broom)
library(dplyr)
library(readr)

# 1. Tidy model results for extraction
tidy_results_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted <- tidy(cox_model_timevaryingCR_all_cause_mortality_3b, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# 2. Extract number of events from the model summary and round
num_events_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted <- summary(cox_model_timevaryingCR_all_cause_mortality_3b)$nevent
rounded_events_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted <- 5 * round(num_events_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted / 5)

# 3. Add as a separate row
events_row_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted <- tibble(
  term = "Number of events (rounded to nearest 5)",
  HR = NA, std.error = NA, statistic = NA,
  p_value = NA, lower_CI = NA, upper_CI = NA,
  extra_info = as.character(rounded_events_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted)
)

# 4. Add column for extra info and combine
tidy_results_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted <- tidy_results_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted %>%
  mutate(extra_info = NA_character_) %>%
  select(term, HR, std.error, statistic, p_value, lower_CI, upper_CI, extra_info)

# 5. Combine and export
output_HR_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted <- bind_rows(events_row_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted,
                                                                              tidy_results_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted)


# Save the rounded table as a CSV file for extraction
output_file_HR_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_mortality/mortality_time_varying_CR_cox_model/output_HR_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted.csv"
write.csv(output_HR_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted, file = output_file_HR_cox_model_timevaryingCR_all_cause_mortality_fully_adjusted, row.names = FALSE)








#########################################################################################################################################
#####                            SUMMARY OF all cause mortality VOLUMES                                      ###############
#########################################################################################################################################


##NOTE: this is based off the filtered out patients from above (n=1295 (rounded to nearest 5)  patients die within 180 days of tavi discharge) - reset the data if you do not want the filtered baseline cohort. 


##raw volumes of all cause mortality of entire cohort since TAVI + 180 days until study end
##mortality are only post 180days of TAVI discharge

all_cause_mort_raw_vols <- clean_cohort_all_cause_mortality |>
  mutate(tavi_cips_disdate_plus180_days = tavi_cips_disdate + days(180),  
         mortality_indicator_all_cause = ifelse((date_of_death >= tavi_cips_disdate_plus180_days), 1,  0) 
  )
all_cause_mort_raw_vols$mortality_indicator_all_cause[is.na(all_cause_mort_raw_vols$mortality_indicator_all_cause)] <- 0

all_cause_mortality_raw_vols_summary <- all_cause_mort_raw_vols |>
  group_by(rehab_indicator) |>    # Group by rehab exposure
  summarise(
    total_patients = round(n_distinct(person_id) / 5) * 5,      # Unique patients in each group and Round to nearest 5
    total_all_cause_mortality = round(sum(mortality_indicator_all_cause) / 5) * 5,  # Total all-cause rehospitalizations and # Round to nearest 5
    avg_all_cause_mortality_per_patient = total_all_cause_mortality / total_patients, # Avg all-cause rehosp per patient
  )
print(all_cause_mortality_raw_vols_summary)  


# Save the rounded table as a CSV file for extraction
output_file_all_cause_mortality_raw_vols_summary <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_mortality/all_cause_mortality_raw_vols_summary.csv"
write.csv(all_cause_mortality_raw_vols_summary, file = output_file_all_cause_mortality_raw_vols_summary, row.names = FALSE)





#########################################################################################################################################
















