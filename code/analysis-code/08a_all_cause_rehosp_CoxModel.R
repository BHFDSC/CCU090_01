#########################################################################################################################################
#####                            OBJECTIVE 4. Effectiveness of cardiac rehabilitation post TAVI                           ###############
#####                                 secondary outcome - all cause readmissions, cox model                               ###############
#########################################################################################################################################


##Cox model with cardiac rehabilitation (CR) as independent variable  using calendar time as the timescale, 
##The observation period is from TAVI discharge + 180days to 31 March 2024 - adjusting for calendar time 
#The model will evaluate how the patient's hazard (instantaneous rate) for all cause rehospitalization based on their participation in rehab

#####################################################
##Step 1: Define and create Calendar Time Variables##
#####################################################

# Load necessary libraries
library(survival)
library(dplyr)
library(lubridate)

#Create Relevant Variables for:
# Time variable: The time to all cause rehospitalization or censoring. all_cause_rehosp_censor_time: days since TAVI discharge+180 days to all cause rehospitalization or censoring (death or end of study period)
# Event variable: A binary variable (1 = all cause rehosp, 0 = censored).
# Exposure variable: Rehab status (1 = rehab, 0 = no rehab).
# Covariates: Adjust for other variables (e.g., age, comorbidities).


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



# out_readm_all_cause_date - study_start: Days from TAVI discharge_plus180 until all cause rehospitalization.
#date_of_death - study_start: Days from TAVI discharge_plus180 until death.
#study_end - study_start: Days from TAVI discharge_plus180 until the end of the study.
#pmin() takes the minimum of these times as the censoring point, considering the event that happens earliest (i.e., all cause rehospitalization, death, or end of follow-up).


##check if any patients have rehosp or censoring time less than study start 
count_invalid <- sum(clean_cohort$all_cause_rehosp_censor_time <= clean_cohort$days_since_tavi_dis_plus180, na.rm = TRUE)
print(count_invalid)
##NOTES, n=8915 (rounded to nearest 5)  patients have have all cause rehosp or censoring time less than or equal to study start (which is tavi dis date plus 180days), only interested in post discharge from TAVIplus 180 so filter out these patients 

# Filter out patients based on rehosp_censor_time conditions
clean_cohort <- clean_cohort %>%
  filter(!(all_cause_rehosp_censor_time <= days_since_tavi_dis_plus180))

count_invalid <- sum(clean_cohort$all_cause_rehosp_censor_time <= clean_cohort$days_since_tavi_dis_plus180, na.rm = TRUE)
print(count_invalid)





#####################################################################
## Fit Cox Proportional Hazards Model ############
#####################################################################
#Fit the Cox proportional hazards model with the coxph function

# Cox model with calendar time as the timescale
cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted <- coxph(
  Surv(all_cause_rehosp_censor_time, all_cause_rehosp_indicator) ~ rehab_indicator,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted)

#notes: unadjusted cox model, starting from 180 days post tavi discharge, shows that CR is significantly associated with a slower rate of all_cause rehosp, however, this is without adjusting 


##########################
### Check Assumptions 
##########################

#Check proportional hazards assumption

# tests whether each covariate in the Cox model satisfies the proportional hazards assumption
#Tests if the effect of each covariate changes over time (which violates the proportional hazards assumption).
# p < 0.05 ??? The PH assumption is violated (the covariate's effect changes over time).
# p ??? 0.05 ??? No evidence of PH violation

# Test proportional hazards assumption
cox.zph(cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted)


#NOTES: P-values are greater than 0.05 and Global is greater than 0.05 therefore proportional hazards assumptions hold true



#extract results
library(broom)
library(dplyr)
library(readr)

# 1. Tidy model results for extraction
tidy_results_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted <- tidy(cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# 2. Extract number of events from the model summary and round
num_events_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted <- summary(cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted)$nevent
rounded_events_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted <- 5 * round(num_events_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted / 5)

# 3. Add as a separate row
events_row_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted <- tibble(
  term = "Number of events (rounded to nearest 5)",
  HR = NA, std.error = NA, statistic = NA,
  p_value = NA, lower_CI = NA, upper_CI = NA,
  extra_info = as.character(rounded_events_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted)
)

# 4. Add column for extra info and combine
tidy_results_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted <- tidy_results_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted %>%
  mutate(extra_info = NA_character_) %>%
  select(term, HR, std.error, statistic, p_value, lower_CI, upper_CI, extra_info)

# 5. Combine and export
output_df_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted <- bind_rows(events_row_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted, tidy_results_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted)


# Save the rounded table as a CSV file for extraction
output_file_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_rehosp/output_df_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted.csv"
write.csv(output_df_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted, file = output_file_cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted, row.names = FALSE)



#######

# Cox model with calendar time as the timescale - adjusting for age and sex
cox_model_CR_all_cause_adjust_age_sex <- coxph(
  Surv(all_cause_rehosp_censor_time, all_cause_rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_CR_all_cause_adjust_age_sex)

#NOTES: adjusting for age and sex, CR exposure is associated with a significantly slower rate of all_cause rehosp vs not exposed


# Test proportional hazards assumption
cox.zph(cox_model_CR_all_cause_adjust_age_sex)


#NOTES: P-values are less than 0.05 and Global is less than 0.05 therefore proportional hazards assumptions does not hold true for age and sex
#Not exptarcting these results b/c proportionla hazards dont hold true



#############################################



# Plot survival curve of unadjusted modle:
plot(survfit(cox_model_CR_all_cause_rehosp_tavi_plus180_unadjusted), main = "Unadjusted Crude Survival Curve for all cause Rehospitalization",
     xlab = "Days since 180 post TAVI discharge to end of followup (31 March 2024)",
     ylab = "Survival Probability")



# Create survival object for individuals in the dataset


surv_object_all_cause_cox_unadjusted <- survfit(
  Surv(all_cause_rehosp_censor_time, all_cause_rehosp_indicator) ~ rehab_indicator, 
  data = clean_cohort)


summary(surv_object_all_cause_cox_unadjusted)

# Plot survival curve for all cause rehospitalization for 2 groups (rehab vs no rehab)
#NOTES THIS IS UNADJUSTED SURVIVAL CURVE
plot(
  surv_object_all_cause_cox_unadjusted,
  col = c("red", "blue"),  # Colors for the curves
  lwd = 2,                 # Line width
  lty = 1,                 # Line type (solid line)
  mark.time = TRUE,        # Add tick marks at censoring times
  main = "Survival Curve for All Cause Rehospitalisation (Unadjusted)",
  xlab = "Starting from days since TAVI discharge plus 180 days and ending at end of followup (31 March 2024)",
  ylab = "Survival Probability",
  conf.int = FALSE         # Disable confidence intervals
)

# Add a legend
legend(
  "bottomleft",
  legend = c("No Rehab", "Rehab"),
  col = c("red", "blue"),
  lwd = 2,
  bty = "n"  # No border around legend
)



summary(surv_object_all_cause_cox_unadjusted, times = c(0, 5, 7, 10, 20, 30, 60, 90, 180, 1000, 1500, 2000))




#PLot survival curve
library(survminer)
library(survival)
library(dplyr)
library(ggplot2)
library(gridExtra)

surv_plot_all_cause_cox_unadjusted <- ggsurvplot(
  surv_object_all_cause_cox_unadjusted,
  data = clean_cohort,
  conf.int = TRUE,  # Show confidence intervals
  pval = FALSE,      
  risk.table = TRUE, 
  legend.title = "Cardiac Rehabilitation Exposure",
  legend.labs = c("Not Exposed to Rehab", "Exposed to Rehab"),
  palette = c("#E69F00", "#0072B2"),  # Distinct colors for rehab vs non-rehab
  linetype = c("solid", "dashed"),  # Solid for rehab, dashed for non-rehab
  break.time.by = 180,  # Adjust x-axis breakpoints for readability
  title = "Unadjusted Survival Curve for All Cause Rehospitalisation",
  xlab = "Calendar Time, starting from 180 days post TAVI discharge to end of follow-up (31 March 2024)",
  ylab = "Survival Probability",
  ggtheme = theme_minimal(base_size = 16)
)

surv_plot_all_cause_cox_unadjusted$plot <- surv_plot_all_cause_cox_unadjusted$plot +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13)
  )

# Print the plot 
print(surv_plot_all_cause_cox_unadjusted)





# Extract and round all relevant columns to nearest 5
rounded_risk_table_long_all_cause_cox_unadjusted <- surv_plot_all_cause_cox_unadjusted$table$data %>%
  mutate(across(
    c(n.risk, n.event, n.censor, cum.n.event, cum.n.censor,strata_size,llabels),
    ~ 5 * round(. / 5)
  ))


# Save the rounded table as a CSV file for extraction
output_file_risk_table_long_all_cause_cox_unadjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_rehosp/rounded_risk_table_long_all_cause_cox_unadjusted.csv"
write.csv(rounded_risk_table_long_all_cause_cox_unadjusted, file = output_file_risk_table_long_all_cause_cox_unadjusted, row.names = FALSE)






##########################################################################################################################
##Adjusted model demographics only  ####################################################
##########################################################################################################################

# Cox model with calendar time as the timescale - adjusting for demographics: age, sex, region, ethnicity, SES, smoking status - to see if minor demographic adjustments change model
cox_model_CR_all_cause_adjust_demographics <- coxph(
  Surv(all_cause_rehosp_censor_time, all_cause_rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + region_name + SMOKING_STATUS,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_CR_all_cause_adjust_demographics)

#NOTES: adjusting for demographics, CR exposure is associated with a significantly slower rate of all_cause rehosp vs not exposed


# Test proportional hazards assumption
cox.zph(cox_model_CR_all_cause_adjust_demographics)

#NOTES: P-values are greater than 0.05 except for age and Global is greater than 0.05. Need to adjust age b/c proportional hazards assumptions does not hold true for age 



#######
#####Stratify by age groups to account for proportonal hazrdas 

# Cox model with calendar time as the timescale - adjusting for demographics: age, sex, region, ethnicity, SES, smoking status - to see if minor demographic adjustments change model
cox_model_CR_all_cause_adjust_demographics <- coxph(
  Surv(all_cause_rehosp_censor_time, all_cause_rehosp_indicator) ~ rehab_indicator + strata(age_groups) + sex_code + ethnicity_5_group + IMD_2019_QUINTILES +region_name + SMOKING_STATUS,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_CR_all_cause_adjust_demographics)

#NOTES: adjusting for demographics, CR exposure is associated with a significantly slower rate of all_cause rehosp vs not exposed


# Test proportional hazards assumption
cox.zph(cox_model_CR_all_cause_adjust_demographics)


#NOTES: P-values are greater than 0.05 and Global is greater than 0.05.Therefore proportional hazards assumptions hold true 





#extract results
library(broom)
library(dplyr)
library(readr)

# Tidy model output (with exponentiated HRs and CIs)
tidy_results_cox_model_CR_all_cause_rehosps_adjusted_demographics <- tidy(cox_model_CR_all_cause_adjust_demographics, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Get number of events and round to nearest 5
num_events_cox_model_CR_all_cause_rehosps_adjusted_demographics <- summary(cox_model_CR_all_cause_adjust_demographics)$nevent
rounded_events_cox_model_CR_all_cause_rehosps_adjusted_demographics <- 5 * round(num_events_cox_model_CR_all_cause_rehosps_adjusted_demographics / 5)

# Add number of events as a top row
events_row_cox_model_CR_all_cause_rehosps_adjusted_demographics <- tibble(
  term = "Number of events (rounded to nearest 5)",
  HR = NA, std.error = NA, statistic = NA,
  p_value = NA, lower_CI = NA, upper_CI = NA,
  extra_info = as.character(rounded_events_cox_model_CR_all_cause_rehosps_adjusted_demographics)
)

# Add an extra_info column to match structure
tidy_results_cox_model_CR_all_cause_rehosps_adjusted_demographics <- tidy_results_cox_model_CR_all_cause_rehosps_adjusted_demographics %>%
  mutate(extra_info = NA_character_) %>%
  select(term, HR, std.error, statistic, p_value, lower_CI, upper_CI, extra_info)

# Combine and export
output_df_cox_model_CR_all_cause_rehosps_adjusted_demographics <- bind_rows(events_row_cox_model_CR_all_cause_rehosps_adjusted_demographics, tidy_results_cox_model_CR_all_cause_rehosps_adjusted_demographics)

# Save the rounded table as a CSV file for extraction
output_file_cox_model_CR_all_cause_rehosps_adjusted_demographics <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_rehosp/all_cause_rehosp_cox model starting from 180 days post tavi/output_df_cox_model_CR_all_cause_rehosps_adjusted_demographics.csv"
write.csv(output_df_cox_model_CR_all_cause_rehosps_adjusted_demographics, file = output_file_cox_model_CR_all_cause_rehosps_adjusted_demographics, row.names = FALSE)














##########################################################################################################################
##Adjusted model 2a: All Covariates  ####################################################
##########################################################################################################################

#Covariate selection by study cardiologist Professor Tom Marwick
#   Complications after TAVI likely dependent on: age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
#  BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
#  previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications  

#Note: did not include PA_SYSTOLIC_PRESSURE_MMHG b/c many missing observations, and was not a key varaible in the clinical selection process 

#Also need to account for calendar time at baseline "year of baseline", given that the start time is TAVI dischargePLus180, we need to account for calendar time differences due to COVID
# include calendar time at baseline in the Cox model, need to extract the year of TAVI discharge and include it as a covariate

#Extract Year of TAVI Discharge
library(lubridate)  # Ensure lubridate is loaded for date handling

# Extract the year from TAVI discharge date
clean_cohort <- clean_cohort %>%
  mutate(year_tavi_discharge = year(tavi_cips_disdate))  # Extract year


#convert to a factor to compare different years independently
#Cox model will estimate a separate hazard ratio for each year relative to the reference year (2018).
clean_cohort <- clean_cohort %>%
  mutate(year_tavi_discharge = as.factor(year_tavi_discharge))


#############################################################
# Cox Proportional Hazards Model 2a ########
#############################################################
#Fit the Cox proportional hazards model with the coxph function


#  Fit Cox model from tavi discahrge plus 180 days, with all covariates
cox_model_taviplus180_all_cause_adjusted_2a <- coxph(
  Surv(all_cause_rehosp_censor_time, all_cause_rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
    BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications + year_tavi_discharge  # Include year
  ,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_taviplus180_all_cause_adjusted_2a)



# Check Assumptions
# Test proportional hazards assumption
cox.zph(cox_model_taviplus180_all_cause_adjusted_2a)


###NOTES: Model runs with all covariates. However, proportional hazards assumptions test does NOT hold
## age_at_tavi and COMORBIDITY_pulmonary_liver_neuro_excardiacvascu is less than 0.05
##therefore proportional hazards assumptions do not hold for these. 

##Stratify by COMORBIDITY_pulmonary_liver_neuro_excardiacvascu and use age_groups categories




#######################################################################################################################################################################################################
# Fit Cox Proportional Hazards Model 2b
#######################################################################################################################################################################################################
#Fit the Cox proportional hazards model with the coxph function
####Stratify by COMORBIDITY_pulmonary_liver_neuro_excardiacvascu and use age_groups categories


#  Fit Cox model from tavi discahrge plus 180 days, with all covariates
cox_model_taviplus180_all_cause_adjusted_2b <- coxph(
  Surv(all_cause_rehosp_censor_time, all_cause_rehosp_indicator) ~ rehab_indicator + age_groups + #use age_groups categorical (refrence group is <70years old)
    sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
    BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + 
    strata(COMORBIDITY_pulmonary_liver_neuro_excardiacvascu) + #stratify hazard violoting variable
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications + year_tavi_discharge  # Include year
  ,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_taviplus180_all_cause_adjusted_2b)



# Check Assumptions
# Test proportional hazards assumption
cox.zph(cox_model_taviplus180_all_cause_adjusted_2b)


###NOTES: Proportional hazard assumption test still does not holds true for age_groups
##interpretation: Therefore also stratify age_groups



#######################################################################################################################################################################################################
# Fit Cox Proportional Hazards Model 2c
#######################################################################################################################################################################################################
#Fit the Cox proportional hazards model with the coxph function
### stratify age_groups and COMORBIDITY


#  Fit Cox model from tavi discahrge plus 180 days, with all covariates, including region_name (REGION was added on)
cox_model_taviplus180_all_cause_adjusted_2c <- coxph(
  Surv(all_cause_rehosp_censor_time, all_cause_rehosp_indicator) ~ rehab_indicator +
    sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + region_name +
    BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + 
    strata(COMORBIDITY_pulmonary_liver_neuro_excardiacvascu, age_groups) + #stratify hazard violating variables age groups and comorbidity
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications + year_tavi_discharge  # Include year
  ,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_taviplus180_all_cause_adjusted_2c)



# Check Assumptions
# Test proportional hazards assumption
cox.zph(cox_model_taviplus180_all_cause_adjusted_2c)



###NOTES: Proportional hazard assumption test holds true for Global and holds true for all other variables 

#NOTES: P-values are greater than 0.05 for Global, hold true for all other variables


# NOTES INTERPRETATION OF MAIN FINDING: Adjusting for all key covariates, accounting for baseline calendar time CR, starting from days since tavi discharge plus 180days, exposure is associated with slower  rate of all cause rehosp 


#######################################################################################################################################################################################################
# Extract model summary adjusted Model 2c
#######################################################################################################################################################################################################



#Broom package extraction 
library(broom)
library(dplyr)
library(readr)

# Tidy the Cox model with HR and CIs
results_cox_model_taviplus180_all_cause_rehosp_fully_adjusted_2c <- tidy(
  cox_model_taviplus180_all_cause_adjusted_2c,
  exponentiate = TRUE,  # gives HR instead of log(HR)
  conf.int = TRUE
)

# Optional: clean and round columns
results_cox_model_taviplus180_all_cause_rehosp_fully_adjusted_2c_clean <- results_cox_model_taviplus180_all_cause_rehosp_fully_adjusted_2c %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))


# Save the rounded table as a CSV file for extraction
output_file_results_cox_model_taviplus180_all_cause_rehosp_fully_adjusted_2c_clean <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_rehosp/all_cause_rehosp_cox model starting from 180 days post tavi/results_cox_model_taviplus180_all_cause_rehosp_fully_adjusted_2c_clean.csv"
write.csv(results_cox_model_taviplus180_all_cause_rehosp_fully_adjusted_2c_clean, file = output_file_results_cox_model_taviplus180_all_cause_rehosp_fully_adjusted_2c_clean, row.names = FALSE)


#######################################################################################################################################################################################################


#########################################################################################################################################
#####                            SUMMARY OF REHOSP VOLUMES - all cause rehosps                                     ###############
#########################################################################################################################################


##NOTE: this is based on the filtered out patients from above (n=8915 - rounded to nearest 5), these patients had rehop within 180 days of TAVI - reset the data if you do not want the filtered baseline cohort. 


##raw volumes of all cause hospitalizations of entire cohort since TAVI + 180 days until study end
##hospitilisations are only rehosps post 180days of TAVI discharge

all_cause_rehosp_raw_vols <- clean_cohort |>
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

##note this excludes the volume of patients (n=8915 rounded to nearest 5) who had an all cause rehosp within 180 days post TAVI, that's why total patients coloumn is not the full baseline cohort








#########################################################################################################################################








