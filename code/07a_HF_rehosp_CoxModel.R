#########################################################################################################################################
#####                            OBJECTIVE 4. Effectiveness of cardiac rehabilitation post TAVI                           ###############
#####                                 Primary outcome - hospital readmission with HF, cox model                           ###############
#########################################################################################################################################


##Cox model with cardiac rehabilitation (CR) as independent variable  using calendar time as the timescale, 
##The observation period is from TAVI discharge + 180days to 31 March 2024 - adjusting for calendar time 
#The model will evaluate how the patient's hazard (instantaneous rate) for rehospitalization for HF changes based on their participation in rehab

#####################################################
##Step 1: Define and create Calendar Time Variables##
#####################################################

# Load necessary libraries
library(survival)
library(dplyr)
library(lubridate)

#Create Relevant Variables for:
# Time variable: The time to rehospitalization or censoring. rehosp_censor_time: days since TAVI discharge+180 days to rehospitalization for HF or censoring (death or end of study period)
# Event variable: A binary variable (1 = HF rehosp, 0 = censored).
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
    hf_rehosp_censor_time = pmin(as.numeric(out_readm_hf_date - study_start_new),         # Days to HF readmission
                              as.numeric(date_of_death - study_start_new),                # Days to death
                              as.numeric(study_end - study_start_new), na.rm = TRUE),     # Days to study end
    
    # Rehab exposure indicator
    rehab_indicator = ifelse(is.na(exp_cardiac_rehab_6m) | exp_cardiac_rehab_6m == 0, 0, 1),
    
    # Rehospitalization indicator (1 for HF rehospitalization, 0 for censoring)
    hf_rehosp_indicator = ifelse(is.na(out_readm_hf_flag) | out_readm_hf_flag == 0, 0, 1)
    
  )



# out_readm_hf_date - study_start: Days from TAVI discharge_plus180 until rehospitalization for HF.
#date_of_death - study_start: Days from TAVI discharge_plus180 until death.
#study_end - study_start: Days from TAVI discharge_plus180 until the end of the study.
#pmin() takes the minimum of these times as the censoring point, considering the event that happens earliest (i.e., HF rehospitalization, death, or end of follow-up).


##check if any patients have rehosp or censoring time less than or equal to study start 
count_invalid <- sum(clean_cohort$hf_rehosp_censor_time <= clean_cohort$days_since_tavi_dis_plus180, na.rm = TRUE)
print(count_invalid)
##NOTES, n= just over 2k patients have have rehosp or censoring time less than or equal to study start (which is tavi dis date plus 180days), only interested in post discharge from TAVIplus 180 so filter out these patients 

# Filter out patients based on rehosp_censor_time conditions
clean_cohort <- clean_cohort %>%
  filter(!(hf_rehosp_censor_time <= days_since_tavi_dis_plus180))

count_invalid <- sum(clean_cohort$hf_rehosp_censor_time <= clean_cohort$days_since_tavi_dis_plus180, na.rm = TRUE)
print(count_invalid)





#####################################################################
## Fit Cox Proportional Hazards Model ############
#####################################################################
#Fit the Cox proportional hazards model with the coxph function

# Cox model with calendar time as the timescale
cox_model_CR_HF_primary_outcome <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_CR_HF_primary_outcome)


#extract results
library(broom)
library(dplyr)
library(readr)

# 1. Tidy model results for extraction
tidy_results_cox_model_CR_HF_primary_outcome <- tidy(cox_model_CR_HF_primary_outcome, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# 2. Extract number of events from the model summary and round
num_events_cox_model_CR_HF_primary_outcome <- summary(cox_model_CR_HF_primary_outcome)$nevent
rounded_events_cox_model_CR_HF_primary_outcome <- 5 * round(num_events_cox_model_CR_HF_primary_outcome / 5)

# 3. Add as a separate row
events_row_cox_model_CR_HF_primary_outcome <- tibble(
  term = "Number of events (rounded to nearest 5)",
  HR = NA, std.error = NA, statistic = NA,
  p_value = NA, lower_CI = NA, upper_CI = NA,
  extra_info = as.character(rounded_events_cox_model_CR_HF_primary_outcome)
)

# 4. Add column for extra info and combine
tidy_results_cox_model_CR_HF_primary_outcome <- tidy_results_cox_model_CR_HF_primary_outcome %>%
  mutate(extra_info = NA_character_) %>%
  select(term, HR, std.error, statistic, p_value, lower_CI, upper_CI, extra_info)

# 5. Combine and export
output_df_cox_model_CR_HF_primary_outcome_unadjusted <- bind_rows(events_row_cox_model_CR_HF_primary_outcome, tidy_results_cox_model_CR_HF_primary_outcome)


# Save the rounded table as a CSV file for extraction
output_file_cox_model_HF_primary_outcome_results_01 <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/HF Rehosp/output_df_cox_model_CR_HF_primary_outcome_unadjusted.csv"
write.csv(output_df_cox_model_CR_HF_primary_outcome_unadjusted, file = output_file_cox_model_HF_primary_outcome_results_01, row.names = FALSE)





# Cox model with calendar time as the timescale - adjusting for demographics: age, sex
cox_model_CR_HF_primary_outcome_age_sex <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_CR_HF_primary_outcome_age_sex)




####
# Cox model with calendar time as the timescale - adjusting for demographics: age, sex, region, ethnicity, SES, smoking status - to see if minor demographic adjustments change model
cox_model_CR_HF_primary_outcome_demographics <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES +region_name + SMOKING_STATUS,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_CR_HF_primary_outcome_demographics)


#Check proportional hazards assumption

# tests whether each covariate in the Cox model satisfies the proportional hazards assumption
#Tests if the effect of each covariate changes over time (which violates the proportional hazards assumption).
# p < 0.05 ??? The PH assumption is violated (the covariate's effect changes over time).
# p ??? 0.05 ??? No evidence of PH violation

cox.zph(cox_model_CR_HF_primary_outcome_demographics)

##Region doesnt meet proportional hazards assumption test so need to stratify by region




####startify by region ####
# Cox model with calendar time as the timescale - adjusting for demographics: age, sex, region, ethnicity, SES, smoking status - to see if minor demographic adjustments change model
cox_model_CR_HF_primary_outcome_demographics <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES +strata(region_name) + SMOKING_STATUS,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_CR_HF_primary_outcome_demographics)


#Check proportional hazards assumption

# tests whether each covariate in the Cox model satisfies the proportional hazards assumption
#Tests if the effect of each covariate changes over time (which violates the proportional hazards assumption).
# p < 0.05 ??? The PH assumption is violated (the covariate's effect changes over time).
# p ??? 0.05 ??? No evidence of PH violation

cox.zph(cox_model_CR_HF_primary_outcome_demographics)

##proportional hazard assumption test is meet for all varaibles and for global
##Interpretation: after minor adjustment for demographic variables, there is no evidnece of an association b/w CR and lower risk of HF rehosp (p=0.456)




#extract results
library(broom)
library(dplyr)
library(readr)

# Tidy model output (with exponentiated HRs and CIs)
tidy_results_cox_model_CR_HF_primary_outcome_demographics <- tidy(cox_model_CR_HF_primary_outcome_demographics, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Get number of events and round to nearest 5
num_events_cox_model_CR_HF_primary_outcome_demographics <- summary(cox_model_CR_HF_primary_outcome_demographics)$nevent
rounded_events_cox_model_CR_HF_primary_outcome_demographics <- 5 * round(num_events_cox_model_CR_HF_primary_outcome_demographics / 5)

# Add number of events as a top row
events_row_cox_model_CR_HF_primary_outcome_demographics <- tibble(
  term = "Number of events (rounded to nearest 5)",
  HR = NA, std.error = NA, statistic = NA,
  p_value = NA, lower_CI = NA, upper_CI = NA,
  extra_info = as.character(rounded_events_cox_model_CR_HF_primary_outcome_demographics)
)

# Add an extra_info column to match structure
tidy_results_cox_model_CR_HF_primary_outcome_demographics <- tidy_results_cox_model_CR_HF_primary_outcome_demographics %>%
  mutate(extra_info = NA_character_) %>%
  select(term, HR, std.error, statistic, p_value, lower_CI, upper_CI, extra_info)

# Combine and export
output_df_cox_model_CR_HF_primary_outcome_demographics_adjusted <- bind_rows(events_row_cox_model_CR_HF_primary_outcome_demographics, tidy_results_cox_model_CR_HF_primary_outcome_demographics)

# Save the rounded table as a CSV file for extraction
output_file_cox_model_CR_HF_primary_outcome_demographics_adjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/HF Rehosp/HF cox model starting from 180 days post tavi/output_df_cox_model_CR_HF_primary_outcome_demographics_adjusted.csv"
write.csv(output_df_cox_model_CR_HF_primary_outcome_demographics_adjusted, file = output_file_cox_model_CR_HF_primary_outcome_demographics_adjusted, row.names = FALSE)








#####################################################################
##       Check Assumptions and plot age and sex ################################
#####################################################################

#Check proportional hazards assumption

# tests whether each covariate in the Cox model satisfies the proportional hazards assumption
#Tests if the effect of each covariate changes over time (which violates the proportional hazards assumption).
# p < 0.05 ??? The PH assumption is violated (the covariate's effect changes over time).
# p ??? 0.05 ??? No evidence of PH violation

# Test proportional hazards assumption
cox.zph(cox_model_CR_HF_primary_outcome)

cox.zph(cox_model_CR_HF_primary_outcome_age_sex)

#NOTES: P-values are greater than 0.05 and Global is greater than 0.05 therefore proportional hazards assumptions hold true




# Plot survival curve:
plot(survfit(cox_model_CR_HF_primary_outcome), main = "Survival Curve for HF Rehospitalization",
     xlab = "Days since 180 post TAVI discharge to end of followup (31 March 2024)",
     ylab = "Survival Probability")



# Create survival object for individuals in the dataset


surv_object_HF_cox <- survfit(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator, 
  data = clean_cohort)


summary(surv_object_HF_cox)

# Plot survival curve for HF rehospitalization for 2 groups (rehab vs no rehab)
#NOTES THIS IS UNADJUSTED SURVIVAL CURVE
plot(
  surv_object_HF_cox,
  col = c("red", "blue"),  # Colors for the curves
  lwd = 2,                 # Line width
  lty = 1,                 # Line type (solid line)
  mark.time = TRUE,        # Add tick marks at censoring times
  main = "Survival Curve for HF Rehospitalization (Unadjusted)",
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



summary(surv_object_HF_cox, times = c(0, 5, 7, 10, 20, 30, 60, 90, 180, 1000, 1500, 2000))




#PLot survival curve
library(survminer)
library(survival)
library(dplyr)
library(ggplot2)
library(gridExtra)

surv_plot_HF_cox <- ggsurvplot(
  surv_object_HF_cox,
  data = clean_cohort,
  conf.int = TRUE,  # Show confidence intervals
  pval = FALSE,      
  risk.table = TRUE, 
  legend.title = "Cardiac Rehabilitation Exposure",
  legend.labs = c("Not Exposed to Rehab", "Exposed to Rehab"),
  palette = c("#E69F00", "#0072B2"),  # Distinct colors for rehab vs non-rehab
  linetype = c("solid", "dashed"),  # Solid for rehab, dashed for non-rehab
  break.time.by = 180,  # Adjust x-axis breakpoints for readability
  title = "Unadjusted Survival Curve for Heart Failure Rehospitalisation",
  xlab = "Calendar Time, starting from 180 days post TAVI discharge to end of follow-up (31 March 2024)",
  ylab = "Survival Probability",
  ggtheme = theme_minimal(base_size = 16)
)

surv_plot_HF_cox$plot <- surv_plot_HF_cox$plot +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13)
  )

# Print the plot 
print(surv_plot_HF_cox)



# Extract and round all relevant columns to nearest 5
risk_table_long_HF_cox_unadjusted <- surv_plot_HF_cox$table$data %>%
  mutate(across(
    c(n.risk, n.event, n.censor, cum.n.event, cum.n.censor,strata_size,llabels),
    ~ 5 * round(. / 5)
  ))


# Save the rounded table as a CSV file for extraction
output_file_risk_table_long_HF_cox_unadjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/HF Rehosp/rounded_risk_table_long_HF_cox_unadjusted.csv"
write.csv(risk_table_long_HF_cox_unadjusted, file = output_file_risk_table_long_HF_cox_unadjusted, row.names = FALSE)






  





##########################################################################################################################
## Plot Adjusted survival curve for Age and Sex #################################################################
##########################################################################################################################


# Define a new dataset for adjusting survival curve, adjusting for age and sex
surv_HF_cox_new_data <- data.frame(
  rehab_indicator = c(0, 1),          # No rehab vs Rehab
  age_at_tavi = mean(clean_cohort$age_at_tavi, na.rm = TRUE),  # Adjust for mean age
  sex_code = factor(c("Female","Male"), levels = levels(clean_cohort$sex_code)),   # Adjust sex
  # Include other covariates if necessary
  stringsAsFactors = FALSE
)

# Fit the survival model and compute predicted survival
surv_object_HF_cox_adj <- survfit(
  cox_model_CR_HF_primary_outcome, 
  newdata = surv_HF_cox_new_data
)

# Plot the adjusted survival curves
plot(
  surv_object_HF_cox_adj,
  col = c("red", "blue"),            # Colors for No Rehab and Rehab
  lwd = 2,                           # Line width
  lty = 1,                           # Line type
  mark.time = TRUE,                  # Add censoring marks
  main = "Survival Curve for HF Rehospitalization (Adjusted for Age and Sex)",
  xlab = "Calendar Time, starting from days since TAVI discharge plus 180 days and ending at end of followup (31 March 2024)", 
  ylab = "Survival Probability",
  conf.int = FALSE                   # Disable confidence intervals
)

# Add a legend
legend(
  "bottomleft",
  legend = c("No Rehab (Adjusted)", "Rehab (Adjusted)"),
  col = c("red", "blue"),
  lwd = 2,
  bty = "n"                          # No border around the legend
)


###Explanation:
# surv_HF_cox_new_data: This dataset ensures predictions are generated for "No Rehab" and "Rehab" groups while holding covariates constant (e.g., mean age, specific sex).
#  survfit(): Uses the fitted Cox model and new_data to compute adjusted survival probabilities.
#  Plot: The adjusted survival curves are plotted, showing the difference between "No Rehab" and "Rehab" while accounting for age and sex.






#PLot survival curve  -  adjusted survival curves for age and sex
library(survminer)

ggsurvplot(
  surv_object_HF_cox_adj,
  data = clean_cohort,
  conf.int = TRUE,  # Show confidence intervals
  pval = FALSE,      
  risk.table = FALSE, 
  legend.title = "Rehabilitation Exposure",
  legend.labs = c("Not Exposed to Rehab", "Exposed to Rehab"),
  palette = c("#E69F00", "#0072B2"),  # Distinct colors for rehab vs non-rehab
  linetype = c("solid", "dashed"),  # Solid for rehab, dashed for non-rehab
  break.time.by = 180,  # Adjust x-axis breakpoints for readability
  title = "Survival Curve for Heart Failure Rehospitalisation (Adjusted for Age and Sex)",
  xlab = "Calendar Time, starting from 180 days post TAVI discharge to end of follow-up (31 March 2024)",
  ylab = "Survival Probability",
  ggtheme = theme_minimal()
)




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


# Fit Cox model with calendar time as the timescale, with all covariates
cox_model_HF_primary_outcome_model_2a <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
  BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
  previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications + year_tavi_discharge  # Include year
  ,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_HF_primary_outcome_model_2a)



# Check Assumptions
# Test proportional hazards assumption
cox.zph(cox_model_HF_primary_outcome_model_2a)


###NOTES: Model runs with all covariates. However, proportional hazards assumptions test does NOT hold
## binary_known_poor_mobility is less than 0.05
##therefore proportional hazards assumptions do not hold for these. 
##Stratify by this variables to correct for hazrad assumptions or time interaction for binary_known_poor_mobility





#######################################################################################################################################################################################################
# Fit Cox Proportional Hazards Model 2b
#######################################################################################################################################################################################################
#Fit the Cox proportional hazards model with the coxph function
##time interaction for binary_known_poor_mobility



clean_cohort <- clean_cohort %>%
  mutate(
    binary_known_poor_mobility_time_interaction  = binary_known_poor_mobility * days_since_tavi_dis_plus180
  )



# Fit Cox model with calendar time as the timescale, with all covariates
cox_model_HF_primary_outcome_model_2b <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code + IMD_2019_QUINTILES + SMOKING_STATUS + 
    BMI_category + KATZ_INDEX_CAT +  COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + year_tavi_discharge +  # Include year
    binary_known_poor_mobility +  binary_known_poor_mobility_time_interaction +           # Include time interaction effect
    binary_Known_DIABETES + 
    TAVI_procedural_complications + ethnicity_5_group + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + LV_FUNCTION + MITRAL_REGURGITATION
  ,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_HF_primary_outcome_model_2b)



# Test proportional hazards assumption
cox.zph(cox_model_HF_primary_outcome_model_2b)


###NOTES: Proportional hazard assumption test still does not holds true for Global and is violated for binary_known_poor_mobility
##interpretation: Therefore stratify binary_known_poor_mobility



#######################################################################################################################################################################################################
# Fit Cox Proportional Hazards Model 2c
#######################################################################################################################################################################################################
#Fit the Cox proportional hazards model with the coxph function
### stratify binary_known_poor_mobility


# Fit Cox model with calendar time as the timescale, with all covariates, including region_name (REGION HAS BEEN ADDED)
cox_model_HF_primary_outcome_model_2c <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code + IMD_2019_QUINTILES + SMOKING_STATUS + 
    BMI_category + KATZ_INDEX_CAT +  COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + year_tavi_discharge +  # Include year
    strata(binary_known_poor_mobility, region_name) +           # stratify hazard assumption violating variable 
    binary_Known_DIABETES + 
    TAVI_procedural_complications + ethnicity_5_group + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + LV_FUNCTION + MITRAL_REGURGITATION
  ,
  data = clean_cohort
)

# Summary of the Cox model
summary(cox_model_HF_primary_outcome_model_2c)



# Test proportional hazards assumption
cox.zph(cox_model_HF_primary_outcome_model_2c)



###NOTES: Proportional hazard assumption test holds true for Global and holds true for all other variables 

#NOTES: P-values are greater than 0.05 for Global, hold true for all other variables




# Plot survival curve:
plot(survfit(cox_model_HF_primary_outcome_model_2c), main = "Survival Curve for HF Rehospitalization",
     xlab = "Starting from days since TAVI discharge plus 180 days to the date of TAVI discharge Plus 180 days and ending at end of followup (31 March 2024)",
     ylab = "Survival Probability")




#######################################################################################################################################################################################################
# Extract model summary Model 2c
#######################################################################################################################################################################################################


# Extract model summary
cox_summary_HF_primary_outcome_model_2c <- summary(cox_model_HF_primary_outcome_model_2c)

# Extract coefficients and confidence intervals
hr_results_HF_rehosp <- data.frame(
  Variable = rownames(cox_summary_HF_primary_outcome_model_2c$conf.int),
  HR = cox_summary_HF_primary_outcome_model_2c$conf.int[, "exp(coef)"],
  Lower_95CI = cox_summary_HF_primary_outcome_model_2c$conf.int[, "lower .95"],
  Upper_95CI = cox_summary_HF_primary_outcome_model_2c$conf.int[, "upper .95"],
  p_value = cox_summary_HF_primary_outcome_model_2c$coefficients[, "Pr(>|z|)"]
)

# Print results
print(hr_results_HF_rehosp)



# Define file path to save the CSV
output_file_hr_HF_rehosp_start_at_tavi_180_fully_adjusted <- "HR_results_HF_rehosp_start_at_tavi_180_fully_adjusted.csv"

# Save the data frame as a CSV file
write.csv(hr_results_HF_rehosp, file = output_file_hr_HF_rehosp_start_at_tavi_180_fully_adjusted, row.names = TRUE)


# NOTES INTERPRETATION OF MAIN FINDING: Adjusting for all key covariates, accounting for baseline calendar time CR, starting from days since tavi discharge plus 180days, exposure is not associated with the rate of HF rehosp 




#Broom package extraction 
library(broom)
library(dplyr)
library(readr)

# Tidy the Cox model with HR and CIs
cox_results_summary_HF_primary_outcome_model_2c <- tidy(
  cox_model_HF_primary_outcome_model_2c,
  exponentiate = TRUE,  # gives HR instead of log(HR)
  conf.int = TRUE
)

# Optional: clean and round columns
cox_results_summary_HF_primary_outcome_model_2c_clean_fully_adjusted <- cox_results_summary_HF_primary_outcome_model_2c %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))


# Save the rounded table as a CSV file for extraction
output_file_cox_results_summary_HF_primary_outcome_model_2c_clean_fully_adjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/HF Rehosp/cox_results_summary_HF_primary_outcome_model_2c_clean_fully_adjusted.csv"
write.csv(cox_results_summary_HF_primary_outcome_model_2c_clean_fully_adjusted, file = output_file_cox_results_summary_HF_primary_outcome_model_2c_clean_fully_adjusted, row.names = FALSE)








##NO NEED FOR THE BELOW

#################################################
### PLOT THE FULLY ADJUSTED SURVIVAL  CURVES ####
#################################################

###Discuss with Angela how to plot adjusted curve if at all???


# Define a new dataset for adjusting survival curve, adjusting for all covarites and splitting by exposed to rehab vs not exposed
surv_HF_cox_new_data_model_2a <- data.frame(        #generates a simple dataset for adjustment, applying predefined values for covariates
  rehab_indicator = factor(c(0, 1), levels = c(0, 1), 
                           labels = c("Not Exposed to Rehab", "Exposed to Rehab")),
  age_at_tavi = mean(clean_cohort$age_at_tavi, na.rm = TRUE),  # Adjust for mean age
  sex_code = factor("Male", levels = levels(clean_cohort$sex_code)),   # Adjust sex to males 
  ethnicity_5_group = factor(levels(clean_cohort$ethnicity_5_group)[1], levels = levels(clean_cohort$ethnicity_5_group)), #first level is "white"
  IMD_2019_QUINTILES = factor(levels(clean_cohort$IMD_2019_QUINTILES)[1], levels = levels(clean_cohort$IMD_2019_QUINTILES)), #first level is "1" - the most deprived 20% of small areas in the UK
  SMOKING_STATUS = factor(levels(clean_cohort$SMOKING_STATUS)[1], levels = levels(clean_cohort$SMOKING_STATUS)),  #first level is "Current smoker"
  BMI_category = factor("Obesity", levels = levels(clean_cohort$BMI_category)),  # selecting Obesity to adjust with
  NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE = factor("Marked limitation of ordinary physical activity", 
                                               levels = levels(clean_cohort$NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE)), # selecting "Marked limitation of ordinary physical activity" to adjust with
  KATZ_INDEX_CAT = factor(levels(clean_cohort$KATZ_INDEX_CAT)[1], levels = levels(clean_cohort$KATZ_INDEX_CAT)), #first level is "Moderate Impairment, Severe Functional Impairment up to Dependent"
  binary_known_poor_mobility = 1,
  COMORBIDITY_pulmonary_liver_neuro_excardiacvascu = 1,
  binary_Known_DIABETES = 1,
  previous_MI_time_grouped = factor(levels(clean_cohort$previous_MI_time_grouped)[1], levels = levels(clean_cohort$previous_MI_time_grouped)), #first level is "MI > 90 days"
  binary_known_prev_cardiac_surgery = 1,
  LV_FUNCTION = factor("Fair (LVEF = 30-49%)", levels = levels(clean_cohort$LV_FUNCTION)),  # selecting ("Fair (LVEF = 30-49%)" to adjust with
  MITRAL_REGURGITATION = factor("Mild", levels = levels(clean_cohort$MITRAL_REGURGITATION)),  # selecting "Mild" mitral regurgitation to adjust the model with
  TAVI_procedural_complications = 1,
  year_tavi_discharge = factor("2018", levels = levels(clean_cohort$year_tavi_discharge)), # year 2018 (pre COVID) is selected to adjust the cohort with
  stringsAsFactors = FALSE
)



#Compute Survival Estimates for Adjusted Model
surv_fit_HF_cox_fully_adjusted_model_2a <- survfit(
  cox_model_HF_primary_outcome_model_2a, 
  newdata = surv_HF_cox_new_data_model_2a
)


print(surv_fit_HF_cox_fully_adjusted_model_2a)


# Plot the adjusted survival curves 3b
plot(
  surv_fit_HF_cox_fully_adjusted_model_2a,
  col = c("red", "blue"),            # Colors for No Rehab and Rehab
  lwd = 2,                           # Line width
  lty = 1,                           # Line type
  mark.time = TRUE,                  # Add censoring marks
  main = "Fully Adjusted Survival Curve for HF Rehospitalization (Rehospitlisation post 180 days post discharge)",
  xlab = "Calendar Time, starting from days since TAVI discharge plus 180 days and ending at end of followup (31 March 2024)",
  ylab = "Survival Probability",
  conf.int = FALSE                   # Disable confidence intervals
)

# Add a legend
legend(
  "bottomleft",
  legend = c("No Rehab (Adjusted)", "Rehab (Adjusted)"),
  col = c("red", "blue"),
  lwd = 2,
  bty = "n"                          # No border around the legend
)



#PLot survival curve
library(survminer)

ggsurvplot(
  surv_fit_HF_cox_fully_adjusted_model_3d,
  data = split_data2,
  conf.int = TRUE,  # Show confidence intervals
  pval = FALSE,      
  risk.table = TRUE, 
  legend.title = "Rehab Exposure",
  legend.labs = c("Not Exposed to Rehab", "Exposed to Rehab"),
  palette = c("#E69F00", "#0072B2"),  # Distinct colors for rehab vs non-rehab
  linetype = c("solid", "dashed"),  # Solid for rehab, dashed for non-rehab
  break.time.by = 180,  # Adjust x-axis breakpoints for readability
  title = "Fully Adjusted Survival Curve for HF Rehospitalization (HF Rehospitilisations from 180 days post TAVI discharge)",
  xlab = "Calendar Time, starting from days since TAVI discharge plus 180 days and ending at end of followup (31 March 2024)",
  ylab = "Survival Probability",
  ggtheme = theme_minimal()
)


#check it is correct
head(surv_HF_cox_new_data_model_3d)












######################################################################################################################################################################################################################
## New visualization for model 3c: creates and visualizes Kaplan-Meier and Cox-adjusted survival curves for heart failure (HF) rehospitalization, stratified by cardiac rehab exposure.
########################################################################################################################################################################################################################

#A survival analysis showing real-world Kaplan-Meier estimates vs. adjusted Cox regression survival curves, making it clear how cardiac rehab is associated with HF rehospitalization.
# Fits survival models (Kaplan-Meier for observed data and Cox for adjusted data).
library(survival)
library(survminer)

# Check for events in the dataset
table(split_data2$rehosp_indicator)



# Kaplan-Meier Survival Model (Observed Data): creates an observed Kaplan-Meier survival model for HF rehospitalization.
km_fit_HF_rehosp <- survfit(Surv(period_start, period_end, rehosp_indicator) ~ rehab_indicator, data = split_data2) 
#NOTES:
  #Defines survival time from period_start to period_end,  rehosp_indicator = 1 means an event (rehospitalization) occurred.
  #~ rehab_indicator: Stratifies survival by whether the patient attended cardiac rehab


# Cox Model Adjusted Survival: Fit Cox Proportional Hazards Model (Adjusted Survival)
cox_fit_HF_rehosp <- survfit(cox_model_timevaryingCR_HF_primary_outcome_model_3a, newdata = surv_HF_cox_new_data_model_3d)
#NOTES: 
#   This code, Uses a previously fitted Cox proportional hazards model (cox_model_timevaryingCR_HF_primary_outcome_model_3b).
#   newdata = surv_HF_cox_new_data_model_3b: Computes survival estimates using a dataset with predefined covariates (adjusting for baseline differences).
#   This produces an adjusted survival curve accounting for covariates



# Create Kaplan-Meier (Observed) plot: This generates the unadjusted Kaplan-Meier survival curve.
km_plot_HF_rehosp <- ggsurvplot(
  km_fit_HF_rehosp, 
  data = split_data2,
  conf.int = TRUE, 
  pval = FALSE,
  risk.table = FALSE, 
  legend.title = "Rehab Exposure",
  legend.labs = c("No Rehab", "Exposed to Rehab"),
  palette = c("#E69F00", "#0072B2"), 
  linetype = c("solid", "dashed"),
  censor.shape = "|", 
  censor.size = 3,
  break.time.by = 180,
  title = "Unadjusted Kaplan-Meier Survival Curve for HF Rehospitalization (HF Rehospitilisations from 180 days post TAVI discharge)",
  xlab = "Calendar Time, starting from days since TAVI discharge plus 180 days and ending at end of followup (31 March 2024)",
  ylab = "Survival Probability",
  ggtheme = theme_minimal()
)

# Create Adjusted Cox Model survival plot:
#Plots the Cox-adjusted survival curve.
# Uses adjusted data (surv_HF_cox_new_data_model_3c).
# Keeps the same visual style as the Kaplan-Meier plot.
# This produces the adjusted survival curve based on Cox regression

cox_plot_HF_rehosp <- ggsurvplot(
  cox_fit_HF_rehosp, 
  data = surv_HF_cox_new_data_model_3d,
  conf.int = TRUE, 
  pval = FALSE,
  risk.table = FALSE, 
  legend.title = "Rehab Exposure (Adjusted)",
  legend.labs = c("No Rehab (Adjusted)", "Exposed to Rehab (Adjusted)"),
  palette = c("#E69F00", "#0072B2"), 
  linetype = c("solid", "dashed"),
  censor.shape = "|", 
  censor.size = 3,
  break.time.by = 180,
  title = "Cox-Adjusted Survival Curve for HF Rehospitalization (HF Rehospitilisations from 180 days post TAVI discharge)",
  xlab = "Calendar Time, starting from days since TAVI discharge plus 180 days and ending at end of followup (31 March 2024)",
  ylab = "Survival Probability",
  ggtheme = theme_minimal()
)

# Combine both plots
library(ggpubr)
combined_plot_HF_rehosp <- ggarrange(km_plot_HF_rehosp$plot, cox_plot_HF_rehosp$plot, 
                           ncol = 2, 
                           labels = c("A", "B"))

# Display the combined plot: The Kaplan-Meier and Cox-adjusted survival curves are displayed side by side.
#The observed (unadjusted) and adjusted survival curves are displayed together, allowing direct comparison of real-world and adjusted estimates.
print(combined_plot_HF_rehosp)












##########################################################################################################################################################################################

