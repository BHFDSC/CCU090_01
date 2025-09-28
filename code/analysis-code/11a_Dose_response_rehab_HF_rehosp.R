
##############################################################################################################################################################################################
#                                                                          DOSE RESPONSE MODELS
#   Test and visualize if there is a dose-response relationship between the number of cardiac rehab appointments attended (exp_cardiac_rehab_attends_6m) and the reduction in rehospitalizations for HF
#   
###############################################################################################################################################################################################



#####################################################
##prepare data ##
#####################################################

# Load necessary libraries
library(survival)
library(dplyr)
library(lubridate)

#Create Relevant Variables for:
# Time variable: The time to rehospitalization or censoring. rehosp_censor_time: days since TAVI discharge+180 days to rehospitalization for HF or censoring (death or end of study period)
# Event variable: A binary variable (1 = HF rehosp, 0 = censored).
# Exposure variable: Rehab appointments attended (number of sessions attended).
# Covariates: Adjust for other variables (e.g., age, comorbidities).


##Explanation:
# Define study start and end dates
study_start_new <- clean_cohort$tavi_cips_disdate_plus180days
study_end <- as.Date("2024-03-31")

# Add necessary variables
clean_cohort_dose_response <- clean_cohort %>%
  mutate(
    # Days since study start (calendar time)
    days_since_tavi_dis_plus180 = as.numeric(tavi_cips_disdate_plus180days - study_start_new),
    
    # Time of HF rehospitalization or censoring (death or end of study period)
    hf_rehosp_censor_time = pmin(as.numeric(out_readm_hf_date - study_start_new),         # Days to HF readmission
                                 as.numeric(date_of_death - study_start_new),                # Days to death
                                 as.numeric(study_end - study_start_new), na.rm = TRUE),     # Days to study end
    
    # Rehab exposure indicator
    exp_cardiac_rehab_attends_6m = as.numeric(clean_cohort$exp_cardiac_rehab_attends_6m),
    
    # Rehospitalization indicator (1 for HF rehospitalization, 0 for censoring)
    hf_rehosp_indicator = ifelse(is.na(out_readm_hf_flag) | out_readm_hf_flag == 0, 0, 1)
    
  )



# out_readm_hf_date - study_start: Days from TAVI discharge_plus180 until rehospitalization for HF.
#date_of_death - study_start: Days from TAVI discharge_plus180 until death.
#study_end - study_start: Days from TAVI discharge_plus180 until the end of the study.
#pmin() takes the minimum of these times as the censoring point, considering the event that happens earliest (i.e., HF rehospitalization, death, or end of follow-up).


##check if any patients have rehosp or censoring time less than or equal to study start 
count_invalid <- sum(clean_cohort_dose_response$hf_rehosp_censor_time <= clean_cohort_dose_response$days_since_tavi_dis_plus180, na.rm = TRUE)
print(count_invalid)
##NOTES, n= 2190(rounded to nearest 5) patients have have rehosp or censoring time less than or equal to study start (which is tavi dis date plus 180days), only interested in post discharge from TAVIplus 180 so filter out these patients 

# Filter out patients based on rehosp_censor_time conditions
clean_cohort_dose_response <- clean_cohort_dose_response %>%
  filter(!(hf_rehosp_censor_time <= days_since_tavi_dis_plus180))

count_invalid_dose_response <- sum(clean_cohort_dose_response$hf_rehosp_censor_time <= clean_cohort_dose_response$days_since_tavi_dis_plus180, na.rm = TRUE)
print(count_invalid_dose_response)




#### Fit the unadjusted Cox model with continuous rehab attendance ####

hf_cox_model_dose_response_continuous_unadjusted <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ 
    exp_cardiac_rehab_attends_6m,
  data = clean_cohort_dose_response
)

# Summary of the model
summary(hf_cox_model_dose_response_continuous_unadjusted)
# Interpretation: no evidence of a dose response relationship b/w appointments attended and HF rehosp (using continuous scale), although need to adjust for covarites


#Check proportional hazards assumption
cox.zph(hf_cox_model_dose_response_continuous_unadjusted)
# NOTES: Model does not meet the proportional hazards assumptions, Global is less than 0.05 


#extract results for continuous dose response relationship unadjusted model
library(broom)
library(dplyr)
library(readr)

#  model results for extraction
HF_rehosp_results_cox_model_dose_response_continuous_unadjusted <- tidy(hf_cox_model_dose_response_continuous_unadjusted, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Save the rounded table as a CSV file for extraction
output_file_HF_rehosp_results_cox_model_dose_response_continuous_unadjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/HF Rehosp/HF_dose response relationship/HF_rehosp_results_cox_model_dose_response_continuous_unadjusted.csv"
write.csv(HF_rehosp_results_cox_model_dose_response_continuous_unadjusted, file = output_file_HF_rehosp_results_cox_model_dose_response_continuous_unadjusted, row.names = FALSE)




#Plot HR (All-Cause Rehospitalization) vs. Appointments 
library(survival)
library(ggplot2)
library(dplyr)

# Sequence over range of rehab appointments
rehab_range <- seq(min(clean_cohort_dose_response$exp_cardiac_rehab_attends_6m, na.rm = TRUE),
                   max(clean_cohort_dose_response$exp_cardiac_rehab_attends_6m, na.rm = TRUE),
                   length.out = 100)

new_data_dose_response <- data.frame(exp_cardiac_rehab_attends_6m = rehab_range)

# Predict log HR and confidence intervals
pred <- predict(hf_cox_model_dose_response_continuous_unadjusted,
                newdata = new_data_dose_response,
                type = "lp",
                se.fit = TRUE)

# Add to data frame
new_data_dose_response$logHR <- pred$fit
new_data_dose_response$HR <- exp(pred$fit)
new_data_dose_response$lower_CI <- exp(pred$fit - 1.96 * pred$se.fit)
new_data_dose_response$upper_CI <- exp(pred$fit + 1.96 * pred$se.fit)

# Plot
ggplot(new_data_dose_response, aes(x = exp_cardiac_rehab_attends_6m, y = HR)) +
  geom_line(size = 1.2, color = "darkgreen") +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2, fill = "darkgreen") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "Number of Cardiac Rehab Appointments within 6 months",
    y = "Hazard Ratio (HR)",
    title = "Unadjusted Heart failure Rehospitalization: Dose-Response Curve"
  ) +
  theme_minimal(base_size = 16)




###fit the Cox Model with Natural Splines + HR Plot

# Plot the Hazard Ratios for Continuous Dose
# Use a spline regression to model the relationship between the number of appointments attended and the hazard of HF rehospitalization
# Add splines for the continuous variable
#Spline regression is a statistical technique used to model relationships in data that are non-linear
#the relationship between the variables changes at different points along the range of the data.

library(survival)
library(splines)    # for ns()
library(ggplot2)
library(dplyr)

# Fit Cox model with natural spline for rehab attendance
hf_rehosp_dose_response_cox_natural_splines_model <- coxph(Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ 
                                                       ns(exp_cardiac_rehab_attends_6m, df = 3),
                                                     data = clean_cohort_dose_response)


# Summary of the model
summary(hf_rehosp_dose_response_cox_natural_splines_model)
# Interpretation: when using natutral splines to model relationships in data that are non-linear there is no dose response relationship b/w appointments attended and HF rehosp (using continuous scale)


#Check proportional hazards assumption
cox.zph(hf_rehosp_dose_response_cox_natural_splines_model)
# NOTES: Model meets the proportional hazards assumptions, Global is > than 0.05



#extract results for continuous dose response relationship
library(broom)
library(dplyr)
library(readr)

#  model results for extraction
hf_rehosp_results_cox_model_dose_response_continuous_natural_splines_unadjusted <- tidy(hf_rehosp_dose_response_cox_natural_splines_model, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Save the rounded table as a CSV file for extraction
output_file_hf_rehosp_results_cox_model_dose_response_continuous_natural_splines_unadjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/HF Rehosp/HF_dose response relationship/hf_rehosp_results_cox_model_dose_response_continuous_natural_splines_unadjusted.csv"
write.csv(hf_rehosp_results_cox_model_dose_response_continuous_natural_splines_unadjusted, file = output_file_hf_rehosp_results_cox_model_dose_response_continuous_natural_splines_unadjusted, row.names = FALSE)




#Plot HR (HF Rehospitalization) vs. Appointments - natural splines model - non linear relationship

# Generate prediction data
rehab_range <- seq(min(clean_cohort_dose_response$exp_cardiac_rehab_attends_6m, na.rm = TRUE),
                   max(clean_cohort_dose_response$exp_cardiac_rehab_attends_6m, na.rm = TRUE),
                   length.out = 100)

new_data_dose_response <- data.frame(exp_cardiac_rehab_attends_6m = rehab_range)

# Predict log HR and 95% CI
pred <- predict(hf_rehosp_dose_response_cox_natural_splines_model, newdata = new_data_dose_response, type = "lp", se.fit = TRUE)

# Add predictions to data frame
new_data_dose_response$logHR <- pred$fit
new_data_dose_response$HR <- exp(pred$fit)
new_data_dose_response$lower_CI <- exp(pred$fit - 1.96 * pred$se.fit)
new_data_dose_response$upper_CI <- exp(pred$fit + 1.96 * pred$se.fit)

# Plot
ggplot(new_data_dose_response, aes(x = exp_cardiac_rehab_attends_6m, y = HR)) +
  geom_line(color = "purple", size = 1.2) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2, fill = "purple") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "Number of Cardiac Rehab Appointments within 6 months",
    y = "Hazard Ratio (HR)",
    title = "Heart Failure Rehospitalization: Non-linear Dose-Response (Natural Splines)"
  ) +
  theme_minimal(base_size = 16)






# Use a spline regression 2 degrees of freedom
library(splines)

##No Covarites added

cox_model_spline_HF_rehosp_dose_response <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ 
    ns(exp_cardiac_rehab_attends_6m, df = 2),
  data = clean_cohort_dose_response
)

summary(cox_model_spline_HF_rehosp_dose_response)



#extract results for continuous dose response relationship using natural splines
library(broom)
library(dplyr)
library(readr)

#  model results for extraction
results_cox_model_2df_spline_HF_rehosp_dose_response <- tidy(cox_model_spline_HF_rehosp_dose_response, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Save the rounded table as a CSV file for extraction
output_file_results_cox_model_2df_spline_HF_rehosp_dose_response <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/HF Rehosp/HF_dose response relationship/results_cox_model_2df_spline_HF_rehosp_dose_response.csv"
write.csv(results_cox_model_2df_spline_HF_rehosp_dose_response, file = output_file_results_cox_model_2df_spline_HF_rehosp_dose_response, row.names = FALSE)







# Create predicted hazard ratios from the model that used natural splines inorder to plot the HR according to number of appointements attended
library(survival)
new_data_cox_model_spline_HF_rehosp_dose_response <- data.frame(exp_cardiac_rehab_attends_6m = seq(0, max(clean_cohort_dose_response$exp_cardiac_rehab_attends_6m, na.rm = TRUE), by = 1))
new_data_cox_model_spline_HF_rehosp_dose_response$predicted_hr <- predict(cox_model_spline_HF_rehosp_dose_response, newdata = new_data_cox_model_spline_HF_rehosp_dose_response, type = "risk")

# Plot the dose-response relationship
plot(
  new_data_cox_model_spline_HF_rehosp_dose_response$exp_cardiac_rehab_attends_6m, new_data_cox_model_spline_HF_rehosp_dose_response$predicted_hr, 
  type = "l", col = "blue", lwd = 2,
  xlab = "Number of Cardiac Rehab Appointments Attended",
  ylab = "Hazard Ratio (HR)",
  main = "Heart Failure Rehospitilisation: Dose-Response Relationship for Cardiac Rehab"
)
abline(h = 1, col = "red", lty = 2)  # Add a reference line at HR = 1





library(ggplot2)

# Scatter plot of raw rehab attendance vs. HR
ggplot(clean_cohort_dose_response, aes(x = exp_cardiac_rehab_attends_6m, y = hf_rehosp_indicator)) +
  geom_jitter(width = 0.3, height = 0.02, alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  labs(title = "Raw Relationship: Cardiac Rehab vs. HF Readmission",
       x = "Number of Cardiac Rehab Appointments Attended",
       y = "Rehospitalization Indicator (0/1)")




#A U-shaped hazard ratio suggests that both extremes (few or many rehab sessions) may be associated with higher HF rehospitalization risk. Some possible explanations:

#Low Attendance (0-3 sessions):

#  Patients may be too sick to attend more rehab.
# Lack of rehab exposure leads to poorer outcomes.
# High Attendance (>10 sessions):
#    Selection bias: Sicker patients who survive longer may be more likely to attend more sessions.
#    Reverse causality: Those with worse health conditions may be recommended more rehab









#### Fit the adjusted Cox model with continuous rehab attendance - adjusting for covariates - continuous dose response relationship ###
#adjust for covariates

#Also need to account for calendar time at baseline "year of baseline", given that the start time is TAVI dischargePLus180, we need to account for calendar time differences due to COVID
# include calendar time at baseline in the Cox model, need to extract the year of TAVI discharge and include it as a covariate

#Extract Year of TAVI Discharge
library(lubridate)  # Ensure lubridate is loaded for date handling

# Extract the year from TAVI discharge date
clean_cohort_dose_response <- clean_cohort_dose_response %>%
  mutate(year_tavi_discharge = year(tavi_cips_disdate))  # Extract year


#convert to a factor to compare different years independently
#Cox model will estimate a separate hazard ratio for each year relative to the reference year (2018).
clean_cohort_dose_response <- clean_cohort_dose_response %>%
  mutate(year_tavi_discharge = as.factor(year_tavi_discharge))

###cox model adjusting for covariates
cox_model_dose_response_continuous_fully_adjusted <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ 
    exp_cardiac_rehab_attends_6m + age_at_tavi + sex_code + ethnicity_5_group + 
    IMD_2019_QUINTILES + SMOKING_STATUS + BMI_category + 
    NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT + 
    binary_known_poor_mobility + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + 
    MITRAL_REGURGITATION + TAVI_procedural_complications + 
    binary_Known_DIABETES + LV_FUNCTION + year_tavi_discharge,
  data = clean_cohort_dose_response
)

# Summary of the model
summary(cox_model_dose_response_continuous_fully_adjusted)


# Test proportional hazards assumption
cox.zph(cox_model_dose_response_continuous_fully_adjusted)

# NOTES: Model is not running correctly,  


###Refine model by stratifying ethnicity_5_group and SMOKING_STATUS, these were the 2 problematic variables and once stratified the model runs well. And included REGION

cox_model_dose_response_continuous_fully_adjusted <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ 
    exp_cardiac_rehab_attends_6m + age_at_tavi + sex_code + region_name + strata(ethnicity_5_group, SMOKING_STATUS) + KATZ_INDEX_CAT + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + 
    IMD_2019_QUINTILES + BMI_category + 
    binary_known_poor_mobility + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + 
    MITRAL_REGURGITATION + TAVI_procedural_complications + 
    binary_Known_DIABETES + LV_FUNCTION + year_tavi_discharge,
  data = clean_cohort_dose_response
)

# Summary of the model
summary(cox_model_dose_response_continuous_fully_adjusted)

#test hazard assumptions
cox.zph(cox_model_dose_response_continuous_fully_adjusted)

#NOTES: model meets  hazard assumptions test
# Interpretation: after adjustment, there does not appear to be any dose response relationship b/w appointments attended and HF rehosp (using continuous scale)




#extract results for continuous dose response relationship
library(broom)
library(dplyr)
library(readr)

#  model results for extraction
HF_rehosp_results_cox_model_dose_response_continuous_fully_adjusted <- tidy(cox_model_dose_response_continuous_fully_adjusted, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Save the rounded table as a CSV file for extraction
output_file_HF_rehosp_results_cox_model_dose_response_continuous_fully_adjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/HF Rehosp/HF_dose response relationship/HF_rehosp_results_cox_model_dose_response_continuous_fully_adjusted.csv"
write.csv(HF_rehosp_results_cox_model_dose_response_continuous_fully_adjusted, file = output_file_HF_rehosp_results_cox_model_dose_response_continuous_fully_adjusted, row.names = FALSE)








#### Fit the Cox model with rehab attendance categorical ####
#Change variable to number of CR appointments grouped into categories


cox_model_dose_response_categorical <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ 
    rehab_attendance_groups + age_at_tavi + sex_code + strata(ethnicity_5_group, SMOKING_STATUS) +
    IMD_2019_QUINTILES + BMI_category + 
    NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT + 
    binary_known_poor_mobility + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + 
    MITRAL_REGURGITATION + TAVI_procedural_complications + 
    binary_Known_DIABETES + LV_FUNCTION + year_tavi_discharge,
  data = clean_cohort_dose_response
)

# Summary of the model
summary(cox_model_dose_response_categorical)

# Test proportional hazards assumption
cox.zph(cox_model_dose_response_categorical)


#Interpretation: There is no evidence of a dose response relationship using categorical groupings of rehab appointments attended



#extract results for categorical dose response relationship
library(broom)
library(dplyr)
library(readr)

#  model results for extraction
HF_rehosp_results_cox_model_dose_response_categorical_fully_adjusted <- tidy(cox_model_dose_response_categorical, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Save the rounded table as a CSV file for extraction
output_file_HF_rehosp_results_cox_model_dose_response_categorical_fully_adjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/HF Rehosp/HF_dose response relationship/HF_rehosp_results_cox_model_dose_response_categorical_fully_adjusted.csv"
write.csv(HF_rehosp_results_cox_model_dose_response_categorical_fully_adjusted, file = output_file_HF_rehosp_results_cox_model_dose_response_categorical_fully_adjusted, row.names = FALSE)






###Create a Kaplan-Meier Survival Curve by Number of Appointments attended Grouped into categories

# Kaplan-Meier survival curves by Number of Appointments attended Grouped
surv_object_CR_dose_response_categorical <- survfit(Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_attendance_groups, data = clean_cohort_dose_response)

# Plot survival curves
plot(
  surv_object_CR_dose_response_categorical, 
  col = c("red", "blue", "green", "purple", "orange"),  # Colors for each group
  lwd = 2,                                   # Line width
  xlab = "Calendar Time (days), starting from 180 days post TAVI discharge to end of follow-up (31 March 2024)", 
  ylab = "Survival Probability", 
  main = "Unadjusted Survival Curve for Heart Failure Rehospitilisation by Number of Rehab Appointements Attended"
)
legend(
  "bottomleft", 
  legend = levels(factor(clean_cohort$rehab_attendance_groups)), 
  col = c("red", "blue", "green", "purple", "orange"), 
  lwd = 2
)





### Kaplan-Meier Survival Curve: Number of Cardiac Rehab Appointments Attended (Grouped)

# Load necessary packages
library(survival)
library(survminer)

# Fit Kaplan-Meier survival model
surv_object_CR_dose_response_categorical_2 <- survfit(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_attendance_groups, 
  data = clean_cohort_dose_response
)

# Define color palette and line types
color_palette <- c("red", "blue", "green", "purple", "orange")
line_types <- c("solid", "dashed", "dotted", "dotdash", "solid")

# Plot the Kaplan-Meier survival curves
ggsurvplot(
  surv_object_CR_dose_response_categorical_2,
  data = clean_cohort_dose_response,
  conf.int = FALSE,              # Show confidence intervals
  pval = FALSE,                  # Display p-value
  risk.table = FALSE,            # Show risk table below the plot
  risk.table.height = 0.2,      # Adjust risk table height
  risk.table.col = "strata",    # Color risk table by strata
  legend.title = "Rehab Attendance",
  legend.labs = levels(factor(clean_cohort_dose_response$rehab_attendance_groups)), 
  palette = color_palette,      # Assign colors to groups
  linetype = line_types,        # Assign line types to groups
  break.time.by = 180,          # Set time interval for x-axis ticks
  xlab = "Days Since Start of Follow-Up",
  ylab = "Survival Probability",
  title = "Kaplan-Meier Survival Curve for Heart Failure Rehospitilisation by Number of Cardiac Rehab Appointments Attended",
  ggtheme = theme_minimal(base_size = 16)     # Use a clean theme
)





##categorising rehab appointments into quartiles 

clean_cohort_dose_response$rehab_quartiles <- cut(clean_cohort_dose_response$exp_cardiac_rehab_attends_6m,
                                    breaks = quantile(clean_cohort_dose_response$exp_cardiac_rehab_attends_6m, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
                                    include.lowest = TRUE)

cox_model_rehab_quartiles <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_quartiles + age_at_tavi + sex_code + strata (ethnicity_5_group,SMOKING_STATUS) + 
    IMD_2019_QUINTILES + BMI_category + 
    NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT + 
    binary_known_poor_mobility + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + 
    MITRAL_REGURGITATION + TAVI_procedural_complications + 
    binary_Known_DIABETES + LV_FUNCTION + year_tavi_discharge,
  data = clean_cohort_dose_response
)

summary(cox_model_rehab_quartiles)

# Interpretation - no evidence of dose response relationship for number of rehab appointments attended and rate of HF rehospitilisation


