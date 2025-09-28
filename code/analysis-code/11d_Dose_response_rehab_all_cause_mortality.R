##############################################################################################################################################################################################
#                                                                          DOSE RESPONSE MODELS
#   Test and visualize if there is a dose-response relationship between the number of cardiac rehab appointments attended (exp_cardiac_rehab_attends_6m) and rate of all cause mortality
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
# Time variable: The time to mortality or censoring. 
# Event variable: A binary variable (1 = all cause rehosp, 0 = censored).
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
    
    # Time of death or censoring (death or end of study period)
    all_cause_mortality_censor_time = pmin(as.numeric(date_of_death - study_start_new),            # Days to death
                                        as.numeric(study_end - study_start_new), na.rm = TRUE),     # Days to study end
    
    # Rehab exposure indicator
    exp_cardiac_rehab_attends_6m = as.numeric(clean_cohort$exp_cardiac_rehab_attends_6m),
    
    # mortality indicator (1 for all cause mortality, 0 for censoring)
    all_cause_mortality_indicator = ifelse(!is.na(date_of_death), 1, 0)
    
  )



##check if any patients die within 180 days of discharge 
count_invalid <- sum(clean_cohort_dose_response$all_cause_mortality_censor_time <= clean_cohort_dose_response$days_since_tavi_dis_plus180, na.rm = TRUE)
print(count_invalid)
##NOTES, n=1295 (rounded to nearest 5)  patients within 180 days post tavi discharge, only interested in post discharge from TAVIplus 180 so filter out these patients 

# Filter out patients
clean_cohort_dose_response <- clean_cohort_dose_response %>%
  filter(!(all_cause_mortality_censor_time <= days_since_tavi_dis_plus180))

count_invalid <- sum(clean_cohort_dose_response$all_cause_mortality_censor_time <= clean_cohort_dose_response$days_since_tavi_dis_plus180, na.rm = TRUE)
print(count_invalid)

##NOTES, n=1295 (rounded to nearest 5)  patients filtered out abnd not incldued in analysis for mortality b/c not relevant to study due to death within 180 days





#### Fit the unadjusted Cox model with continuous rehab attendance####

all_cause_mortality_cox_model_dose_response_continuous_unadjusted <- coxph(
  Surv(all_cause_mortality_censor_time, all_cause_mortality_indicator) ~ 
    exp_cardiac_rehab_attends_6m,
  data = clean_cohort_dose_response
)

# Summary of the model
summary(all_cause_mortality_cox_model_dose_response_continuous_unadjusted)
# Interpretation: before adjustment, a dose response relationship exists b/w appointments attended and all cause mortality (using continuous scale), although need to adjust for covarites
##for every increase in CR appointment attended, the rate of all cause mortality decreases by 7% (p=0.024)
#this assumes there is a linear relationship - see below results for natural splines (which does not assume linear)

#Check proportional hazards assumption
cox.zph(all_cause_mortality_cox_model_dose_response_continuous_unadjusted)
# NOTES: Model meets the proportional hazards assumptions, Global is > than 0.05 



#extract results for continuous dose response relationship
library(broom)
library(dplyr)
library(readr)

#  model results for extraction
all_cause_mortality_results_cox_model_dose_response_continuous_unadjusted <- tidy(all_cause_mortality_cox_model_dose_response_continuous_unadjusted, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Save the rounded table as a CSV file for extraction
output_file_all_cause_mortality_cox_model_dose_response_continuous_unadjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_mortality/Mortality Dose response/all_cause_mortality_results_cox_model_dose_response_continuous_unadjusted.csv"
write.csv(all_cause_mortality_results_cox_model_dose_response_continuous_unadjusted, file = output_file_all_cause_mortality_cox_model_dose_response_continuous_unadjusted, row.names = FALSE)



#Plot HR (All-Cause mortality) vs. Appointments 
library(survival)
library(ggplot2)
library(dplyr)

# Sequence over range of rehab appointments
rehab_range <- seq(min(clean_cohort_dose_response$exp_cardiac_rehab_attends_6m, na.rm = TRUE),
                   max(clean_cohort_dose_response$exp_cardiac_rehab_attends_6m, na.rm = TRUE),
                   length.out = 100)

new_data_dose_response <- data.frame(exp_cardiac_rehab_attends_6m = rehab_range)

# Predict log HR and confidence intervals
pred <- predict(all_cause_mortality_cox_model_dose_response_continuous_unadjusted,
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
    title = "Unadjusted All-Cause Mortality: Dose-Response Curve"
  ) +
  theme_minimal(base_size = 16)




###Cox Model with Natural Splines + HR Plot

# Plot the Hazard Ratios for Continuous Dose
# Use a spline regression to model the relationship between the number of appointments attended and the hazard of death
# Add splines for the continuous variable
#Spline regression is a statistical technique used to model relationships in data that are non-linear
#the relationship between the variables changes at different points along the range of the data.

library(survival)
library(splines)    # for ns()
library(ggplot2)
library(dplyr)

# Fit Cox model with natural spline for rehab attendance
all_cause_mortality_dose_response_cox_ns_model <- coxph(Surv(all_cause_mortality_censor_time, all_cause_mortality_indicator) ~ 
                                                       ns(exp_cardiac_rehab_attends_6m, df = 3),
                                                     data = clean_cohort_dose_response)


# Summary of the model
summary(all_cause_mortality_dose_response_cox_ns_model)
# Interpretation: when using natutral splines to model relationships in data that are non-linear there is a dose response relationship b/w appointments attended and all cause mortality (using continuous scale)


#Check proportional hazards assumption
cox.zph(all_cause_mortality_dose_response_cox_ns_model)
# NOTES: Model meets the proportional hazards assumptions, Global is > than 0.05



#extract results for continuous dose response relationship
library(broom)
library(dplyr)
library(readr)

#  model results for extraction
all_cause_mortality_results_cox_model_dose_response_continuous_natural_splines_unadjusted <- tidy(all_cause_mortality_dose_response_cox_ns_model, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Save the rounded table as a CSV file for extraction
output_file_all_cause_mortality_results_cox_model_dose_response_continuous_natural_splines_unadjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_mortality/Mortality Dose response/all_cause_mortality_results_cox_model_dose_response_continuous_natural_splines_unadjusted.csv"
write.csv(all_cause_mortality_results_cox_model_dose_response_continuous_natural_splines_unadjusted, file = output_file_all_cause_mortality_results_cox_model_dose_response_continuous_natural_splines_unadjusted, row.names = FALSE)




#Plot HR (All-Cause mortality) vs. Appointments - natutral splines model - non linear relationship

# Generate prediction data
rehab_range <- seq(min(clean_cohort_dose_response$exp_cardiac_rehab_attends_6m, na.rm = TRUE),
                   max(clean_cohort_dose_response$exp_cardiac_rehab_attends_6m, na.rm = TRUE),
                   length.out = 100)

new_data_dose_response <- data.frame(exp_cardiac_rehab_attends_6m = rehab_range)

# Predict log HR and 95% CI
pred <- predict(all_cause_mortality_dose_response_cox_ns_model, newdata = new_data_dose_response, type = "lp", se.fit = TRUE)

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
    title = "All-Cause Mortality: Non-linear Dose-Response (Natural Splines)"
  ) +
  theme_minimal(base_size = 16)






#### Fit the adjusted Cox model with continuous rehab attendance - adjusting for covariates - continuous dose response relationship ###

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
all_cause_mortality_cox_model_dose_response_continuous_fully_adjusted <- coxph(
  Surv(all_cause_mortality_censor_time, all_cause_mortality_indicator) ~ 
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
summary(all_cause_mortality_cox_model_dose_response_continuous_fully_adjusted)


# Test proportional hazards assumption
cox.zph(all_cause_mortality_cox_model_dose_response_continuous_fully_adjusted)

# NOTES: Model runs but tavi complications and year of discharge do not meet hazard assumption test, therefore stratify these


###Refine model by stratifying tavi complications and year of discharge and include REGION

all_cause_mortality_cox_model_dose_response_continuous_fully_adjusted <- coxph(
  Surv(all_cause_mortality_censor_time, all_cause_mortality_indicator) ~ 
    exp_cardiac_rehab_attends_6m + age_at_tavi + sex_code + ethnicity_5_group + region_name + SMOKING_STATUS + KATZ_INDEX_CAT + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + 
    IMD_2019_QUINTILES + BMI_category + binary_known_poor_mobility + MITRAL_REGURGITATION +
    strata(TAVI_procedural_complications, year_tavi_discharge) + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + 
     + 
    binary_Known_DIABETES + LV_FUNCTION,
  data = clean_cohort_dose_response
)

# Summary of the model
summary(all_cause_mortality_cox_model_dose_response_continuous_fully_adjusted)

#test hazard assumptions
cox.zph(all_cause_mortality_cox_model_dose_response_continuous_fully_adjusted)

#NOTES: model meets  hazard assumptions test
# Interpretation: even after adjustment, there does  appear to be any dose response relationship b/w appointments attended and all-cause mortality (using continuous scale). before and after adjutsment there is a dose response




#extract results for continuous dose response relationship fully adjusted
library(broom)
library(dplyr)
library(readr)

#  model results for extraction
all_cause_mortality_results_cox_model_dose_response_continuous_fully_adjusted <- tidy(all_cause_mortality_cox_model_dose_response_continuous_fully_adjusted, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Save the rounded table as a CSV file for extraction
output_file_all_cause_mortality_results_cox_model_dose_response_continuous_fully_adjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_mortality/Mortality Dose response/all_cause_mortality_results_cox_model_dose_response_continuous_fully_adjusted.csv"
write.csv(all_cause_mortality_results_cox_model_dose_response_continuous_fully_adjusted, file = output_file_all_cause_mortality_results_cox_model_dose_response_continuous_fully_adjusted, row.names = FALSE)







#### Fit the Cox model with rehab attendance categorical ####
#Change variable to number of CR appointments grouped into categories

all_cause_mortality_cox_model_dose_response_categorical_unadjusted <- coxph(
  Surv(all_cause_mortality_censor_time, all_cause_mortality_indicator) ~ 
    rehab_attendance_groups,
  data = clean_cohort_dose_response
)

# Summary of the model
summary(all_cause_mortality_cox_model_dose_response_categorical_unadjusted)

# Test proportional hazards assumption
cox.zph(all_cause_mortality_cox_model_dose_response_categorical_unadjusted)

##model meets hazard assumption test
#Interpretation: number of CR appointments grouped into categories shows evidence of dose response beyond 6 appointemnts

#  model results for extraction
all_cause_mortality_results_cox_model_dose_response_categorical_unadjusted <- tidy(all_cause_mortality_cox_model_dose_response_categorical_unadjusted, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Save the rounded table as a CSV file for extraction
output_file_all_cause_mortality_results_cox_model_dose_response_categorical_unadjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_mortality/Mortality Dose response/all_cause_mortality_results_cox_model_dose_response_categorical_unadjusted.csv"
write.csv(all_cause_mortality_results_cox_model_dose_response_categorical_unadjusted, file = output_file_all_cause_mortality_results_cox_model_dose_response_categorical_unadjusted, row.names = FALSE)







##Fully adjusted model categorial dose response
all_cause_mortality_cox_model_dose_response_categorical_fully_adjusted <- coxph(
  Surv(all_cause_mortality_censor_time, all_cause_mortality_indicator) ~ 
    rehab_attendance_groups + age_at_tavi + sex_code + ethnicity_5_group + SMOKING_STATUS +
    IMD_2019_QUINTILES + BMI_category + binary_known_poor_mobility + MITRAL_REGURGITATION +
    NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT + 
    strata (TAVI_procedural_complications, year_tavi_discharge) + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery +
    binary_Known_DIABETES + LV_FUNCTION,
  data = clean_cohort_dose_response
)

# Summary of the model
summary(all_cause_mortality_cox_model_dose_response_categorical_fully_adjusted)

# Test proportional hazards assumption
cox.zph(all_cause_mortality_cox_model_dose_response_categorical_fully_adjusted)
 

#Interpretation: after adjustment, there is no evidence of a dose response relationship using categorical groupings of rehab appointments attended



#extract results for categorical dose response relationship
library(broom)
library(dplyr)
library(readr)

#  model results for extraction
all_cause_mortality_results_cox_model_dose_response_categorical_fully_adjusted <- tidy(all_cause_mortality_cox_model_dose_response_categorical_fully_adjusted, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Save the rounded table as a CSV file for extraction
output_file_all_cause_mortality_results_cox_model_dose_response_categorical_fully_adjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_mortality/Mortality Dose response/all_cause_mortality_results_cox_model_dose_response_categorical_fully_adjusted.csv"
write.csv(all_cause_mortality_results_cox_model_dose_response_categorical_fully_adjusted, file = output_file_all_cause_mortality_results_cox_model_dose_response_categorical_fully_adjusted, row.names = FALSE)





###Create a Kaplan-Meier Survival Curve by Number of Appointments attended Grouped into categories

# Kaplan-Meier survival curves by Number of Appointments attended Grouped
all_cause_mortality_surv_object_CR_dose_response_categorical <- survfit(Surv(all_cause_mortality_censor_time, all_cause_mortality_indicator) ~ rehab_attendance_groups, data = clean_cohort_dose_response)

# Plot survival curves
plot(
  all_cause_mortality_surv_object_CR_dose_response_categorical, 
  col = c("red", "blue", "green", "purple", "orange"),  # Colors for each group
  lwd = 2,                                   # Line width
  xlab = "Calendar Time (days), starting from 180 days post TAVI discharge to end of follow-up (31 March 2024)", 
  ylab = "Survival Probability", 
  main = "Unadjusted Survival Curve for All-Cause Mortality by Number of Rehab Appointements Attended"
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
all_cause_mortality_surv_object_CR_dose_response_categorical_2 <- survfit(
  Surv(all_cause_mortality_censor_time, all_cause_mortality_indicator) ~ rehab_attendance_groups, 
  data = clean_cohort_dose_response
)

# Define color palette and line types
color_palette <- c("red", "blue", "green", "purple", "orange")
line_types <- c("solid", "dashed", "dotted", "dotdash", "solid")

# Plot the Kaplan-Meier survival curves
ggsurvplot(
  all_cause_mortality_surv_object_CR_dose_response_categorical_2,
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
  xlab = "Calendar Time (days), starting from 180 days post TAVI discharge to end of follow-up (31 March 2024)",
  ylab = "Survival Probability",
  title = "Kaplan-Meier Survival Curve for All-Cause Mortality by Number of Cardiac Rehab Appointments Attended",
  ggtheme = theme_minimal(base_size = 16)     # Use a clean theme
)




