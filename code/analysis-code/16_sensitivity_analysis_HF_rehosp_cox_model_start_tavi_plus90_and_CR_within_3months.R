#########################################################################################################################################
#####              Sensitivity analysis - changing rehab definition to within 3 months (vs the base case of 6 months) and HF rehosps within 90days ###############
#####                                 Primary outcome - hospital readmission with HF, cox model                           ###############
#########################################################################################################################################



# Summary volumes b/w rehab within 6 months and rehab within 3months
summary(clean_cohort$binary_exposed_rehab_3m)
summary(clean_cohort$binary_exposed_rehab) #rehab within 6months (base case analysis for whole project)

# Tabulate the counts
CR_tab_3m <- table(clean_cohort$binary_exposed_rehab_3m)
CR_tab_6m <- table(clean_cohort$binary_exposed_rehab)

# Round to nearest 5
round_to_5 <- function(x) round(x / 5) * 5

# Combine into a data frame
CR_definition_summary_table <- data.frame(
  Rehab_Within_3m = round_to_5(CR_tab_3m),
  Rehab_Within_6m = round_to_5(CR_tab_6m)
)



# Label rows
rownames(CR_definition_summary_table) <- names(CR_tab_3m)

CR_definition_summary_table

### 91% of all CR exposures occurred within 3months, only 9% (95) occurred b/w 3-6months - this gives confidence to the base case analysis results for the whole study



# Save the rounded table as a CSV file for extraction
output_file_CR_definition_summary_table_summary <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Sensitivity analysis_3months/CR_definition_summary_table.csv"
write.csv(CR_definition_summary_table, file = output_file_CR_definition_summary_table_summary, row.names = FALSE)




###SENSITIVITY ANALYSIS CR within 3 months and HF rehosps beyond 90days post TAVI
###MAKE SURE RESET THE DATA SET BEFORE RUNNING THIS CODE

##Cox model with cardiac rehabilitation (CR) as independent variable  using calendar time as the timescale, 
##The observation period is from TAVI discharge + 90days to 31 March 2024 - adjusting for calendar time 
#The model will evaluate how the patient's hazard (instantaneous rate) for rehospitalization for HF changes based on their participation in rehab

#####################################################
##Step 1: Define and create Calendar Time Variables##
#####################################################

# Load necessary libraries
library(survival)
library(dplyr)
library(lubridate)

#Create Relevant Variables for:
# Time variable: The time to rehospitalization or censoring. rehosp_censor_time: days since TAVI discharge+90 days to rehospitalization for HF or censoring (death or end of study period)
# Event variable: A binary variable (1 = HF rehosp, 0 = censored).
# Exposure variable: Rehab status (1 = rehab, 0 = no rehab).
# Covariates: Adjust for other variables (e.g., age, comorbidities).


##Explanation:
# Define study start and end dates
study_start_new <- clean_cohort$tavi_cips_disdate_plus90_days
study_end <- as.Date("2024-03-31")

# Add necessary variables
clean_cohort <- clean_cohort %>%
  mutate(
    # Days since study start (calendar time)
    days_since_tavi_dis_plus90 = as.numeric(tavi_cips_disdate_plus90_days - study_start_new),
    
    # Time of HF rehospitalization or censoring (death or end of study period)
    hf_rehosp_censor_time = pmin(as.numeric(out_readm_hf_date - study_start_new),         # Days to HF readmission
                                 as.numeric(date_of_death - study_start_new),                # Days to death
                                 as.numeric(study_end - study_start_new), na.rm = TRUE),     # Days to study end
    
    # Rehab exposure indicator
    rehab_indicator = ifelse(is.na(exp_cardiac_rehab_3m) | exp_cardiac_rehab_3m == 0, 0, 1),
    
    # Rehospitalization indicator (1 for HF rehospitalization, 0 for censoring)
    hf_rehosp_indicator = ifelse(is.na(out_readm_hf_flag) | out_readm_hf_flag == 0, 0, 1)
    
  )



# out_readm_hf_date - study_start: Days from TAVI discharge_plus90 until rehospitalization for HF.
#date_of_death - study_start: Days from TAVI discharge_plus90 until death.
#study_end - study_start: Days from TAVI discharge_plus90 until the end of the study.
#pmin() takes the minimum of these times as the censoring point, considering the event that happens earliest (i.e., HF rehospitalization, death, or end of follow-up).


##check if any patients have rehosp or censoring time less than or equal to study start 
count_invalid <- sum(clean_cohort$hf_rehosp_censor_time <= clean_cohort$days_since_tavi_dis_plus90, na.rm = TRUE)
print(count_invalid)
##NOTES, n= 1450 (rounded to nearest 5) patients have have rehosp or censoring time less than or equal to study start (which is tavi dis date plus 90days), only interested in post discharge from TAVIplus 90days so filter out these patients 
##NOTES, BASE CASE ANALYSIS was n= 2190 (rounded to nearest 5) patients who had rehosp or censoring time less than or equal to study start (which is tavi dis date plus 180days) 
##therefore 740 (rounded to nearest 5) patients had HF rehosp b/w 90days and 180days - this is the difference in the sensitivity analysis


# Filter out patients based on rehosp_censor_time conditions
clean_cohort <- clean_cohort %>%
  filter(!(hf_rehosp_censor_time <= days_since_tavi_dis_plus90))

count_invalid <- sum(clean_cohort$hf_rehosp_censor_time <= clean_cohort$days_since_tavi_dis_plus90, na.rm = TRUE)
print(count_invalid)





#####################################################################
## Fit Cox Proportional Hazards Model ############
#####################################################################
#Fit the Cox proportional hazards model with the coxph function

# Cox model with calendar time as the timescale
sensitivity_3months_cox_model_CR_HF_primary_outcome <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator,
  data = clean_cohort
)

# Summary of the Cox model
summary(sensitivity_3months_cox_model_CR_HF_primary_outcome)
##INTERPRETATION: Sensitivity analysis, whereby CR exposure is within 3months of tavi and HF rehospititilisations post 90days of discharge - shows that before adjustment, exposure to rehab was associated with a lower rate of HF rehosp



# Test proportional hazards assumption
cox.zph(sensitivity_3months_cox_model_CR_HF_primary_outcome)
#NOTES: P-values are greater than 0.05 and Global is greater than 0.05 therefore proportional hazards assumptions hold true



#extract results
library(broom)
library(dplyr)
library(readr)

# 1. Tidy model results for extraction
results_sensitivity_3months_cox_model_CR_HF_primary_outcome <- tidy(sensitivity_3months_cox_model_CR_HF_primary_outcome, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# 2. Extract number of events from the model summary and round
num_events_cox_model_CR_HF_primary_outcome <- summary(sensitivity_3months_cox_model_CR_HF_primary_outcome)$nevent
rounded_events_cox_model_CR_HF_primary_outcome <- 5 * round(num_events_cox_model_CR_HF_primary_outcome / 5)

# 3. Add as a separate row
events_row_cox_model_CR_HF_primary_outcome <- tibble(
  term = "Number of events (rounded to nearest 5)",
  HR = NA, std.error = NA, statistic = NA,
  p_value = NA, lower_CI = NA, upper_CI = NA,
  extra_info = as.character(rounded_events_cox_model_CR_HF_primary_outcome)
)

# 4. Add column for extra info and combine
results_sensitivity_3months_cox_model_CR_HF_primary_outcome <- results_sensitivity_3months_cox_model_CR_HF_primary_outcome %>%
  mutate(extra_info = NA_character_) %>%
  select(term, HR, std.error, statistic, p_value, lower_CI, upper_CI, extra_info)

# 5. Combine and export
output_sensitivity_3months_cox_model_CR_HF_primary_outcome_unadjusted <- bind_rows(events_row_cox_model_CR_HF_primary_outcome, results_sensitivity_3months_cox_model_CR_HF_primary_outcome)


# Save the rounded table as a CSV file for extraction
output_file_sensitivity_3months_cox_model_CR_HF_primary_outcome_unadjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Sensitivity analysis_3months/output_sensitivity_3months_cox_model_CR_HF_primary_outcome_unadjusted.csv"
write.csv(output_sensitivity_3months_cox_model_CR_HF_primary_outcome_unadjusted, file = output_file_sensitivity_3months_cox_model_CR_HF_primary_outcome_unadjusted, row.names = FALSE)





# Cox model with calendar time as the timescale - adjusting for age and sex
sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code,
  data = clean_cohort
)

# Summary of the Cox model
summary(sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex)
##INTERPRETATION: After adjusting for age and sex, there is no association b/w CR exposure and the rate of HF rehosp (p=0.052)

# Test proportional hazards assumption
cox.zph(sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex)
#NOTES: P-values are greater than 0.05 and Global is greater than 0.05 therefore proportional hazards assumptions hold true



#extract results
library(broom)
library(dplyr)
library(readr)

# Tidy model output (with exponentiated HRs and CIs)
results_sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex <- tidy(sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))

# Get number of events and round to nearest 5
num_events_cox_model_CR_HF_primary_outcome_age_sex <- summary(sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex)$nevent
rounded_events_cox_model_CR_HF_primary_outcome_age_sex <- 5 * round(num_events_cox_model_CR_HF_primary_outcome_age_sex / 5)

# Add number of events as a top row
events_row_cox_model_CR_HF_primary_outcome_age_sex <- tibble(
  term = "Number of events (rounded to nearest 5)",
  HR = NA, std.error = NA, statistic = NA,
  p_value = NA, lower_CI = NA, upper_CI = NA,
  extra_info = as.character(rounded_events_cox_model_CR_HF_primary_outcome_age_sex)
)

# Add an extra_info column to match structure
results_sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex <- results_sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex %>%
  mutate(extra_info = NA_character_) %>%
  select(term, HR, std.error, statistic, p_value, lower_CI, upper_CI, extra_info)

# Combine and export
output_results_sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex_adjusted <- bind_rows(events_row_cox_model_CR_HF_primary_outcome_age_sex, results_sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex)

# Save the rounded table as a CSV file for extraction
output_file_results_sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex_adjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Sensitivity analysis_3months/output_results_sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex_adjusted.csv"
write.csv(output_results_sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex_adjusted, file = output_file_results_sensitivity_3months_cox_model_CR_HF_primary_outcome_age_sex_adjusted, row.names = FALSE)





# Plot survival curve:
plot(survfit(sensitivity_3months_cox_model_CR_HF_primary_outcome), main = "Survival Curve for HF Rehospitalization",
     xlab = "Days since 90 post TAVI discharge to end of followup (31 March 2024)",
     ylab = "Survival Probability")



# Create survival object for individuals in the dataset


sensitivity_3months_surv_object_HF_cox <- survfit(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator, 
  data = clean_cohort)


summary(sensitivity_3months_surv_object_HF_cox)

# Plot survival curve for HF rehospitalization for 2 groups (rehab vs no rehab)
#NOTES THIS IS UNADJUSTED SURVIVAL CURVE
plot(
  sensitivity_3months_surv_object_HF_cox,
  col = c("red", "blue"),  # Colors for the curves
  lwd = 2,                 # Line width
  lty = 1,                 # Line type (solid line)
  mark.time = TRUE,        # Add tick marks at censoring times
  main = "Sensitivity Analysis: Survival Curve for HF Rehospitalization beyond 90days from discharge (Unadjusted)",
  xlab = "Calendar Time (days), starting from 90 days post TAVI discharge to end of follow-up (31 March 2024)",
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



summary(sensitivity_3months_surv_object_HF_cox, times = c(0, 5, 7, 10, 20, 30, 60, 90, 180, 1000, 1500, 2000))




#PLot survival curve
library(survminer)
library(survival)
library(dplyr)
library(ggplot2)
library(gridExtra)

sensitivity_3months_surv_plot_HF_cox <- ggsurvplot(
  sensitivity_3months_surv_object_HF_cox,
  data = clean_cohort,
  conf.int = TRUE,  # Show confidence intervals
  pval = FALSE,      
  risk.table = FALSE, 
  legend.title = "Cardiac Rehabilitation Exposure",
  legend.labs = c("Not Exposed to Rehab", "Exposed to Rehab"),
  palette = c("#E69F00", "#0072B2"),  # Distinct colors for rehab vs non-rehab
  linetype = c("solid", "dashed"),  # Solid for rehab, dashed for non-rehab
  break.time.by = 180,  # Adjust x-axis breakpoints for readability
  title = "Sensitivity Analysis: Survival Curve for Heart Failure Rehospitalization beyond 90days from discharge (Unadjusted)",
  xlab = "Calendar Time (days), starting from 90 days post TAVI discharge to end of follow-up (31 March 2024)",
  ylab = "Survival Probability",
  ggtheme = theme_minimal(base_size = 16)
)

sensitivity_3months_surv_plot_HF_cox$plot <- sensitivity_3months_surv_plot_HF_cox$plot +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13)
  )

# Print the plot 
print(sensitivity_3months_surv_plot_HF_cox)



# Extract and round all relevant columns to nearest 5
risk_table_long_HF_cox_unadjusted_sensitivity_3months <- sensitivity_3months_surv_plot_HF_cox$table$data %>%
  mutate(across(
    c(n.risk, n.event, n.censor, cum.n.event, cum.n.censor,strata_size,llabels),
    ~ 5 * round(. / 5)
  ))


# Save the rounded table as a CSV file for extraction
output_file_risk_table_long_HF_cox_unadjusted_sensitivity_3months <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Sensitivity analysis_3months/risk_table_long_HF_cox_unadjusted_sensitivity_3months.csv"
write.csv(risk_table_long_HF_cox_unadjusted_sensitivity_3months, file = output_file_risk_table_long_HF_cox_unadjusted_sensitivity_3months, row.names = FALSE)












##########################################################################################################################
## Plot Adjusted survival curve for Age and Sex #################################################################
##########################################################################################################################


# Define a new dataset for adjusting survival curve, adjusting for age and sex
sensitivity_3months_surv_HF_cox_new_data <- data.frame(
  rehab_indicator = c(0, 1),          # No rehab vs Rehab
  age_at_tavi = mean(clean_cohort$age_at_tavi, na.rm = TRUE),  # Adjust for mean age
  sex_code = factor(c("Female","Male"), levels = levels(clean_cohort$sex_code)),   # Adjust sex
  # Include other covariates if necessary
  stringsAsFactors = FALSE
)

# Fit the survival model and compute predicted survival
sensitivity_3months_surv_object_HF_cox_adj <- survfit(
  sensitivity_3months_cox_model_CR_HF_primary_outcome, 
  newdata = sensitivity_3months_surv_HF_cox_new_data
)

# Plot the adjusted survival curves
plot(
  sensitivity_3months_surv_object_HF_cox_adj,
  col = c("red", "blue"),            # Colors for No Rehab and Rehab
  lwd = 2,                           # Line width
  lty = 1,                           # Line type
  mark.time = TRUE,                  # Add censoring marks
  main = "Sensitivity Analysis: Survival Curve for HF Rehospitalization (Adjusted for Age and Sex)",
  xlab = "Calendar Time (days), starting from days since TAVI discharge plus 90 days and ending at end of followup (31 March 2024)", 
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

sensitivity_3months_surv_object_HF_cox_adj_plot <- ggsurvplot(
  sensitivity_3months_surv_object_HF_cox_adj,
  data = clean_cohort,
  conf.int = TRUE,  # Show confidence intervals
  pval = FALSE,      
  risk.table = FALSE, 
  legend.title = "Rehabilitation Exposure",
  legend.labs = c("Not Exposed to Rehab", "Exposed to Rehab"),
  palette = c("#E69F00", "#0072B2"),  # Distinct colors for rehab vs non-rehab
  linetype = c("solid", "dashed"),  # Solid for rehab, dashed for non-rehab
  break.time.by = 90,  # Adjust x-axis breakpoints for readability
  title = "Sensitivity Analysis 90 days post discharge: Survival Curve for Heart Failure Rehospitalization (Adjusted for Age and Sex)",
  xlab = "Calendar Time, starting from 90 days post TAVI discharge to end of follow-up (31 March 2024)",
  ylab = "Survival Probability",
  ggtheme = theme_minimal()
)



sensitivity_3months_surv_object_HF_cox_adj_plot$plot <- sensitivity_3months_surv_object_HF_cox_adj_plot$plot +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13)
  )

# Print the plot 
print(sensitivity_3months_surv_object_HF_cox_adj_plot)





##########################################################################################################################
##Adjusted model: All Covariates  ####################################################
##########################################################################################################################

#Covariate selection by study cardiologist Professor Tom Marwick
#   Complications after TAVI likely dependent on: age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
#  BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
#  previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications  


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
# Cox Proportional Hazards Model  ########
#############################################################
#Fit the Cox proportional hazards model with the coxph function


# Fit time-varying Cox model with calendar time as the timescale, with all covariates
sensitivity_3months_cox_model_HF_primary_outcome_model_adjusted <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + SMOKING_STATUS + 
    BMI_category + NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + KATZ_INDEX_CAT +  binary_known_poor_mobility + binary_Known_DIABETES + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + LV_FUNCTION + MITRAL_REGURGITATION + TAVI_procedural_complications + year_tavi_discharge  # Include year
  ,
  data = clean_cohort
)

# Summary of the Cox model
summary(sensitivity_3months_cox_model_HF_primary_outcome_model_adjusted)



# Check Assumptions
# Test proportional hazards assumption
cox.zph(sensitivity_3months_cox_model_HF_primary_outcome_model_adjusted)


###NOTES: Model runs with all covariates. However, proportional hazards assumptions test does NOT hold
## NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE; binary_known_poor_mobility and TAVI_procedural_complications is less than 0.05
##therefore proportional hazards assumptions do not hold for these. 
##Stratify by these variables to correct for hazrad assumptions 





#######################################################################################################################################################################################################
# Fit Cox Proportional Hazards Model with stratification 
#######################################################################################################################################################################################################
#Fit the Cox proportional hazards model with the coxph function
## stratify NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE; binary_known_poor_mobility and TAVI_procedural_complications



# Fit Cox model with calendar time as the timescale, with all covariates, also included REGION in this model
sensitivity_3months_cox_model_HF_primary_outcome_model_adjusted <- coxph(
  Surv(hf_rehosp_censor_time, hf_rehosp_indicator) ~ rehab_indicator + age_at_tavi + sex_code  + IMD_2019_QUINTILES + SMOKING_STATUS + region_name +
    BMI_category + KATZ_INDEX_CAT +  COMORBIDITY_pulmonary_liver_neuro_excardiacvascu +
    previous_MI_time_grouped + binary_known_prev_cardiac_surgery + year_tavi_discharge +  # Include year
    strata(NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE , binary_known_poor_mobility, TAVI_procedural_complications) +        # stratify hazard assumption violating variable 
    binary_Known_DIABETES + ethnicity_5_group + LV_FUNCTION + MITRAL_REGURGITATION
  ,
  data = clean_cohort
)

# Summary of the Cox model
summary(sensitivity_3months_cox_model_HF_primary_outcome_model_adjusted)
##Interpretation: in the sensitivity analysis, after adjustment, there is no association between CR exposure and HF rehosps beyond 90days. Before adjustment there is an association 


# Test proportional hazards assumption
cox.zph(sensitivity_3months_cox_model_HF_primary_outcome_model_adjusted)



###NOTES: Proportional hazard assumption test holds true for Global and holds true for all other variables 

#NOTES: P-values are greater than 0.05 for Global, hold true for all other variables







#######################################################################################################################################################################################################
# Extract model summary Model 
#######################################################################################################################################################################################################


# NOTES INTERPRETATION OF MAIN FINDING: Adjusting for all key covariates, accounting for baseline calendar time CR, starting from days since tavi discharge plus 90days, exposure is not associated with the rate of HF rehosp 

# Extract model summary

#Broom package extraction 
library(broom)
library(dplyr)
library(readr)

# Tidy the Cox model with HR and CIs
cox_results_sensitivity_3months_cox_model_HF_primary_outcome_model_fully_adjusted <- tidy(
  sensitivity_3months_cox_model_HF_primary_outcome_model_adjusted,
  exponentiate = TRUE,  # gives HR instead of log(HR)
  conf.int = TRUE
)

# Optional: clean 
cox_results_sensitivity_3months_cox_model_HF_primary_outcome_model_fully_adjusted <- cox_results_sensitivity_3months_cox_model_HF_primary_outcome_model_fully_adjusted %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(across(c(HR, lower_CI, upper_CI, std.error, statistic, p_value), round, 3))


# Save the rounded table as a CSV file for extraction
output_file_cox_results_sensitivity_3months_cox_model_HF_primary_outcome_model_fully_adjusted <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Sensitivity analysis_3months/cox_results_sensitivity_3months_cox_model_HF_primary_outcome_model_fully_adjusted.csv"
write.csv(cox_results_sensitivity_3months_cox_model_HF_primary_outcome_model_fully_adjusted, file = output_file_cox_results_sensitivity_3months_cox_model_HF_primary_outcome_model_fully_adjusted, row.names = FALSE)





#########################################################################################################################################
#####                            SUMMARY OF HF REHOSP VOLUMES - HF rehosps with 90day definition                          ###############
#########################################################################################################################################



##NOTE: this includes the filtered out patients from above - reset the data if you do not want the filtered baseline cohort. 


##raw volumes of HF hospitalizations of entire cohort since TAVI + 90 days until study end
## hospitilisations are only HF rehosps post 90days of TAVI discharge

HF_rehosp_raw_vols_90days <- clean_cohort |>
  mutate(tavi_cips_disdate_plus90_days = tavi_cips_disdate + days(90),  
         rehosp_indicator_HF = ifelse((out_readm_hf_date >= tavi_cips_disdate_plus90_days), 1,  0) 
  )
HF_rehosp_raw_vols_90days$rehosp_indicator_HF[is.na(HF_rehosp_raw_vols_90days$rehosp_indicator_HF)] <- 0

HF_rehosp_raw_vols_summary_post_90days <- HF_rehosp_raw_vols_90days |>
  group_by(rehab_indicator) |>    # Group by rehab exposure
  summarise(
    total_patients = round(n_distinct(person_id) / 5) * 5,      # Unique patients in each group and Round to nearest 5
    total_HF_rehosp = round(sum(rehosp_indicator_HF) / 5) * 5,  # Total HF rehospitalizations and # Round to nearest 5
    avg_HF_rehosp_per_patient = total_HF_rehosp / total_patients, # Avg HF rehosp per patient
  )
print(HF_rehosp_raw_vols_summary_post_90days)  





# Save the rounded table as a CSV file for extraction
output_file_HF_rehosp_raw_vols_summary_post_90days <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Sensitivity analysis_3months/HF_rehosp_raw_vols_summary_post_90days.csv"
write.csv(HF_rehosp_raw_vols_summary_post_90days, file = output_file_HF_rehosp_raw_vols_summary_post_90days, row.names = FALSE)

##note this excludes the volume of patients who had a rehosp within 90 days post TAVI, that's why total patients coloumn is not the full cohort







##########################################################################################################################################################################################

