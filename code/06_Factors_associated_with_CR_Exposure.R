########################################################################################################################################################
#####                            OBJECTIVE 3. Factors associated with cardiac  rehabilitation exposure                          ########################
#####        Quantify the factors associated with cardiac rehabilitation attendance - demographic, socioeconomic, and clinical factors  ################
########################################################################################################################################################


########perform Cox regression for time to cardiac rehab exposure while accounting for censoring due to death or the end of follow-up##############

######  Step 1: Install and Load Necessary Packages ######
# Install if not already installed
if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
if (!requireNamespace("survminer", quietly = TRUE)) install.packages("survminer")

# Load packages
library(survival)
library(survminer)  # for better visualization and summaries
library(dplyr)
library(tidyr)


######  Step 2: Prepare the Data ######
  #time variable (cardiac_rehab_follow_up_days) - time from discharge to cardiac rehab, death, or censoring (same variable created during script6 for rehab rate)
  #status variable indicating whether rehab exposure occurred (1 for rehab exposure, 0 for censoring). (numeric variable 'exp_cardiac_rehab_6m')
  #Covariates representing patient demographic, socioeconomic, and clinical factors (e.g., age, sex, IMD_2019_QUINTILES, etc.).

# Create survival object for time to cardiac rehab exposure while accounting for censoring due to death or the end of follow-up
clean_cohort <- clean_cohort %>%
  mutate(
    rehab_status = ifelse(is.na(exp_cardiac_rehab_6m), 0, ifelse(exp_cardiac_rehab_6m == 1, 1, 0)),  # 1 for exposed, 0 for censored
    cardiac_rehab_time = cardiac_rehab_follow_up_days        # Time in days (time to rehab or censoring)
  )



######  Step 3: Run Cox Proportional Hazards Regression  ######
#NOTE: Covariates were selected by study cardiologist Professor Tom Marwick. These covariates are based on the factors that attending CR is likely dependent on. 

#################################################################################################
###########################   MODEL 1   #########################################################
# Cox proportional hazards model for Factors associated with cardiac  rehabilitation exposure ###


cox_model_CR_exposure_model_1 <- coxph(Surv(cardiac_rehab_time, rehab_status) ~ age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + 
                                                       BMI_category + binary_known_poor_mobility + KATZ_INDEX_CAT + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
                                                       NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + SMOKING_STATUS + binary_Known_DIABETES +
                                                       binary_known_prev_cardiac_surgery + TAVI_procedural_complications + region_name,
                   data = clean_cohort)

# View summary of the Cox model
summary(cox_model_CR_exposure_model_1)




###Step 4: Check Proportional Hazards Assumption
# Test for proportional hazards
ph_test_cox_model_CR_exposure_model_1 <- cox.zph(cox_model_CR_exposure_model_1)
print(ph_test_cox_model_CR_exposure_model_1)


###Explanation###
#   cox.zph()  tests whether the proportional hazards assumption holds for each covariate.
#   If the p-value is small (<0.05), it indicates a violation. suggests that the assumption is violated, meaning the effect of a covariate changes over time


### NOTE: system is computationally singular: reciprocal condition number = 2.32943e-18
###MODEL NOT RUNNING CORRECTLY. Baseline hazard rates vary significantly by region, therefore need to run random effects model by stratification by region.



###step 5: Check for Multicollinearity ###
library(car)
vif_values_cox_model_CR_exposure_model_1 <- vif(cox_model_CR_exposure_model_1)
print(vif_values_cox_model_CR_exposure_model_1)

##Look for variables with VIF > 5 or 10, which indicate potential multicollinearity
##ALL VARIABLES have VIF <5





############################################################################################################
###########################   MODEL 2   ####################################################################
# step3: run Cox proportional hazards model for Factors associated with cardiac  rehabilitation exposure ###
#strata(region_name) allows different baseline hazards per region, avoiding direct estimation of a regional effect.

cox_model_CR_exposure_model_2 <- coxph(Surv(cardiac_rehab_time, rehab_status) ~ age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + 
                                                       BMI_category + binary_known_poor_mobility + KATZ_INDEX_CAT + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
                                                       NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + SMOKING_STATUS + binary_Known_DIABETES +
                                                       binary_known_prev_cardiac_surgery + TAVI_procedural_complications + strata(region_name),
                                                     data = clean_cohort)

# View summary of the Cox model
summary(cox_model_CR_exposure_model_2)




######   Step 4: Check Proportional Hazards Assumption  ###### 
# Test for proportional hazards
ph_test_cox_model_CR_exposure_model_2 <- cox.zph(cox_model_CR_exposure_model_2)
print(ph_test_cox_model_CR_exposure_model_2)

###Explanation###
#   cox.zph()  tests whether the proportional hazards assumption holds for each covariate.
#   If the p-value is small (<0.05), it indicates a violation. suggests that the assumption is violated, meaning the effect of a covariate changes over time

##NOTE: for ethnicity_5_group Proportional Hazards Assumption does not hold true.
### Therefore ethnicity is suspected to have different baseline hazard functions, stratifying by ethnicity_5_group allows for separate baseline hazards



###step 5: Check for Multicollinearity ###
library(car)
vif_values_cox_model_CR_exposure_model_2 <- vif(cox_model_CR_exposure_model_2)
print(vif_values_cox_model_CR_exposure_model_2)

##Look for variables with VIF > 5 or 10, which indicate potential multicollinearity

##ALL VARIABLES have VIF <5






############################################################################################
###########################   MODEL 3   ####################################################
###Cox proportional hazards model for Factors associated with cardiac  rehabilitation exposure stratification by region and ethnicity_5_group 
##strata(region_name) and strata(ethnicity_5_group) allows different baseline hazards per region, avoiding direct estimation of a regional and ethnicity_5_group effect.

cox_model_CR_exposure_model_3 <- coxph(Surv(cardiac_rehab_time, rehab_status) ~ age_at_tavi + sex_code + IMD_2019_QUINTILES + 
                                                       BMI_category + binary_known_poor_mobility + KATZ_INDEX_CAT + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
                                                       NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + SMOKING_STATUS + binary_Known_DIABETES +
                                                       binary_known_prev_cardiac_surgery + TAVI_procedural_complications + strata(region_name) + strata (ethnicity_5_group),
                                                     data = clean_cohort)



# View summary of the Cox model with region and ethnicity stratified
summary(cox_model_CR_exposure_model_3)


###Step 4: Check Proportional Hazards Assumption
# Test for proportional hazards
ph_test_cox_model_CR_exposure_model_3 <- cox.zph(cox_model_CR_exposure_model_3)
print(ph_test_cox_model_CR_exposure_model_3)

###Explanation###
#   cox.zph()  tests whether the proportional hazards assumption holds for each covariate.
#   If the p-value is small (<0.05), it indicates a violation. suggests that the assumption is violated, meaning the effect of a covariate changes over time

##NOTES: Proportional Hazards Assumption holds true now that region and ethnicity_5_group are stratified. However, TAVI_procedural_complications does not hold true (p=0.030), 
##therefore for TAVI_procedural_complications, include an Interaction with Time (Time-Varying Covariate) -see model 4 below



####Step5: Check for Multicollinearity: Use car::vif() (Variance Inflation Factor) to detect collinearity:
##VIF measures how much variance of a regression coefficient is inflated due to collinearity with other predictors.
##Look for variables with VIF > 5 or 10, which indicate potential multicollinearity
##If VIF > 10, this suggests high multicollinearity. Drop one of the highly correlated variables.

###Interpreting VIF Values
## VIF < 5  Low collinearity (Generally acceptable)
## VIF 5-10 Moderate collinearity (Consider investigating)
## VIF > 10 High collinearity (Serious issue; may need correction)


library(car)
vif_values_cox_model_CR_exposure_model_3 <- vif(cox_model_CR_exposure_model_3)
print(vif_values_cox_model_CR_exposure_model_3)


###NOTE: ALL VIF VALUES ARE VERY LOW THUS NO multicollinearity observed  
##ALL VARIABLES have VIF <5




####                                            NOTES: REGION and ethnicity_5_group is in strata                                                            ####
### NOTE: included "region_name" and "ethnicity_5_group" with starta to make it random effects modle with region and ethnicity_5_group, given that there is huge variation b/w region and ethnicity_5_group
### do not assume proportional hazards between regions and ethnicity_5_group.
### baseline hazard rates vary significantly by region.
#What Stratification Does:
#   Allows different baseline hazard functions for different regions.
#   Controls for confounding by ensuring that differences in region do not affect the estimated hazard ratios.
#   prevents estimation of a region effect directly (i.e., region_name is not included as a coefficient in the model).







############################################################################################
###########################   MODEL 4   ####################################################
###Cox proportional hazards model for Factors associated with cardiac  rehabilitation exposure stratification by region and ethnicity_5_group AND
## include an Interaction with Time (Time-Varying Covariate) for TAVI_procedural_complications
##strata(region_name) and strata(ethnicity_5_group) allows different baseline hazards per region, avoiding direct estimation of a regional and ethnicity_5_group effect.
##If the effect of TAVI_procedural_complications changes over time, include a time-dependent interaction


clean_cohort$TAVI_procedural_complications_time <- 
  clean_cohort$TAVI_procedural_complications * clean_cohort$cardiac_rehab_time

cox_model_CR_exposure_model_4 <- coxph(Surv(cardiac_rehab_time, rehab_status) ~ 
                                                       age_at_tavi + sex_code + IMD_2019_QUINTILES + 
                                                       BMI_category + binary_known_poor_mobility + KATZ_INDEX_CAT + 
                                                       COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
                                                       NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + SMOKING_STATUS + binary_Known_DIABETES + 
                                                       binary_known_prev_cardiac_surgery + TAVI_procedural_complications + 
                                                       TAVI_procedural_complications_time + strata(region_name) + strata (ethnicity_5_group),
                                                     data = clean_cohort)





# View summary of the Cox model
summary(cox_model_CR_exposure_model_4)


###Step 4: Check Proportional Hazards Assumption
# Test for proportional hazards
ph_test_cox_model_CR_exposure_model_4 <- cox.zph(cox_model_CR_exposure_model_4)
print(ph_test_cox_model_CR_exposure_model_4)



##NOTES: Proportional Hazards Assumption does NOT hold true for global but holds true now that region and ethnicity_5_group are stratified and TAVI_procedural_complications includes an Interaction with Time, 
##NEED TO FIX TAVI_procedural_complications


####Step5: Check for Multicollinearity: Use car::vif() (Variance Inflation Factor) to detect collinearity:
##VIF measures how much variance of a regression coefficient is inflated due to collinearity with other predictors.
##Look for variables with VIF > 5 or 10, which indicate potential multicollinearity
##If VIF > 10, this suggests high multicollinearity. Drop one of the highly correlated variables.

###Interpreting VIF Values
## VIF < 5  Low collinearity (Generally acceptable)
## VIF 5-10 Moderate collinearity (Consider investigating)
## VIF > 10 High collinearity (Serious issue; may need correction)


library(car)
vif_values_cox_model_CR_exposure_model_4 <- vif(cox_model_CR_exposure_model_4)
print(vif_values_cox_model_CR_exposure_model_4)


###NOTE: ALL VIF VALUES ARE VERY LOW THUS NO multicollinearity observed  
##ALL VARIABLES have VIF <5 excpet TAVI_complications with time effect







############################################################################################
###########################   MODEL 5   ####################################################
###Cox proportional hazards model for Factors associated with cardiac  rehabilitation exposure stratification by region and ethnicity_5_group AND
## include a Time-Varying Coefficient Model (TT Interaction) for TAVI_procedural_complications. 
##Instead of treating TAVI_procedural_complications as constant over time, allow its effect to change over time using a time-dependent interaction term. 
##strata(region_name) and strata(ethnicity_5_group) allows different baseline hazards per region, avoiding direct estimation of a regional and ethnicity_5_group effect.


#Create a Time-Interaction Term - log(time)
#This allows the effect of TAVI_procedural_complications to vary over time while avoiding extreme values at early time points.

clean_cohort <- clean_cohort %>%
  mutate(TAVI_procedural_complications_time_interaction = TAVI_procedural_complications * log(cardiac_rehab_time + 1))

cox_model_CR_exposure_model_5 <- coxph(Surv(cardiac_rehab_time, rehab_status) ~ 
                                         age_at_tavi + sex_code + IMD_2019_QUINTILES + 
                                         BMI_category + binary_known_poor_mobility + KATZ_INDEX_CAT + 
                                         COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
                                         NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + SMOKING_STATUS + binary_Known_DIABETES + 
                                         binary_known_prev_cardiac_surgery + TAVI_procedural_complications + 
                                         TAVI_procedural_complications_time_interaction   ## Time-interaction term
                                          + strata(region_name) + strata (ethnicity_5_group),  ##stratify violing PH variables
                                       data = clean_cohort)





# View summary of the Cox model 
summary(cox_model_CR_exposure_model_5)


###Step 4: Check Proportional Hazards Assumption
# Test for proportional hazards
ph_test_cox_model_CR_exposure_model_5 <- cox.zph(cox_model_CR_exposure_model_5)
print(ph_test_cox_model_CR_exposure_model_5)



##NOTES: Proportional Hazards Assumption still does not hold true for global but holds true now that region and ethnicity_5_group are stratified and TAVI_procedural_complications includes an Interaction with Time, 
##NEED TO FIX TAVI_procedural_complications


####Step5: Check for Multicollinearity: Use car::vif() (Variance Inflation Factor) to detect collinearity:
##VIF measures how much variance of a regression coefficient is inflated due to collinearity with other predictors.
##Look for variables with VIF > 5 or 10, which indicate potential multicollinearity
##If VIF > 10, this suggests high multicollinearity. Drop one of the highly correlated variables.

###Interpreting VIF Values
## VIF < 5  Low collinearity (Generally acceptable)
## VIF 5-10 Moderate collinearity (Consider investigating)
## VIF > 10 High collinearity (Serious issue; may need correction)


library(car)
vif_values_cox_model_CR_exposure_model_5 <- vif(cox_model_CR_exposure_model_5)
print(vif_values_cox_model_CR_exposure_model_5)








############################################################################################
###########################   MODEL 6   ####################################################
###Cox proportional hazards model for Factors associated with cardiac  rehabilitation exposure stratification by region and ethnicity_5_group AND
## Use Splines for a Flexible Effect for TAVI_procedural_complications 
## model TAVI_procedural_complications using restricted cubic splines
##strata(region_name) and strata(ethnicity_5_group) allows different baseline hazards per region, avoiding direct estimation of a regional and ethnicity_5_group effect.


library(splines)

cox_model_CR_exposure_model_6 <- coxph(Surv(cardiac_rehab_time, rehab_status) ~ 
                                         age_at_tavi + sex_code + IMD_2019_QUINTILES + 
                                         BMI_category + binary_known_poor_mobility + KATZ_INDEX_CAT + 
                                         COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
                                         NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + SMOKING_STATUS + binary_Known_DIABETES + 
                                         binary_known_prev_cardiac_surgery + ns(TAVI_procedural_complications, df = 3)   ## restricted cubic splines
                                       + strata(region_name) + strata (ethnicity_5_group),  ##stratify violing PH variables
                                       data = clean_cohort)





# View summary of the Cox model 
summary(cox_model_CR_exposure_model_6)


###Step 4: Check Proportional Hazards Assumption
# Test for proportional hazards
ph_test_cox_model_CR_exposure_model_6 <- cox.zph(cox_model_CR_exposure_model_6)
print(ph_test_cox_model_CR_exposure_model_6)



##NOTES: Proportional Hazards Assumption holds true for global but doesnt hold true for TAVI_procedural_complications with spline  
##NEED TO stratify by tavi_compications and wont have HR


####Step5: Check for Multicollinearity: Use car::vif() (Variance Inflation Factor) to detect collinearity:
##VIF measures how much variance of a regression coefficient is inflated due to collinearity with other predictors.
##Look for variables with VIF > 5 or 10, which indicate potential multicollinearity
##If VIF > 10, this suggests high multicollinearity. Drop one of the highly correlated variables.


library(car)
vif_values_cox_model_CR_exposure_model_6 <- vif(cox_model_CR_exposure_model_6)
print(vif_values_cox_model_CR_exposure_model_6)







############################################################################################
###########################   MODEL 7   ####################################################
###Cox proportional hazards model for Factors associated with cardiac  rehabilitation exposure stratification by region and ethnicity_5_group AND
## stratify by TAVI_procedural_complications given all the above still doesnt not meet assumptions for PH assumption test
##strata(region_name) and strata(ethnicity_5_group) and strata(tavi_complications) allows different baseline hazards , avoiding direct estimation of a regional and ethnicity_5_group effect and tavi_complications.



cox_model_CR_exposure_model_7 <- coxph(Surv(cardiac_rehab_time, rehab_status) ~ 
                                         age_at_tavi + sex_code + IMD_2019_QUINTILES + 
                                         BMI_category + binary_known_poor_mobility + KATZ_INDEX_CAT + 
                                         COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
                                         NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + SMOKING_STATUS + binary_Known_DIABETES + 
                                         binary_known_prev_cardiac_surgery + 
                                         strata(TAVI_procedural_complications) + strata(region_name) + strata (ethnicity_5_group),  ##stratify violating PH variables
                                       data = clean_cohort)





# View summary of the Cox model 
summary(cox_model_CR_exposure_model_7)


###Step 4: Check Proportional Hazards Assumption
# Test for proportional hazards
ph_test_cox_model_CR_exposure_model_7 <- cox.zph(cox_model_CR_exposure_model_7)
print(ph_test_cox_model_CR_exposure_model_7)



##NOTES: Proportional Hazards Assumption holds true for global and holds true for all the other variables 


####Step5: Check for Multicollinearity: Use car::vif() (Variance Inflation Factor) to detect collinearity:
##VIF measures how much variance of a regression coefficient is inflated due to collinearity with other predictors.
##Look for variables with VIF > 5 or 10, which indicate potential multicollinearity
##If VIF > 10, this suggests high multicollinearity. Drop one of the highly correlated variables.


library(car)
vif_values_cox_model_CR_exposure_model_7 <- vif(cox_model_CR_exposure_model_7)
print(vif_values_cox_model_CR_exposure_model_7)

##No evidence of multicollinearity.


##Interpretation: Modle 7 shows the varaibles that are associated with attending CR, namely: age, sex (females less likley), IMD, BMI, katz, cOMORBIDITY_pulmonary_liver_neuro_excardiacvascu, NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLe
##(seems like more stable and more healthy is less likely to attend?? - check each comparison level to make proper analysis)




###Extract results of model 7 and save as CSV file

# Extract summary
cox_model_CR_exposure_model_7_summary <- summary(cox_model_CR_exposure_model_7)

# Extract coefficients and relevant statistics
results_df_cox_model_CR_exposure_model_7 <- as.data.frame(cox_model_CR_exposure_model_7_summary$coefficients)

# Add HR and confidence intervals
confint_df_cox_model_CR_exposure_model_7 <- as.data.frame(cox_model_CR_exposure_model_7_summary$conf.int)
results_df_cox_model_CR_exposure_model_7$HR <- confint_df_cox_model_CR_exposure_model_7$`exp(coef)`
results_df_cox_model_CR_exposure_model_7$lower_CI <- confint_df_cox_model_CR_exposure_model_7$`lower .95`
results_df_cox_model_CR_exposure_model_7$upper_CI <- confint_df_cox_model_CR_exposure_model_7$`upper .95`

# Get stratification variables
strata_vars <- attr(cox_model_CR_exposure_model_7$terms, "specials")$strata
strata_names <- attr(cox_model_CR_exposure_model_7$terms, "term.labels")[strata_vars]

# Create metadata frame
meta_df <- data.frame(
  Variable = "Stratification variables",
  Value = paste(strata_names, collapse = ", ")
)

# Write both results and metadata to a CSV (append metadata on top)
write.table(meta_df, "cox_model_results_with_strata_predictors_of_CR.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
write.table(results_df_cox_model_CR_exposure_model_7, "cox_model_results_with_strata_predictors_of_CR.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE, append = TRUE)




###do the same extraction but with Broom package 
library(broom)
library(dplyr)
library(readr)

# Tidy the Cox model
tidy_results_cox_model_CR_exposure_model_7 <- tidy(cox_model_CR_exposure_model_7, exponentiate = TRUE, conf.int = TRUE)

# Rename columns (optional)
tidy_results_cox_model_CR_exposure_model_7 <- tidy_results_cox_model_CR_exposure_model_7 %>%
  rename(
    HR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  )

# Extract stratification variable names
strata_vars <- attr(cox_model_CR_exposure_model_7$terms, "specials")$strata
strata_names <- attr(cox_model_CR_exposure_model_7$terms, "term.labels")[strata_vars]
strata_string <- paste(strata_names, collapse = ", ")

# Add a metadata column 
tidy_results_cox_model_CR_exposure_model_7 <- tidy_results_cox_model_CR_exposure_model_7 %>%
  mutate(strata_info = NA_character_)

# Add the strata info as a separate row
strata_row <- tibble(
  term = "Stratification variables",
  HR = NA_real_,
  std.error = NA_real_,
  statistic = NA_real_,
  p_value = NA_real_,
  lower_CI = NA_real_,
  upper_CI = NA_real_,
  strata_info = strata_string
)

# Combine and export
output_df_cox_model_CR_exposure_model_7 <- bind_rows(strata_row, tidy_results_cox_model_CR_exposure_model_7)
write_csv(output_df_cox_model_CR_exposure_model_7, "cox_model_results_broom_CR_exposure_model_7.csv")


# Save the rounded table as a CSV file for extraction
output_file_cox_model_CR_exposure_model_7 <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Predictors of rehab/output_df_cox_model_CR_exposure_model_7.csv"
write.csv(output_df_cox_model_CR_exposure_model_7, file = output_file_cox_model_CR_exposure_model_7, row.names = FALSE)






##################################################################################################################################################################
#######################################                     Logistic Regression Model          ###################################################################
### Factors associated with cardiac  rehabilitation exposure - a simple binary analysis (exposed or not exposed to rehab within the follow-up period)#############
##################################################################################################################################################################

############################################################################################
###########################   Logistic Regression Model  1   ###############################
############################################################################################



# Fit the logistic regression model
logistic_model_CR_exposure_1 <- glm(rehab_status ~ 
                                      age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + region_name +
                                      BMI_category + binary_known_poor_mobility + KATZ_INDEX_CAT + COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
                                      NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + SMOKING_STATUS + binary_Known_DIABETES +
                                      binary_known_prev_cardiac_surgery + TAVI_procedural_complications, 
                      data = clean_cohort, family = binomial)

# Summarize the model results
summary(logistic_model_CR_exposure_1)


###Check for Model Assumptions
###For Logistic regression: Check for multicollinearity using the Variance Inflation Factor (VIF):
library(car)
vif(logistic_model_CR_exposure_1)


###NOTES ON OUTPUTS: coeffeicients are too small, if remove region_name it corrects the output, need to do random effects for region_name

# summary(logistic_model_CR_exposure) Output:
  
#  Coefficients Table: Each covariate has an estimate (log-odds), standard error (SE), z-value, and p-value.
# Interpretation:
#   Positive coefficients indicate an increased probability of cardiac rehab exposure.
#   Negative coefficients suggest a decreased probability.
#   Significance levels (p-value < 0.05) indicate statistically significant relationships.
#   Null deviance and Residual deviance: These indicate model fit. A large reduction from null deviance to residual deviance suggests good explanatory power.
#   Variance Inflation Factor (VIF) Output (vif(logistic_model_CR_exposure)):
  
#   VIF values indicate multicollinearity:
#   VIF > 5 suggests moderate collinearity.
#   VIF > 10 indicates high collinearity that may affect the stability of coefficient estimates.


# All VIF values are below 5, 


#Factors associated with cardiac  rehabilitation exposure - a simple binary analysis (exposed or not exposed to rehab within the follow-up period)
############################################################################################
###########################    Logistic Regression Model  2   ##############################
############################################################################################


# Fit the logistic regression model, Use region_name as a Random Effect (Mixed-Effects Model)
#Regional variability is important, so use a Generalized Linear Mixed Model (GLMM) instead
#controls for region without estimating a coefficient for each region
#region affects baseline rehab likelihood
#Use if regional differences are important but coefficients are unstable

install.packages("Matrix", type = "binary")
install.packages("lme4")
library(lme4)
library(Matrix)


logistic_model_CR_exposure_2 <- glmer(rehab_status ~ 
                                        age_at_tavi + sex_code + ethnicity_5_group + IMD_2019_QUINTILES + 
                                        BMI_category + binary_known_poor_mobility + KATZ_INDEX_CAT + 
                                        COMORBIDITY_pulmonary_liver_neuro_excardiacvascu + 
                                        NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE + SMOKING_STATUS + 
                                        binary_Known_DIABETES + binary_known_prev_cardiac_surgery + 
                                        TAVI_procedural_complications + (1 | region_name),  # Random effect
                                      data = clean_cohort, family = binomial,nAGQ = 0)  # Faster computation

# Summarize the model results
summary(logistic_model_CR_exposure_2)

###Check for Model Assumptions
###For Logistic regression: Check for multicollinearity using the Variance Inflation Factor (VIF):
library(car)
vif(logistic_model_CR_exposure_2)

## OUTPUT = All VIF values are below 5, therefore no further changes required


##NOTES: MODEL RUNS WELL - key predictors are identified



###Extract results and save as CSV file
install.packages("broom.mixed", type = "binary")
install.packages("future", type = "binary")
library(future)
library(broom.mixed)
library(dplyr)
library(readr)

# Fixed effects with ORs and CIs
fixed_effects_logistic_model_CR_exposure_2 <- tidy(logistic_model_CR_exposure_2, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
  rename(
    OR = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    p_value = p.value
  ) %>%
  mutate(effect_type = "fixed")

# Random effects (raw values)
random_effects_logistic_model_CR_exposure_2 <- tidy(logistic_model_CR_exposure_2, effects = "ran_vals") %>%
  mutate(
    OR = exp(estimate),
    lower_CI = NA,
    upper_CI = NA,
    p_value = NA,
    std.error = NA,
    statistic = NA,
    effect_type = "random"
  ) %>%
  select(term, group, level, OR, lower_CI, upper_CI, p_value, effect_type)

# Combine
combined_logistic_model_CR_exposure_2 <- bind_rows(
  fixed_effects_logistic_model_CR_exposure_2 %>% mutate(group = NA, level = NA) %>% select(term, group, level, OR, lower_CI, upper_CI, p_value, effect_type),
  random_effects_logistic_model_CR_exposure_2
)

# Export
write_csv(combined_logistic_model_CR_exposure_2, "logistic_model_results_with_random_effects_predictors_of_CR.csv")

# Save the rounded table as a CSV file for extraction
output_file_logistic_model_CR_exposure_2 <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Predictors of rehab/logistic_model_results_with_random_effects_predictors_of_CR.csv"
write.csv(combined_logistic_model_CR_exposure_2, file = output_file_logistic_model_CR_exposure_2, row.names = FALSE)




############################################################################################################################################################################








