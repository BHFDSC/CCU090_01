##############################################################################################################################################################
#################################################################   TABLE 1: DESCRIPTIVE STATS    ############################################################
##############################################################################################################################################################

#Install Necessary Packages
if (!requireNamespace("tableone", quietly = TRUE)) install.packages("tableone")
library(tableone)



###Select Variables: Specify the variables for Table 1, separating continuous and categorical variables####
table1_vars <- c(
  "age_at_tavi", "sex_code", "ethnicity_5_group", "IMD_2019_QUINTILES", "region_name", "SMOKING_STATUS", 
  "BMI", "BMI_category","obesity_binary", "NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE", "CCS_ANGINA_STATUS_PRE_PRO_STABLE","CSHA_CLINICAL_FRAILTY_SCALE_SCORE",
  "binary_known_poor_mobility", "binary_Known_DIABETES", "binary_known_on_dialysis", "BINARY_known_SEVERE_LIVER_DISEASE", "binary_known_pulmonary_disease",
  "BINARY_known_HISTORY_OF_NEUROLOGICAL_DISEASE", "KATZ_INDEX_CLEAN", "KATZ_INDEX_CAT", "binary_known_extracardiac_arteriopathy", "COMORBIDITY_pulmonary_liver_neuro_excardiacvascu",
  "EXTENT_OF_CORONARY_VESSEL_DISEASE","LV_FUNCTION",
  "binary_known_prev_mi", "previous_MI_time_grouped", "binary_known_prev_cardiac_surgery", "binary_known_prev_pci","binary_previous_known_CABG", "binary_previous_known_valve_operation",
   "binary_previous_known_other_operation_pericardium",
  "MITRAL_REGURGITATION",
  "PA_SYSTOLIC_PRESSURE_MMHG",
  "BINARY_known_PERMANENT_PACING", "BINARY_known_BLEEDING", "BINARY_known_USE_OF_CARDIOPULMONARY_BYPASS","TAVI_procedural_complications",
  "exp_cardiac_rehab_appts_6m", "exp_cardiac_rehab_attends_6m", "rehab_attendance_groups"
)

table1_cat_vars <- c(
  "sex_code", "ethnicity_5_group", "IMD_2019_QUINTILES", "region_name", "SMOKING_STATUS", 
  "BMI_category","obesity_binary", "NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE", "CCS_ANGINA_STATUS_PRE_PRO_STABLE","CSHA_CLINICAL_FRAILTY_SCALE_SCORE",
  "binary_known_poor_mobility", "binary_Known_DIABETES", "binary_known_on_dialysis", "BINARY_known_SEVERE_LIVER_DISEASE", "binary_known_pulmonary_disease",
  "BINARY_known_HISTORY_OF_NEUROLOGICAL_DISEASE", "KATZ_INDEX_CAT", "binary_known_extracardiac_arteriopathy", "COMORBIDITY_pulmonary_liver_neuro_excardiacvascu",
  "EXTENT_OF_CORONARY_VESSEL_DISEASE","LV_FUNCTION",
  "binary_known_prev_mi", "previous_MI_time_grouped", "binary_known_prev_cardiac_surgery", "binary_known_prev_pci","binary_previous_known_CABG", "binary_previous_known_valve_operation",
  "binary_previous_known_other_operation_pericardium", 
  "MITRAL_REGURGITATION",
  "BINARY_known_PERMANENT_PACING", "BINARY_known_BLEEDING", "BINARY_known_USE_OF_CARDIOPULMONARY_BYPASS","TAVI_procedural_complications",
  "rehab_attendance_groups"
)




###Generate Table 1: Create the summary table using CreateTableOne and with only showing the "1" level for binary variables.###
table_1 <- CreateTableOne(
  vars = table1_vars,
  strata = "binary_exposed_rehab",
  data = clean_cohort,
  factorVars = table1_cat_vars,
  includeNA = TRUE,
  addOverall = TRUE # Adds a column for all patients combined
)



### Print Table 1 with standardized format
print(table_1, showAllLevels = TRUE)


# Print the Table 1 with only the "1" level for binary variables
print(
  table_1,
  showAllLevels = FALSE, # Only show the specified level of binary variables
  catDigits = 1, # Adjust the number of decimal places for percentages
  contDigits = 2 # Adjust for continuous variables
)





#####################UPDATE table_1 VAIRABLES NAMES################

### Rename variables in the dataset to match the desired labels ## create new dataset so dont change orginal variable names in data set
### make variable Labels that are cleaner and more refined

clean_cohort_renamed_table_1 <- clean_cohort %>%
  rename(
    Age = age_at_tavi,
    Sex = sex_code,
    Ethnicity = ethnicity_5_group,
    "Socioeconomic status quintile" = IMD_2019_QUINTILES,
    Region = region_name,
    "Smoking Status" = SMOKING_STATUS,
    "Body Mass Index (BMI)" = BMI,
    "Body Mass Index Category" = BMI_category,
    Obesity = obesity_binary,
    "NYHA Dyspnea Status (at baseline)" = NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE,
    "CCS Angina Status (at baseline)" = CCS_ANGINA_STATUS_PRE_PRO_STABLE,
    "CSHA Clinical Frailty Score" = CSHA_CLINICAL_FRAILTY_SCALE_SCORE,
    "Poor Mobility" = binary_known_poor_mobility,
    Diabetes = binary_Known_DIABETES,
    "On Dialysis" = binary_known_on_dialysis,
    "Severe Liver Disease" = BINARY_known_SEVERE_LIVER_DISEASE,
    "Pulmonary Disease" = binary_known_pulmonary_disease,
    "History of Neurological Disease" = BINARY_known_HISTORY_OF_NEUROLOGICAL_DISEASE,
    "Katz Index of Independence Daily Living" = KATZ_INDEX_CLEAN, 
    "Katz Index of Independence Daily Living Categories" = KATZ_INDEX_CAT,
    "Extracardiac Arteriopathy" = binary_known_extracardiac_arteriopathy,
    "Comorbidity (pulmonary, liver, neuro, extracardiac arteriopathy)" = COMORBIDITY_pulmonary_liver_neuro_excardiacvascu,
    "Extent of Coronary Vessel Disease" = EXTENT_OF_CORONARY_VESSEL_DISEASE,
    "LV Function" = LV_FUNCTION,
    "Previous MI" = binary_known_prev_mi,
    "Previous MI time" = previous_MI_time_grouped,
    "Previous Cardiac surgery" = binary_known_prev_cardiac_surgery,
    "Previous PCI" = binary_known_prev_pci,
    "Previous CABG"= binary_previous_known_CABG,
    "Previous Valve Operation" = binary_previous_known_valve_operation,
    "Previous Pericardium Opening Operation" = binary_previous_known_other_operation_pericardium,
    "Mitral Regurgitation" = MITRAL_REGURGITATION,
    "PA systolic pressure (mmHg)" =  "PA_SYSTOLIC_PRESSURE_MMHG",
    "Permanent Pacing" = BINARY_known_PERMANENT_PACING,
    Bleeding = BINARY_known_BLEEDING,
    "Vascular Access Site Complications" = BINARY_known_VASCULAR_ACCESS_SITE_COMPLICATIONS,
    "Use of Cardiopulmonary Bypass" = BINARY_known_USE_OF_CARDIOPULMONARY_BYPASS,
    "Acute Kidney Injury" = BINARY_known_ACUTE_KIDNEY_INJURY,
    "TAVI procedural complications" = TAVI_procedural_complications,
    "Cardiac Rehabilitation Appointments Scheduled" = exp_cardiac_rehab_appts_6m,
    "Cardiac Rehabilitation Appointments Attended" = exp_cardiac_rehab_attends_6m,
    "Cardiac Rehabilitation Appointments Attended in categories" = rehab_attendance_groups
)



###Select Variables with new names: Specify the variables for Table 1, separating continuous and categorical variables####
table1_vars_renamed <- c(
  "Age", "Sex",  "Ethnicity", "Socioeconomic status quintile",  "Region", "Smoking Status", 
  "Body Mass Index (BMI)",  "Body Mass Index Category", "Obesity",  "NYHA Dyspnea Status (at baseline)", 
  "CCS Angina Status (at baseline)",  "CSHA Clinical Frailty Score",  "Poor Mobility",  "Diabetes", 
  "On Dialysis", "Severe Liver Disease", "Pulmonary Disease", "History of Neurological Disease", 
  "Katz Index of Independence Daily Living", "Katz Index of Independence Daily Living Categories", "Extracardiac Arteriopathy", "Comorbidity (pulmonary, liver, neuro, extracardiac arteriopathy)",
  "Extent of Coronary Vessel Disease", "LV Function", "Previous MI",  "Previous MI time",  "Previous Cardiac Surgery", 
  "Previous PCI", "Previous CABG", "Previous Valve Operation", "Previous Pericardium Opening Operation", "Mitral Regurgitation",
  "PA systolic pressure (mmHg)",
  "Permanent Pacing", "Bleeding", "Use of Cardiopulmonary Bypass", "TAVI procedural complications",
  "Cardiac Rehabilitation Appointments Scheduled", "Cardiac Rehabilitation Appointments Attended", "Cardiac Rehabilitation Appointments Attended in categories"
)

table1_cat_vars_renamed <- c(
  "Sex",  "Ethnicity", "Socioeconomic status quintile",  "Region", "Smoking Status", 
  "Body Mass Index Category", "Obesity",  "NYHA Dyspnea Status (at baseline)", 
  "CCS Angina Status (at baseline)",  "CSHA Clinical Frailty Score",  "Poor Mobility",  "Diabetes", 
  "On Dialysis", "Severe Liver Disease", "Pulmonary Disease", "History of Neurological Disease", 
  "Katz Index of Independence Daily Living Categories", "Extracardiac Arteriopathy", "Comorbidity (pulmonary, liver, neuro, extracardiac arteriopathy)",
  "Extent of Coronary Vessel Disease", "LV Function", "Previous MI",  "Previous MI time",  "Previous Cardiac Surgery", 
  "Previous PCI", "Previous CABG", "Previous Valve Operation", "Previous Pericardium Opening Operation", "Mitral Regurgitation",
  "Permanent Pacing", "Bleeding", "Use of Cardiopulmonary Bypass",  "TAVI procedural complications", "Cardiac Rehabilitation Appointments Attended in categories"
)



### After renaming, use the renamed dataset to create the table:
### Generate Table 1 with renamed variables
###Generate Table 1: Create the summary table using CreateTableOne and with only showing the "1" level for binary variables.###
table_1_renamed_variables <- CreateTableOne(
  vars = table1_vars_renamed,
  strata = "binary_exposed_rehab",
  data = clean_cohort_renamed_table_1,
  factorVars = table1_cat_vars_renamed,
  includeNA = TRUE,
  addOverall = TRUE # Adds a column for all patients combined
)




###TO EXTRACT OUTPUT FILE and ROUND COUNTS TO NEAREST 5

# Isolated the categorical summary table from the TableOne object
cat_table <- table_1_renamed_variables$CatTable %>%
  purrr::map(\(x) bind_rows(x, .id = "var")) %>%
  bind_rows(.id = "strata") %>%
  mutate(freq = round(freq / 5) * 5) %>% # Rounded counts to the nearest 5
  mutate(freq = ifelse(freq <10, as.character("<10"), as.character(freq))) %>% # Suppressed counts under 10
  mutate(miss = round(miss / 5) * 5) %>% # Rounded missing counts to the nearest 5
  mutate(miss = ifelse(miss <10, as.character("<10"), as.character(miss))) %>% # Suppressed missing counts under 10
  mutate(stat = paste0(freq, " (", round(percent, digits = 2), ")")) %>% # Combining statistics of interest in one column (n[%])
  mutate(miss_stat = paste0(miss, " (", round(p.miss, digits = 2), ")")) %>% # Combining missingness statistics
  select(strata, var, n , level, stat, miss_stat)  # Only selecting essential variables

# Extract the overall n numbers, and convert them in a row-wise format, to join to the overall summary table
new_rows_n <- as.data.frame(matrix(NA, nrow = 0, ncol = 0))
j <- 1
for(i in unique(cat_table$strata)){
  new_rows_n[j,"strata"] <- i
  new_rows_n[j,"var"] <- "n"
  new_rows_n[j,"stat"] <- round(unique(cat_table[which(cat_table$strata == i),"n"]) / 5) * 5
  new_rows_n[j,"miss_stat"] <- NA
  new_rows_n[j,"level"] <- NA
  j <- j+1
}

# Remove the n column from the cat_table as this will be added in a row-wise format later
cat_table <- cat_table %>% 
  select(strata, var, level, stat, miss_stat)


library(purrr)

# Isolating the continuous summary table from the TableOne object
cont_table <- table_1_renamed_variables$ContTable %>%
  map(\(x) as_tibble(x, rownames = "var")) %>%
  bind_rows(.id = "strata") %>%
  mutate(miss = round(miss / 5) * 5) %>% # Rounding missingness counts to the nearest 5
  mutate(miss = ifelse(miss <10, as.character("<10"), as.character(miss))) %>% # Suppressing missing counts under 10
  mutate(miss_stat = paste0(miss, " (", round(p.miss, digits = 2), ")")) %>% # Combining missingness statistics in one column
  mutate(stat = paste0( round(mean, digits = 2), " (", round(sd, digits = 2), ")")) %>% # Combining statistics of interest (mean[SD]) into one column
  select(strata, var, stat, miss_stat) %>% # Selecting variables of interest
  mutate(level = NA) # Adding an empty level column so that format matches to cat_table

# Combing overall count, categorical summary and continuous summary into one table
overall_table1_output <- rbind(new_rows_n, cat_table, cont_table) %>%
  pivot_wider(names_from = strata, 
              values_from = c(stat, miss_stat))

print(overall_table1_output) #save this as a CSV file so can take output


# Save the rounded table as a CSV file for extraction 
output_file_table1 <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/Table_one/overall_table1_output.csv"
write.csv(overall_table1_output, file = output_file_table1, row.names = TRUE)

message("Table 1 with rounded counts has been saved to: ", output_file_table1)



# Print the table 1 with renamed variables
# Print the Table 1 with only the "1" level for binary variables
#NOTES: This is not the rounded values
table_1_renamed_variables <- print(
  table_1_renamed_variables,
  showAllLevels = FALSE, # Only show the specified level of binary variables
  catDigits = 1, # Adjust the number of decimal places for percentages
  contDigits = 2, # Adjust for continuous variables
  noSpaces = TRUE        # Ensure clean column names
) %>%
  as.data.frame()








##OPTION TO EXTRACT WITH BELOW CODE
####################################################################################

# finalfit method with addon functions
# add finalfit package and addon functions created by Jamie
setwd("C:\\ProgramData\\UserDataFolders\\S-1-5-21-1534755828-1230176972-3053122617-1012\\My Files\\Home Folder\\Justin Braver\\Scripts\\Table_one_JN")
if (!requireNamespace("finalfit", quietly = TRUE)) install.packages("finalfit")
library(finalfit)
library(tidyverse)
source("finalfit addon functions.R")


# Can add or change labels ---- Example - only include the labels with spaces for now
clean_cohort_table1 = clean_cohort %>% 
  mutate(
    ethnicity_5_group = ethnicity_5_group %>% ff_label("Ethnicity Grouping (5 levels)"), 
    
  )



# Original finalfit Table 1
table_1_original = clean_cohort_renamed_table_1 %>% 
  summary_factorlist(
    dependent = "binary_exposed_rehab",
    explanatory = table1_vars_renamed_test,
    total_col = TRUE,      # Add column with row count
    add_col_totals = TRUE,  # Add row with column count
  )

table_1_original

# Rounded and redacted Table 1 (optional: binary reference levels removed)
table_1_safe_output = table_1_original %>% 
  ff_round_counts() %>%    # Round counts
  ff_redact_counts() %>%   # Redact low counts
  ff_remove_ref_custom()   # Optional: remove reference levels, useful when many binary variables

table_1_safe_output

