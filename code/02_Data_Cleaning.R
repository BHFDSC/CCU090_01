
####Step 1: load the data##
cohort_dat

#new clean data set
clean_cohort <- cohort_dat



####Step 2: rename variables with numbers in their names####

names(clean_cohort)

# Identify and rename only the variables with numbers in their names ##$REPEAT THIS TWICE
names(clean_cohort) <- gsub("^\\d+_", "", names(clean_cohort))

# Verify the updated column names
names(clean_cohort)


###NEED TO REPEAT STEP 2 TWICE TO CLEAN THE VARIABLES


##############################################################################
#########Step 3:	Remove the numbers from the character strings AND###########
#########Step 4:	Convert character variables to factors######################
##############################################################################

# Load dplyr
library(dplyr)

structure(clean_cohort)
summary(clean_cohort)

##Clean character values by removing the number in front of each word and convert to factors


weird_columns <- c("DIABETES", "SMOKING_STATUS", "ON_DIALYSIS", "PREV_MI_AND_INT_BTW_PRO_AND_LAST_MI", "HISTORY_OF_PULMONARY_DISEASE", "SEVERE_LIVER_DISEASE",
                   "HISTORY_OF_NEUROLOGICAL_DISEASE", "EXTRACARDIAC_ARTERIOPATHY", "POOR_MOBILITY", "EXT_CALCIFICATION_OF_ASCENDIN_AORTA", "PRE_OPERATIVE_HEART_RHYTHM",
                   "PREVIOUS_CARDIAC_SURGERY", "PREVIOUS_TAVI", "PREVIOUS_PCI", "CRITICAL_PRE_OPERATIVE_STATUS_V4", "CCS_ANGINA_STATUS_PRE_PRO_STABLE", "NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE",
                  "CSHA_CLINICAL_FRAILTY_SCALE_SCORE", "PA_SYSTOLIC_PRESSURE_MEASURED", "MITRAL_REGURGITATION", "LV_FUNCTION", "EXTENT_OF_CORONARY_VESSEL_DISEASE",
                  "LEFT_MAIN_STEM_DISEASE", "PROCTORED_CASE","PER_PROCEDURAL_IMAGING", "PROCEDURE_URGENCY", "ANAESTHESIA_INTENDED_TREATMENT", "UNPLAN_CONVERSION_GEN_ANAESTHETIC",
                  "CEREBRAL_CIRCULATION_PROTEC_DEVICE", "DELIVERY_APPROACH", "USE_OF_CARDIOPULMONARY_BYPASS", "DEVICE_FAILURE_REFERS_TO_VALVE_ONLY","VASCULAR_CLOSURE_TECHNIQUE",
                  "VALVE_SUCCESSFULLY_DEPLOYED", "VALVE_NOT_DEPLOYED_SUCCESS_REASON", "VALVE_MALPOSITIONING", "BAIL_OUT_VALVE_IN_VALVE","POST_IMPLA_BAL_DILATATION_OF_VALVE",
                  "FURTHER_VALVE_INTERVEN_B4_DISCHARG", "TAMPONADE_DURING_POST_PROCEDURE", "MAJ_APICAL_CANNULATION_COMPLICATION", "CONV_TO_FULL_STERNOTOMY_DURING_PRO", "BAILOUT_PCI", 
                 "PERI_PROCEDURAL_MI_72HRS_AFTER_PRO", "PERMANENT_PACING", "CVA_UP_TO_DX","VASCULAR_ACCESS_SITE_COMPLICATIONS",
                 "PERCUTANEOUS_CLOSURE_DEVICE_FAIL", "BLEEDING", "ACUTE_KIDNEY_INJURY_WITHIN_7_DAYS", "NEW_RENAL_REPLACEMENT_THERAPY", "DISCHARGE_DESTINATION_CARDIO_WARD",
                 "DRUGS_AT_DISCHARGE_ANTITHROMBOTIC", "DRUGS_AT_DISCHARGE_ANTI_PLATELET", "IF_ALIVE_CCS_ANGINA_STATUS_1Y", "IF_ALIVE_NYHA_DYSPNOEA_STATUS_1Y", "IF_ALIVE_CCS_ANGINA_STATUS_3Y",
               "IF_ALIVE_NYHA_DYSPNOEA_STATUS_3Y", "LATE_VALVE_STENOSIS", "LATE_INTRINSIC_VALVE_REGURGITATN", "VALVE_FAILURE_MODE", "LATE_PARAVALVULAR_REGURGITATION",
               "INTERVENTION_PARAVALVULA_REG")


# Process the columns: remove first 3 characters and convert to factor
clean_cohort <- clean_cohort %>%
  mutate(across(all_of(weird_columns), ~ as.factor(substr(., 4, nchar(.)))))


# Check structure of the processed columns
str(clean_cohort[weird_columns])

# Verify the levels for one/two of the processed columns
levels(clean_cohort$DIABETES)
levels((clean_cohort$SEVERE_LIVER_DISEASE))
levels((clean_cohort$SMOKING_STATUS))


#######################################################################################################
#################### Re-code and Convert other variables to correct types###############################
######################################################################################################

unique(clean_cohort$sex_code)
#recode sex_code
clean_cohort$sex_code <- ifelse(clean_cohort$sex_code == 2, 0, 1)
table(clean_cohort$sex_code)
#convert to factor
clean_cohort$sex_code <- factor(clean_cohort$sex_code, levels = c(0, 1), labels = c("Female", "Male"))
summary(clean_cohort$sex_code)

#Check missing variables
colSums(is.na(clean_cohort))

# Count missing values in sex_code
sum(is.na(clean_cohort$sex_code))

##re-code ethnicity_5_group to a factor
class(clean_cohort$ethnicity_5_group)
unique(clean_cohort$ethnicity_5_group)

### Convert ethnicity_5_group into a factor###
clean_cohort$ethnicity_5_group <- factor(clean_cohort$ethnicity_5_group)
class(clean_cohort$ethnicity_5_group)

# Check the levels of the factor
levels(clean_cohort$ethnicity_5_group)


# Convert to factor and specify levels in desired order
clean_cohort$ethnicity_5_group <- factor(
  clean_cohort$ethnicity_5_group,
  levels = c(
    "White", 
    "Mixed or multiple ethnic groups", 
    "Asian or Asian British", 
    "Black, Black British, Caribbean or African", 
    "Other ethnic group"
  )
)

# Verify the levels
levels(clean_cohort$ethnicity_5_group)

# Count missing values in ethnicity_5_group
sum(is.na(clean_cohort$ethnicity_5_group))


#Summary ethnicity_5_group after conversion
# Summarize the factor variable
summary(clean_cohort$ethnicity_5_group)



###convert ethncity_18_group variable to factor###
clean_cohort$ethnicity_18_group <- as.factor(clean_cohort$ethnicity_18_group)
levels(clean_cohort$ethnicity_18_group)
class(clean_cohort$ethnicity_18_group)
sum(is.na(clean_cohort$ethnicity_18_group))

#Summary ethnicity_18_group after conversion
# Summarize the factor variable
summary(clean_cohort$ethnicity_18_group)

levels(clean_cohort$ethnicity_18_group)
structure((clean_cohort$ethnicity_18_group))


###CONVERT REGION TO A FACTOR###
clean_cohort <- clean_cohort %>%
  mutate(region_name = as.factor(region_name))

class(clean_cohort$region_name)
unique(clean_cohort$region_name)
summary(clean_cohort$region_name)



######CONVERT IMD_2019_DECILES and IMD_2019_QUINTILES in ORDERED FACTORS########

class(clean_cohort$IMD_2019_DECILES)
class(clean_cohort$IMD_2019_QUINTILES)


#### Convert IMD_2019_DECILES to an ordered factor#####################
###Deciles range from 1 (most deprived 10%) to 10 (least deprived 10%)###

clean_cohort <- clean_cohort %>%
  mutate(IMD_2019_DECILES = factor(
    IMD_2019_DECILES,
    levels = 1:10,       # Explicitly specify the order
    ordered = TRUE       # Make it an ordered factor
  ))



#### Convert IMD_2019_QUINTILES to an ordered factor######################
###Quintiles range from 1 (most deprived 20%) to 5 (least deprived 20%)####
clean_cohort <- clean_cohort %>%
  mutate(IMD_2019_QUINTILES = factor(
    IMD_2019_QUINTILES,
    levels = 1:5,       # Explicitly specify the order
    ordered = TRUE      # Make it an ordered factor
  ))

###VERIFY THE RESULTS###
# Check the structure and levels
str(clean_cohort$IMD_2019_DECILES)
levels(clean_cohort$IMD_2019_DECILES)
str(clean_cohort$IMD_2019_QUINTILES)
levels(clean_cohort$IMD_2019_QUINTILES)

# Summarize the data
table(clean_cohort$IMD_2019_DECILES)
table(clean_cohort$IMD_2019_QUINTILES)




####STEP 5. Grouping levels to Create binary variables##########

#binary known diabetes (missing values assigned 0)
summary(clean_cohort$DIABETES)

clean_cohort <- clean_cohort %>%
  mutate(binary_Known_DIABETES = if_else(
    DIABETES %in% c("Diabetes (dietary control)", "Diabetes (insulin)", 
                    "Diabetes (oral medicine)", "Newly diagnosed diabetes"), 1, 
    0, missing = 0
  ))

table(clean_cohort$binary_Known_DIABETES)
head(clean_cohort$binary_Known_DIABETES)


##REPEAT FOR OTHER VARIABLES: CONVERT TO BINARY VARABLES###



# Convert ON_DIALYSIS into a binary variable (missing values assigned 0)
summary(clean_cohort$ON_DIALYSIS)

clean_cohort <- clean_cohort %>%
  mutate(binary_known_on_dialysis = if_else(
    ON_DIALYSIS == "Yes", 1, 0, missing = 0
  ))


# Summarize the new binary variable
table(clean_cohort$binary_known_on_dialysis)


#Binary Previous known MI (missing values assigned 0)
summary(clean_cohort$PREV_MI_AND_INT_BTW_PRO_AND_LAST_MI)

clean_cohort <- clean_cohort %>%
  mutate(binary_known_prev_mi = if_else(
    PREV_MI_AND_INT_BTW_PRO_AND_LAST_MI %in% c("MI < 6 hours", "MI > 90 days", 
                                               "MI 1-30 days", "MI 31-90 days", 
                                               "MI 6-24 hours"), 1, 0, missing = 0))
table(clean_cohort$binary_known_prev_mi)


#Binary History of pulmonary disease (missing values assigned 0)
summary(clean_cohort$HISTORY_OF_PULMONARY_DISEASE)

clean_cohort <- clean_cohort %>%
  mutate(binary_known_pulmonary_disease = if_else(
    HISTORY_OF_PULMONARY_DISEASE %in% c("Asthma", "COAD/emphysema", 
                                        "Other significant pulmonary disease"), 1, 0, missing = 0))

table(clean_cohort$binary_known_pulmonary_disease)


#Binary poor mobility (missing values assigned 0)
summary(clean_cohort$POOR_MOBILITY)

clean_cohort <- clean_cohort %>%
  mutate(binary_known_poor_mobility = if_else(
    POOR_MOBILITY == "Yes", 1, 0, missing = 0))

table(clean_cohort$binary_known_poor_mobility)


#Binary Previous known PCI (missing values assigned 0)
summary(clean_cohort$PREVIOUS_PCI)

clean_cohort <- clean_cohort %>%
  mutate(binary_known_prev_pci = if_else(
    PREVIOUS_PCI %in% c("Yes - as part of a staged or hybrid procedure", 
                        "Yes - previous standalone PCI (NOT as part of staged or hybrid procedure)"), 1, 0, missing = 0))

table(clean_cohort$binary_known_prev_pci)



############################################################################
##########STEP 6: Cleaning PREVIOUS CARDIAC SURGERY VARIABLE#####################
###########################################################################

#####EXPLANATION OF THE CODE ################
#Edge Cases (e.g., "No;1. Previous CABG"):Accounts for compound levels that imply both "No" and "Previous CABG". These are considered 0 for "previous cardiac surgery".
#If any of the 19 levels have the word "NO" in it, we assigned a 0. 
#missing = 0: Ensures that NA values in PREVIOUS_CARDIAC_SURGERY are assigned 0 in the binary variable.
#0 for "No":Assigns 0 to rows where "No" is the only value.


#Binary Previous known cardiac surgery (missing values assigned 0)
summary(clean_cohort$PREVIOUS_CARDIAC_SURGERY)

clean_cohort <- clean_cohort %>%
  mutate(binary_known_prev_cardiac_surgery = if_else(
    PREVIOUS_CARDIAC_SURGERY %in% c("Previous CABG", 
                                    "Previous valve operation", 
                                    "Previous CABG;2. Previous valve operation", 
                                    "Other operation requiring opening of the pericardium",
                                    "Other operation requiring opening of the pericardium;1. Previous CABG",
                                    "Other operation requiring opening of the pericardium;1. Previous CABG;2. Previous valve operation",
                                    "Other operation requiring opening of the pericardium;2. Previous valve operation",
                                    "Previous CABG;3. Other operation requiring opening of the pericardium",
                                    "Previous valve operation;1. Previous CABG",
                                    "Previous valve operation;3. Other operation requiring opening of the pericardium"
                                    ), 1, 0, missing = 0))

#Check summary of binary known previous cardiac surgery
table(clean_cohort$binary_known_prev_cardiac_surgery)




# Create binary variables for a)Previous CABG and separately b)Previous Valve operation by grouping the relevant levels in the PREVIOUS_CARDIAC_SURGERY for KNOWN type of surgery (missing values assigned 0)

#Edge Cases (e.g., "No;1. Previous CABG"):Accounts for compound levels that imply both "No" and "Previous CABG". These are considered 0 for "previous cardiac surgery".
#If any of the original 19 levels have the word "NO" in it, we assigned a 0.

#Previous KNOWN CABG - includes CABG only and/or any combination of CABG with another cardiac procedure. Included previous cardiac surgery CABG levels are:
#  "Previous CABG", 
#  "Previous CABG;2. Previous valve operation",
#  "Other operation requiring opening of the pericardium;1. Previous CABG",
#  "Other operation requiring opening of the pericardium;1. Previous CABG;2. Previous valve operation",
#  "Previous CABG;3. Other operation requiring opening of the pericardium",
#  "Previous valve operation;1. Previous CABG"


# Create binary variable for previous KNOWN CABG
clean_cohort <- clean_cohort %>%
  mutate(binary_previous_known_CABG = if_else(
    PREVIOUS_CARDIAC_SURGERY %in% c(
      "Previous CABG", 
      "Previous CABG;2. Previous valve operation", 
      "Other operation requiring opening of the pericardium;1. Previous CABG", 
      "Other operation requiring opening of the pericardium;1. Previous CABG;2. Previous valve operation", 
      "Previous CABG;3. Other operation requiring opening of the pericardium",
      "Previous valve operation;1. Previous CABG"
    ), 1, 0, missing = 0))

table(clean_cohort$binary_previous_known_CABG)


#Previous KNOWN VALVE OPERATION - includes VALVE OPERATION only and/or any combination of previous cardiac surgery levels with VALVE OPERATION. Included previous cardiac surgery VALVE OPERATION levels are:
#  "Previous valve operation", 
#  "Previous CABG;2. Previous valve operation",
#  "Other operation requiring opening of the pericardium;2. Previous valve operation",
#  "Other operation requiring opening of the pericardium;1. Previous CABG;2. Previous valve operation",
#  "Previous valve operation;3. Other operation requiring opening of the pericardium",
#  "Previous valve operation;1. Previous CABG"


# Create binary variable for previous VALVE OPERATION
clean_cohort <- clean_cohort %>%
  mutate(binary_previous_known_valve_operation = if_else(
    PREVIOUS_CARDIAC_SURGERY %in% c(
      "Previous valve operation", 
      "Previous CABG;2. Previous valve operation", 
      "Other operation requiring opening of the pericardium;2. Previous valve operation", 
      "Other operation requiring opening of the pericardium;1. Previous CABG;2. Previous valve operation", 
      "Previous valve operation;3. Other operation requiring opening of the pericardium", 
      "Previous valve operation;1. Previous CABG"
    ), 1, 0, missing = 0))

table(clean_cohort$binary_previous_known_valve_operation)



# Create binary variable for previous "Other operation requiring opening of the pericardium"
clean_cohort <- clean_cohort %>%
  mutate(binary_previous_known_other_operation_pericardium = if_else(
    PREVIOUS_CARDIAC_SURGERY %in% c(
      "Other operation requiring opening of the pericardium", 
      "Other operation requiring opening of the pericardium;1. Previous CABG", 
      "Other operation requiring opening of the pericardium;2. Previous valve operation", 
      "Other operation requiring opening of the pericardium;1. Previous CABG;2. Previous valve operation",
      "Previous CABG;3. Other operation requiring opening of the pericardium",
      "Previous valve operation;3. Other operation requiring opening of the pericardium"
      ), 1, 0, missing = 0))


table(clean_cohort$binary_previous_known_other_operation_pericardium)





# Create a new PREVIOUS CARDIAC SURGERY variable (CLEAN_PREVIOUS_CARDIAC_SURGERY) with updated levels, reducing the original 19 levels down to 6 levels.
#Missing values are included in "No known previous cardiac surgery"
#Edge Cases (e.g., "No;1. Previous CABG"):Accounts for compound levels that imply both "No" and "Previous CABG". These are assigned to "No known previous cardiac surgery".
#If any of the original 19 levels have the word "NO" in it, we assigned to "No known previous cardiac surgery".

clean_cohort <- clean_cohort %>%
  mutate(CLEAN_PREVIOUS_CARDIAC_SURGERY = case_when(
    # Assign "Previous valve operation"
    PREVIOUS_CARDIAC_SURGERY %in% c( "Previous valve operation" 
                                     ) ~ "Previous valve operation",
    
    # Assign "Previous CABG"
    PREVIOUS_CARDIAC_SURGERY %in% c("Previous CABG" 
                                    ) ~ "Previous CABG",
    
    # Assign "Previous other opening of the pericardium"
    PREVIOUS_CARDIAC_SURGERY %in% c("Other operation requiring opening of the pericardium" 
                                    ) ~ "Previous other opening of the pericardium",
    
    # Assign "Previous CABG and previous valve operation"
    PREVIOUS_CARDIAC_SURGERY %in% c("Previous CABG;2. Previous valve operation",
                                    "Previous valve operation;1. Previous CABG") ~ "Previous CABG and previous valve operation",
    
    # Assign "Previous CABG and Previous other opening of the pericardium" AND/or
    # Assign "Previous valve operation and Previous other opening of the pericardium" AND/or
    # Assign "Previous CABG and previous Valve and previous pericardium operation"
    PREVIOUS_CARDIAC_SURGERY %in% c("Other operation requiring opening of the pericardium;1. Previous CABG",
                                    "Previous CABG;3. Other operation requiring opening of the pericardium",
                                    "Other operation requiring opening of the pericardium;2. Previous valve operation",
                                    "Previous valve operation;3. Other operation requiring opening of the pericardium",
                                    "Other operation requiring opening of the pericardium;1. Previous CABG;2. Previous valve operation"
                                    ) ~ "Previous CABG and/or VALVE and/or previous other opening of the pericardium",
    
    
    # Assign "No known previous cardiac surgery" to all other cases, including "No" and missing values
    TRUE ~ "No known previous cardiac surgery"
  ))


# Convert the new CLEAN_PREVIOUS_CARDIAC_SURGERY variable to a factor with the desired levels
clean_cohort <- clean_cohort %>%
  mutate(CLEAN_PREVIOUS_CARDIAC_SURGERY = factor(
    CLEAN_PREVIOUS_CARDIAC_SURGERY, 
    levels = c("Previous valve operation", 
               "Previous CABG", 
               "Previous other opening of the pericardium",
               "Previous CABG and previous valve operation",
               "Previous CABG and/or VALVE and/or previous other opening of the pericardium",
               "No known previous cardiac surgery")))

summary(clean_cohort$CLEAN_PREVIOUS_CARDIAC_SURGERY)
table(clean_cohort$PREVIOUS_CARDIAC_SURGERY, clean_cohort$CLEAN_PREVIOUS_CARDIAC_SURGERY)



#############################################################
###Step 7: Check for extreme values for HEIGHT and WEIGHT#####
##############################################################

summary(clean_cohort$WEIGHT)
summary(clean_cohort$HEIGHT)

######################################
### Check for extreme HEIGHT values###
######################################
extreme_height <- clean_cohort %>%
  filter(HEIGHT < 1.1 | HEIGHT > 2.5) # thresholds (in meters)


# View the extremes
# Show only HEIGHT & WEIGHT columns for extreme HEIGHT cases
extreme_height %>%
  select(HEIGHT, WEIGHT)

print(extreme_height %>%
        select(HEIGHT, WEIGHT), n=Inf)



###NUMBER EXCLUDEDE DUE TO EXTREME HEIGHT VALUES = 8 people had a height of 0 or around 1.01 etc. - these should be made NAs. There were no to range extreme values only bottom range ie below 1.1
####the remaining extreme values (approx n=445) have units in cm so these just need to be converted to metres 





####################################################################################################
###Step 8:  Convert HEIGHT to meters if in cm and if HEIGHT is less than 1.1 or greater than 2.5 assign NA####
###################################################################################################
library(dplyr)

# Convert HEIGHT to meters if in cm and handle HEIGHT == 0
clean_cohort <- clean_cohort %>%
  mutate(HEIGHT = if_else(
    HEIGHT > 10,          # Assume HEIGHT > 10 indicates cm
    HEIGHT / 100,         # Convert cm to meters
    HEIGHT                # Keep existing values for meters
  )) %>%
  mutate(HEIGHT = if_else(
    HEIGHT < 1.1 | HEIGHT > 2.5,  # Set to NA if HEIGHT is less than 1.1 or greater than 2.5
    NA_real_,                   # Assign NA for invalid values
    HEIGHT                      # Keep existing values
  ))



###VERIFY THE CHANGES###
# Check summary of HEIGHT
summary(clean_cohort$HEIGHT)

# View rows with updated HEIGHT values
print(clean_cohort %>%
  filter(is.na(HEIGHT) | HEIGHT > 2.5 | HEIGHT < 1.1) %>%  # Example range check for extreme values
    select(person_id, HEIGHT, WEIGHT),n=Inf)



#############################################################
###############  Check for extreme WEIGHT values ############
#############################################################

extreme_weight <- clean_cohort %>%
  filter(WEIGHT < 25 | WEIGHT > 210) # thresholds (in KGs)

# View the extremes
# Show only HEIGHT & WEIGHT columns for extreme WEIGHT cases
extreme_weight %>%
  select(HEIGHT, WEIGHT)

print(extreme_weight %>%
        select(HEIGHT, WEIGHT))



###NUMBER EXCLUDEDE DUE TO EXTREME WEIGHT VALUES = ONLY 4 people had an extreme weight
####all 4 weights were below 15kgs ### these people were assigned NA's



#######################################################################
##### If WEIGHT is less than 25kgs or greater than 210kgs assign NA####
#######################################################################


# Update WEIGHT with new conditions
clean_cohort <- clean_cohort %>%
  mutate(WEIGHT = if_else(
    WEIGHT < 25 | WEIGHT > 210,  # Check if WEIGHT is out of bounds
    NA_real_,                   # Assign NA for invalid values
    WEIGHT                      # Keep existing valid values
  ))



###Verify the changes###

# Check summary of WEIGHT
summary(clean_cohort$WEIGHT)

# Inspect rows where WEIGHT is NA
clean_cohort %>%
  filter(is.na(WEIGHT)) %>%
  select(person_id, WEIGHT, HEIGHT)




###########################################################
##################STEP 9: CREATE BMI#######################
###########################################################
# Create BMI variable in the clean_cohort dataset

clean_cohort <- clean_cohort %>%
  mutate(BMI = if_else(
    !is.na(HEIGHT) & !is.na(WEIGHT) & HEIGHT > 0,  # Ensure HEIGHT and WEIGHT are non-missing and HEIGHT > 0
    WEIGHT / (HEIGHT^2),                           # BMI formula: WEIGHT / HEIGHT^2
    NA_real_                                       # Assign NA for invalid cases
  ))


# View the first few rows of the BMI column
head(clean_cohort %>% select(HEIGHT, WEIGHT, BMI))

# Summarize BMI values
summary(clean_cohort$BMI)

# Inspect rows with BMI greater than 60
clean_cohort %>%
  filter(BMI > 60) %>%
  select(person_id, HEIGHT, WEIGHT, BMI)

# Inspect rows with BMI less than 15
clean_cohort %>%
  filter(BMI <15) %>%
  select(person_id, HEIGHT, WEIGHT, BMI)




######For BMI GREATER THAN 90 CONSIDER MAKING NA??############# DISCUSS WITH JOHN??






###################################################################################
########### STEP 10:   CREATE BMI CATEGORY VARIABLE FROM BMI  ###################################
###################################################################################

# Convert BMI into broader categories
clean_cohort <- clean_cohort %>%
  mutate(BMI_category = case_when(
    is.na(BMI) ~ "Missing",                       # Category for missing values
    BMI < 18.5 ~ "Underweight",                   # BMI < 18.5
    BMI >= 18.5 & BMI < 25 ~ "Normal Weight",     # 18.5 ≤ BMI < 25
    BMI >= 25 & BMI < 30 ~ "Overweight",          # 25 ≤ BMI < 30
    BMI >= 30 ~ "Obesity"                         # BMI ≥ 30
  ))

# Convert BMI_category to a factor with ordered levels
clean_cohort <- clean_cohort %>%
  mutate(BMI_category = factor(
    BMI_category,
    levels = c("Underweight", "Normal Weight", "Overweight", "Obesity", "Missing"),
    ordered = TRUE
  ))


###VERIFY RESULTS###
# View the distribution of BMI categories
table(clean_cohort$BMI_category)

# Inspect a sample of rows with BMI and BMI_category
head(clean_cohort %>% select(BMI, BMI_category))



###################################################################################
########### STEP 11:   CREATE BINARY OBESITY VARIABLE FROM BMI  ###########################
###################################################################################

# Create a binary obesity variable
clean_cohort <- clean_cohort %>%
  mutate(obesity_binary = if_else(
    BMI >= 30,  # Obesity threshold (BMI ≥ 30)
    1,          # Assign 1 for obese
    0,          # Assign 0 for non-obese
    missing = NA_real_  # Keep missing values as NA
  ))


###VERIFY RESULTS###
# View the distribution of the obesity_binary variable
table(clean_cohort$obesity_binary, useNA = "ifany")

# Inspect a sample of rows with BMI and obesity_binary
head(clean_cohort %>% select(BMI, obesity_binary))


class(clean_cohort$BMI_category)


# Create a binary 'known' obesity variable with missing variables assigned 0
clean_cohort <- clean_cohort %>%
  mutate(known_obesity_binary = if_else(
    BMI >= 30,  # Obesity threshold (BMI ≥ 30)
    1,          # Assign 1 for obese
    0,          # Assign 0 for non-obese
    missing = 0  # do not keep missing values as NA
  ))

summary(clean_cohort$known_obesity_binary)
table(clean_cohort$known_obesity_binary)



#########################################################################################################################
########### STEP 12: ADDITIONAL CLEANING, REPLACE "Unknown" category with NA for the following categorical variables ####
#########################################################################################################################

### List of variables to update
variables_to_update_Replace_Unknown <- c(
  "SMOKING_STATUS",
  "CCS_ANGINA_STATUS_PRE_PRO_STABLE",
  "NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE",
  "CSHA_CLINICAL_FRAILTY_SCALE_SCORE",
  "CEREBRAL_CIRCULATION_PROTEC_DEVICE",
  "MITRAL_REGURGITATION",
  "DRUGS_AT_DISCHARGE_ANTITHROMBOTIC"
)


##### Replace "Unknown" with NA and retain factor class and drop the "Unknown" level###
clean_cohort <- clean_cohort %>%
  mutate(across(
    all_of(variables_to_update_Replace_Unknown),
    ~ if_else(as.character(.) == "Unknown", NA_character_, as.character(.)) %>%
      factor()  # Recreate the factor to drop unused levels
  ))



###VERIFY THE CHANGES###
### Check summary of variables after the update
summary(clean_cohort$SMOKING_STATUS)
summary(clean_cohort$CCS_ANGINA_STATUS_PRE_PRO_STABLE)
summary(clean_cohort$NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE)
summary(clean_cohort$CSHA_CLINICAL_FRAILTY_SCALE_SCORE)
summary(clean_cohort$CEREBRAL_CIRCULATION_PROTEC_DEVICE)
summary(clean_cohort$MITRAL_REGURGITATION)
summary(clean_cohort$DRUGS_AT_DISCHARGE_ANTITHROMBOTIC)


#########################################################################################################################
########### STEP 13: ADDITIONAL CLEANING, make the following categorical variables binary ###############################
#########################################################################################################################

### Convert ON_DIALYSIS and SEVERE_LIVER_DISEASE to binary variables (Missing values assigned 0)
clean_cohort <- clean_cohort %>%
  mutate(
    BINARY_known_ON_DIALYSIS = if_else(ON_DIALYSIS == "Yes", 1, 0, missing = 0),
    BINARY_known_SEVERE_LIVER_DISEASE = if_else(SEVERE_LIVER_DISEASE == "Yes", 1, 0, missing = 0)
  )


##VERIFY The changes
summary(clean_cohort$BINARY_known_ON_DIALYSIS)
summary(clean_cohort$BINARY_known_SEVERE_LIVER_DISEASE)
sum(clean_cohort$BINARY_known_ON_DIALYSIS)
sum(clean_cohort$BINARY_known_SEVERE_LIVER_DISEASE)



#########################################################################################################################
########### STEP 14: CLEANING, Re-categorize the following categorical variables #########################################
#########################################################################################################################

### Create a binary variable for CVA presence
clean_cohort <- clean_cohort %>%
  mutate(
    BINARY_CVA_UP_TO_DX_known = if_else(
      CVA_UP_TO_DX %in% c("Yes Haemorrhagic", "Yes Ischaemic", "Yes undetermined"),
      1,  # Assign 1 if any type of CVA is present
      0,  # Assign 0 otherwise (including "No", "Unknown", or NA)
      missing = 0  # Handle missing values as 0
    )
  )

summary(clean_cohort$BINARY_CVA_UP_TO_DX_known)
sum(clean_cohort$BINARY_CVA_UP_TO_DX_known)


############################################################
### Create a binary variable for USE_OF_CARDIOPULMONARY_BYPASS
clean_cohort <- clean_cohort %>%
  mutate(
    BINARY_known_USE_OF_CARDIOPULMONARY_BYPASS = if_else(
      USE_OF_CARDIOPULMONARY_BYPASS %in% c("Yes - elective", "Yes - emergency"),
      1,  # Assign 1 if bypass was used (elective or emergency)
      0,  # Assign 0 for "No", "Unknown", or NA
      missing = 0  # Handle missing values as 0
    )
  )


summary(clean_cohort$BINARY_known_USE_OF_CARDIOPULMONARY_BYPASS) 
sum(clean_cohort$BINARY_known_USE_OF_CARDIOPULMONARY_BYPASS)



##################################################################
### Create a binary variable for HISTORY_OF_NEUROLOGICAL_DISEASE
clean_cohort <- clean_cohort %>%
  mutate(
    BINARY_known_HISTORY_OF_NEUROLOGICAL_DISEASE = if_else(
      HISTORY_OF_NEUROLOGICAL_DISEASE %in% c(
        "CVA with full recovery", 
        "CVA with residual deficit", 
        "TIA or RIND", 
        "Other history of neurological dysfunction"
      ),
      1,  # Assign 1 for any history of neurological disease
      0,  # Assign 0 for "No history of neurological disease", "Unknown", or NA
      missing = 0  # Handle missing values as 0
    )
  )

summary(clean_cohort$BINARY_known_HISTORY_OF_NEUROLOGICAL_DISEASE)



##########################################################
# Create a binary variable for PERMANENT_PACING
clean_cohort <- clean_cohort %>%
  mutate(
    BINARY_known_PERMANENT_PACING = if_else(
      PERMANENT_PACING %in% c(
        "Yes-post-procedure", 
        "Yes-pre-procedure prophylactic", 
        "Yes-pre-procedure therapeutic (including distant past)", 
        "Yes-per-procedure"
      ),
      1,  # Assign 1 for any type of permanent pacing
      0,  # Assign 0 for "No", "Unknown", or NA
      missing = 0  # Handle missing values as 0
    )
  )

summary(clean_cohort$BINARY_known_PERMANENT_PACING)
sum(clean_cohort$BINARY_known_PERMANENT_PACING)


##################################################################
##### Create binary variables for the three complications#######
clean_cohort <- clean_cohort %>%
  mutate(
    # VASCULAR_ACCESS_SITE_COMPLICATIONS binary variable
    BINARY_known_VASCULAR_ACCESS_SITE_COMPLICATIONS = if_else(
      VASCULAR_ACCESS_SITE_COMPLICATIONS %in% c("Major", "Minor"),
      1,  # Assign 1 if there's any vascular access site complication
      0,  # Assign 0 for "None", "Unknown", or NA
      missing = 0  # Handle missing values as 0
    ),
    
    # BLEEDING binary variable
    BINARY_known_BLEEDING = if_else(
      BLEEDING %in% c(
        "Yes - Life threatening or disabling", 
        "Yes - Major", 
        "Yes - Minor"
      ),
      1,  # Assign 1 if there is any bleeding event
      0,  # Assign 0 for "No", "Unknown", or NA
      missing = 0  # Handle missing values as 0
    ),
    
    # ACUTE_KIDNEY_INJURY_WITHIN_7_DAYS binary variable
    BINARY_known_ACUTE_KIDNEY_INJURY = if_else(
      ACUTE_KIDNEY_INJURY_WITHIN_7_DAYS %in% c("Stage 1", "Stage 2", "Stage 3"),
      1,  # Assign 1 if there is any AKI
      0,  # Assign 0 for "No AKI", "Unknown", or NA
      missing = 0  # Handle missing values as 0
    )
  )


sum(clean_cohort$BINARY_known_VASCULAR_ACCESS_SITE_COMPLICATIONS)
sum(clean_cohort$BINARY_known_BLEEDING)
sum(clean_cohort$BINARY_known_ACUTE_KIDNEY_INJURY)




#############################################################
###CLEANING: LV_FUNCTION ################
library(dplyr)
install.packages("forcats")
library(forcats)

# Combine "Not measured", "Unknown", and NA into a single category "Not measured or unknown"
clean_cohort <- clean_cohort %>%
  mutate(
    LV_FUNCTION = fct_explicit_na(
      fct_collapse(
        LV_FUNCTION,
        `Not measured or unknown` = c("Not measured", "Unknown")
      ),
      na_level = "Not measured or unknown"
    )
  )

summary(clean_cohort$LV_FUNCTION)
class(clean_cohort$LV_FUNCTION)
levels(clean_cohort$LV_FUNCTION)


#############################################################
#####CLEANING: EXTENT_OF_CORONARY_VESSEL_DISEASE###########
library(dplyr)
library(forcats)

# Combine "Not investigated" and NA into a single category "Not investigated or missing"
clean_cohort <- clean_cohort %>%
  mutate(
    EXTENT_OF_CORONARY_VESSEL_DISEASE = fct_explicit_na(
      fct_collapse(
        EXTENT_OF_CORONARY_VESSEL_DISEASE,
        `Not investigated or missing` = c("Not investigated")
      ),
      na_level = "Not investigated or missing"
    )
  )

summary(clean_cohort$EXTENT_OF_CORONARY_VESSEL_DISEASE)
class(clean_cohort$EXTENT_OF_CORONARY_VESSEL_DISEASE)




##################################################################################################
##### create a binary variable for "Current smoker" (1 for current smokers and 0 for others) #####
##################################################################################################

# Create a binary variable for known current smokers
library(dplyr)
library(tidyr)

# Ensure NAs are set to 0 for binary variable
clean_cohort <- clean_cohort %>%
  mutate(
    binary_known_current_smoker = ifelse(SMOKING_STATUS == "Current smoker", 1, 0),
    binary_known_current_smoker = replace_na(binary_known_current_smoker, 0)  # Replace NA with 0
  )


###Explanation:
#   If SMOKING_STATUS is "Current smoker", the value is 1.
#   If SMOKING_STATUS is NA, it is explicitly set to 0.
#   Otherwise (for any other smoking status), it is also set to 0


# Check the summary to verify
summary(clean_cohort$binary_known_current_smoker)
table(clean_cohort$binary_known_current_smoker)




###################################################################################################################
##### add a new level called "Missing" to the SMOKING_STATUS variable and replace all NAs with this new level #####
###################################################################################################################

# Add "Missing" as a level
levels(clean_cohort$SMOKING_STATUS) <- c(
  levels(clean_cohort$SMOKING_STATUS), "Missing"
)

# Replace NA values with "Missing"
clean_cohort$SMOKING_STATUS[is.na(clean_cohort$SMOKING_STATUS)] <- "Missing"

summary(clean_cohort$SMOKING_STATUS)



#########################################################################################################################
########### STEP 15: CLEANING, exp_cardiac_rehab_6m                             #########################################
#########################################################################################################################

###exp_cardiac_rehab_6m###
# Convert the existing integer variable to numeric
clean_cohort$exp_cardiac_rehab_6m <- as.numeric(clean_cohort$exp_cardiac_rehab_6m)
class(clean_cohort$exp_cardiac_rehab_6m)
sum(clean_cohort$exp_cardiac_rehab_6m)
summary(clean_cohort$exp_cardiac_rehab_6m)

# Create a new binary variable for rehab exposure (1 = exposed to rehab, 0 = not exposed to rehab)
clean_cohort$binary_exposed_rehab <- ifelse(is.na(clean_cohort$exp_cardiac_rehab_6m), 0, 1)
sum(clean_cohort$binary_exposed_rehab)
summary(clean_cohort$binary_exposed_rehab)

# Convert the new binary variable into a factor with labeled levels
clean_cohort$binary_exposed_rehab <- factor(clean_cohort$binary_exposed_rehab, 
                                            levels = c(0, 1), 
                                            labels = c("Not Exposed to Rehab", "Exposed to Rehab"))

# Verify the changes
summary(clean_cohort$binary_exposed_rehab)


summary(clean_cohort$exp_cardiac_rehab_attends_6m)
class(clean_cohort$exp_cardiac_rehab_appt_date_6m)
class(clean_cohort$exp_cardiac_rehab_completion_date_6m)




###CONVERT "exp_cardiac_rehab_appts_6m" and "exp_cardiac_rehab_attends_6m" from integer64 class to integer class ##########

# Convert integer64 to integer
clean_cohort$exp_cardiac_rehab_appts_6m <- as.integer(clean_cohort$exp_cardiac_rehab_appts_6m)
clean_cohort$exp_cardiac_rehab_attends_6m <- as.integer(clean_cohort$exp_cardiac_rehab_attends_6m)

# Verify the conversion
str(clean_cohort$exp_cardiac_rehab_appts_6m)
str(clean_cohort$exp_cardiac_rehab_attends_6m)
summary(clean_cohort$exp_cardiac_rehab_appts_6m)
summary(clean_cohort$exp_cardiac_rehab_attends_6m)





#########################################################################################################################
###########  CLEANING, exp_cardiac_rehab_3m                                     #########################################
#########################################################################################################################



###exp_cardiac_rehab_3m###
# Convert the existing integer variable to numeric
clean_cohort$exp_cardiac_rehab_3m <- as.numeric(clean_cohort$exp_cardiac_rehab_3m)
class(clean_cohort$exp_cardiac_rehab_3m)
sum(clean_cohort$exp_cardiac_rehab_3m)
summary(clean_cohort$exp_cardiac_rehab_3m)

# Create a new binary variable for rehab exposure (1 = exposed to rehab, 0 = not exposed to rehab)
clean_cohort$binary_exposed_rehab_3m <- ifelse(is.na(clean_cohort$exp_cardiac_rehab_3m), 0, 1)
sum(clean_cohort$binary_exposed_rehab_3m)
summary(clean_cohort$binary_exposed_rehab_3m)

# Convert the new binary variable into a factor with labeled levels
clean_cohort$binary_exposed_rehab_3m <- factor(clean_cohort$binary_exposed_rehab_3m, 
                                               levels = c(0, 1), 
                                               labels = c("Not Exposed to Rehab", "Exposed to Rehab within 3 months"))

# Verify the changes
summary(clean_cohort$binary_exposed_rehab_3m)


summary(clean_cohort$exp_cardiac_rehab_attends_3m)
class(clean_cohort$exp_cardiac_rehab_appt_date_3m)
class(clean_cohort$exp_cardiac_rehab_completion_date_3m)




###CONVERT "exp_cardiac_rehab_appts_3m" and "exp_cardiac_rehab_attends_3m" from integer64 class to integer class ##########

# Convert integer64 to integer
clean_cohort$exp_cardiac_rehab_appts_3m <- as.integer(clean_cohort$exp_cardiac_rehab_appts_3m)
clean_cohort$exp_cardiac_rehab_attends_3m <- as.integer(clean_cohort$exp_cardiac_rehab_attends_3m)

# Verify the conversion
str(clean_cohort$exp_cardiac_rehab_appts_3m)
str(clean_cohort$exp_cardiac_rehab_attends_3m)
summary(clean_cohort$exp_cardiac_rehab_appts_3m)
summary(clean_cohort$exp_cardiac_rehab_attends_3m)




#########################################################################################################################
########### STEP 16: GROUPING, exp_cardiac_rehab_attends_6m                    #########################################
#########################################################################################################################

# Categorize exp_cardiac_rehab_attends_6m into meaningful groups (e.g., 0, 1–5, 6–10, >10 appointments) 

library(dplyr)

# Categorize the number of rehab appointments into groupings according to the number of appointements attended
clean_cohort <- clean_cohort %>%
  mutate(
    rehab_attendance_groups = case_when(
      exp_cardiac_rehab_attends_6m == 0 ~ "0 appointments",
      exp_cardiac_rehab_attends_6m > 0 & exp_cardiac_rehab_attends_6m <= 3 ~ "1-3 appointments",
      exp_cardiac_rehab_attends_6m > 3 & exp_cardiac_rehab_attends_6m <= 6 ~ "3-6 appointments",
      exp_cardiac_rehab_attends_6m > 6 & exp_cardiac_rehab_attends_6m <= 10 ~ "6-10 appointments",
      exp_cardiac_rehab_attends_6m > 10 ~ ">10 appointments"
    ),
    # Convert to a factor with levels in logical order
    rehab_attendance_groups = factor(
      rehab_attendance_groups,
      levels = c("0 appointments", "1-3 appointments", "3-6 appointments", "6-10 appointments", ">10 appointments")
    )
  )


# View the unique groups
unique(clean_cohort$rehab_attendance_groups)

# If the variable is a factor
summary(clean_cohort$rehab_attendance_groups)
# Count the number of patients in each group
table(clean_cohort$rehab_attendance_groups)





#########################################################################################################################
########### GROUPING, exp_cardiac_rehab_attends_3m                              #########################################
#########################################################################################################################

# Categorize exp_cardiac_rehab_attends_3m into meaningful groups (e.g., 0, 1–5, 6–10, >10 appointments) 

library(dplyr)

# Categorize the number of rehab appointments within 3 months into groupings according to the number of appointements attended
clean_cohort <- clean_cohort %>%
  mutate(
    rehab_attendance_groups_3m = case_when(
      exp_cardiac_rehab_attends_3m == 0 ~ "0 appointments",
      exp_cardiac_rehab_attends_3m > 0 & exp_cardiac_rehab_attends_3m <= 3 ~ "1-3 appointments",
      exp_cardiac_rehab_attends_3m > 3 & exp_cardiac_rehab_attends_3m <= 6 ~ "3-6 appointments",
      exp_cardiac_rehab_attends_3m > 6 & exp_cardiac_rehab_attends_3m <= 10 ~ "6-10 appointments",
      exp_cardiac_rehab_attends_3m > 10 ~ ">10 appointments"
    ),
    # Convert to a factor with levels in logical order
    rehab_attendance_groups_3m = factor(
      rehab_attendance_groups_3m,
      levels = c("0 appointments", "1-3 appointments", "3-6 appointments", "6-10 appointments", ">10 appointments")
    )
  )


# View the unique groups
unique(clean_cohort$rehab_attendance_groups_3m)

# If the variable is a factor
summary(clean_cohort$rehab_attendance_groups_3m)
# Count the number of patients in each group
table(clean_cohort$rehab_attendance_groups_3m)





#########################################################################################################################
########### STEP 17: Cleaning KATZ_INDEX_OF_INDEPENDENCE_DAILY       ################################################
#########################################################################################################################

###################################
# KATZ_INDEX_OF_INDEPENDENCE_DAILY#
##################################

# Check variable distribution
summary(clean_cohort$KATZ_INDEX_OF_INDEPENDENCE_DAILY)

# structures correctly with min of 0 and max 6 

# Convert to numeric variable
clean_cohort <- clean_cohort %>%
  mutate(KATZ_INDEX_CLEAN = case_when(
    is.na(KATZ_INDEX_OF_INDEPENDENCE_DAILY) ~ NA_real_,  # Keep missing values
    TRUE ~ as.numeric(KATZ_INDEX_OF_INDEPENDENCE_DAILY)  # Ensure numeric
  ))

# Check for missing values
sum(is.na(clean_cohort$KATZ_INDEX_CLEAN))

class(clean_cohort$KATZ_INDEX_CLEAN)
summary(clean_cohort$KATZ_INDEX_CLEAN)


#Recategorize the KATZ_INDEX_OF_INDEPENDENCE_DAILY variable based on classification:

# Recode KATZ_INDEX_OF_INDEPENDENCE_DAILY into new categories
# Ensure all values are categorized properly
clean_cohort <- clean_cohort %>%
  mutate(KATZ_INDEX_CAT = case_when(
    KATZ_INDEX_OF_INDEPENDENCE_DAILY >= 5 & KATZ_INDEX_OF_INDEPENDENCE_DAILY <= 6 ~ "Full Function",
    KATZ_INDEX_OF_INDEPENDENCE_DAILY >= 4 & KATZ_INDEX_OF_INDEPENDENCE_DAILY < 5 ~ "Mild Impairment",
    KATZ_INDEX_OF_INDEPENDENCE_DAILY >= 0 & KATZ_INDEX_OF_INDEPENDENCE_DAILY < 4 ~ "Moderate Impairment, Severe Functional Impairment up to Depednent", 
        TRUE ~ NA_character_  # Ensures missing values remain NA
  )) %>%
  mutate(KATZ_INDEX_CAT = factor(KATZ_INDEX_CAT, levels = c("Moderate Impairment, Severe Functional Impairment up to Depednent", "Mild Impairment",
                                                            "Full Function")))

# Verify that all values have been categorized properly
table(clean_cohort$KATZ_INDEX_CAT, useNA = "ifany")




# Add "Missing" as a level
levels(clean_cohort$KATZ_INDEX_CAT) <- c(
  levels(clean_cohort$KATZ_INDEX_CAT), "Missing"
)

# Replace NA values with "Missing"
clean_cohort$KATZ_INDEX_CAT[is.na(clean_cohort$KATZ_INDEX_CAT)] <- "Missing"

summary(clean_cohort$KATZ_INDEX_CAT)


##########################################################################################################################################################
########### STEP 18: CREATING comorbidity variable, which is a composite of (pulmonary disease, liver disease, neuro, extracardiac vasculopathy)  #########
##########################################################################################################################################################

##comorbidity (pulmonary, liver, neuro, extracardiac vasculopathy
# binary_known_pulmonary_disease
# BINARY_known_SEVERE_LIVER_DISEASE
# BINARY_known_HISTORY_OF_NEUROLOGICAL_DISEASE
# EXTRACARDIAC_ARTERIOPATHY



#first clean EXTRACARDIAC_ARTERIOPATHY
summary(clean_cohort$EXTRACARDIAC_ARTERIOPATHY)
class(clean_cohort$EXTRACARDIAC_ARTERIOPATHY)

# Create a binary variable in clean_cohort

clean_cohort$binary_known_extracardiac_arteriopathy <- ifelse(
  clean_cohort$EXTRACARDIAC_ARTERIOPATHY == "Yes", 1, 0
)

# Convert NA or unknown values to 0
clean_cohort$binary_known_extracardiac_arteriopathy[is.na(clean_cohort$EXTRACARDIAC_ARTERIOPATHY) | 
                                                      clean_cohort$EXTRACARDIAC_ARTERIOPATHY == "Unknown"] <- 0

# Confirm the transformation
summary(clean_cohort$binary_known_extracardiac_arteriopathy)
sum(clean_cohort$binary_known_extracardiac_arteriopathy)
class(clean_cohort$binary_known_extracardiac_arteriopathy)


#####################################################################################################
#CREATE  composite of (pulmonary disease, liver disease, neuro, extracardiac vasculopathy)  #########
#####################################################################################################


##comorbidity (pulmonary, liver, neuro, extracardiac vasculopathy
# binary_known_pulmonary_disease
# BINARY_known_SEVERE_LIVER_DISEASE
# BINARY_known_HISTORY_OF_NEUROLOGICAL_DISEASE
#binary_known_extracardiac_arteriopathy


# Create the composite binary variable
clean_cohort$COMORBIDITY_pulmonary_liver_neuro_excardiacvascu <- ifelse(
  clean_cohort$binary_known_pulmonary_disease == 1 |
    clean_cohort$BINARY_known_SEVERE_LIVER_DISEASE == 1 |
    clean_cohort$BINARY_known_HISTORY_OF_NEUROLOGICAL_DISEASE == 1 |
    clean_cohort$binary_known_extracardiac_arteriopathy == 1, 
  1, 0
)

# Verify the transformation
summary(clean_cohort$COMORBIDITY_pulmonary_liver_neuro_excardiacvascu)
sum(clean_cohort$COMORBIDITY_pulmonary_liver_neuro_excardiacvascu)





##########################################################################################################################################################
########### STEP 19: Cleaning: NYHA pre-procedure #########################################################################################################
##########################################################################################################################################################

# Add "Missing" as a level
levels(clean_cohort$NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE) <- c(
  levels(clean_cohort$NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE), "Missing"
)

# Replace NA with "Missing"
clean_cohort$NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE[is.na(clean_cohort$NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE)] <- "Missing"

# Reorder the levels
clean_cohort$NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE <- factor(
  clean_cohort$NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE,
  levels = c(
    "No limitation of physical activity",
    "Slight limitation of ordinary physical activity",
    "Marked limitation of ordinary physical activity",
    "Symptoms at rest or minimal activity",
    "Missing"
  )
)


summary(clean_cohort$NYHA_DYSPNOEA_STATUS_PRE_PRO_STABLE)





##########################################################################################################################################################
########### STEP 20: Cleaning: PREV_MI_AND_INT_BTW_PRO_AND_LAST_MI ########################################################################################
##########################################################################################################################################################

#Rename the variable, add a "Missing" level, make NAs correspond to "Missing," and reorder the levels 

# Create a new variable called previous_MI_time
clean_cohort$previous_MI_time <- clean_cohort$PREV_MI_AND_INT_BTW_PRO_AND_LAST_MI

# Add "Missing" as a level
levels(clean_cohort$previous_MI_time) <- c(
  levels(clean_cohort$previous_MI_time), "Missing"
)

# Replace NA with "Missing"
clean_cohort$previous_MI_time[is.na(clean_cohort$previous_MI_time)] <- "Missing"

# Reorder the levels
clean_cohort$previous_MI_time <- factor(
  clean_cohort$previous_MI_time,
  levels = c(
    "No previous MI",
    "MI < 6 hours",
    "MI 6-24 hours",
    "MI 1-30 days",
    "MI 31-90 days",
    "MI > 90 days",
    "Missing"  # Keep this as the last level for clarity
  )
)

summary(clean_cohort$previous_MI_time)


###Not enough volumes in some of teh categories therefore regroup into >90days and <90days
# Create the new categorized variable
clean_cohort$previous_MI_time_grouped <- ifelse(
  clean_cohort$previous_MI_time == "MI > 90 days", "MI > 90 days",
  ifelse(clean_cohort$previous_MI_time %in% c("MI < 6 hours", "MI 6-24 hours", "MI 1-30 days", "MI 31-90 days"), 
         "MI < 90 days",
         ifelse(clean_cohort$previous_MI_time == "No previous MI", "No previous MI", "Missing")))

# Convert to a factor with the desired order
clean_cohort$previous_MI_time_grouped <- factor(
  clean_cohort$previous_MI_time_grouped,
  levels = c("MI > 90 days", "MI < 90 days", "No previous MI", "Missing")
)

# Verify the new variable
summary(clean_cohort$previous_MI_time_grouped)





##########################################################################################################################################################
########### STEP 21: Reorder: LV Function  ################################################################################################################
##########################################################################################################################################################

# Reorder the levels
clean_cohort$LV_FUNCTION <- factor(
  clean_cohort$LV_FUNCTION,
  levels = c(
    "Good (LVEF >=50%)",
    "Fair (LVEF = 30-49%)",
    "Poor (LVEF <30%)",
    "Not measured or unknown"
  )
)

summary(clean_cohort$LV_FUNCTION)



##########################################################################################################################################################
########### STEP 22: Clean: Mitral Regurgitation  #########################################################################################################
##########################################################################################################################################################


# Add "Missing" as a level
levels(clean_cohort$MITRAL_REGURGITATION) <- c(
  levels(clean_cohort$MITRAL_REGURGITATION), "Missing"
)

# Replace NA with "Missing"
clean_cohort$MITRAL_REGURGITATION[is.na(clean_cohort$MITRAL_REGURGITATION)] <- "Missing"

# Reorder the levels
clean_cohort$MITRAL_REGURGITATION <- factor(
  clean_cohort$MITRAL_REGURGITATION,
  levels = c(
    "None",
    "Mild",
    "Moderate",
    "Severe",
    "Missing"
  )
)


summary(clean_cohort$MITRAL_REGURGITATION)




##########################################################################################################################################################
########### STEP 23: Clean: PA_SYSTOLIC_PRESSURE_MMHG and procedure time ##################################################################################
##########################################################################################################################################################

# Change the class to numeric
clean_cohort$PA_SYSTOLIC_PRESSURE_MMHG <- as.numeric(clean_cohort$PA_SYSTOLIC_PRESSURE_MMHG)

summary(clean_cohort$PA_SYSTOLIC_PRESSURE_MMHG)


clean_cohort$PROCEDURE_TIME_MINS <- as.numeric(clean_cohort$PROCEDURE_TIME_MINS)
summary(clean_cohort$PROCEDURE_TIME_MINS)

class(clean_cohort$PROCEDURE_TIME_MINS)


##########################################################################################################################################################
########### STEP 24: Clean: DEVICE_FAILURE_REFERS_TO_VALVE_ONLY ###########################################################################################
##########################################################################################################################################################


# Replace "Unknown" with "Missing or unknown"
clean_cohort$DEVICE_FAILURE_REFERS_TO_VALVE_ONLY <- as.character(clean_cohort$DEVICE_FAILURE_REFERS_TO_VALVE_ONLY)
clean_cohort$DEVICE_FAILURE_REFERS_TO_VALVE_ONLY[
  clean_cohort$DEVICE_FAILURE_REFERS_TO_VALVE_ONLY == "Unknown"
] <- "Missing or unknown"

# Replace NA values with "Missing or unknown"
clean_cohort$DEVICE_FAILURE_REFERS_TO_VALVE_ONLY[is.na(clean_cohort$DEVICE_FAILURE_REFERS_TO_VALVE_ONLY)] <- "Missing or unknown"

# Convert back to a factor and reorder levels
clean_cohort$DEVICE_FAILURE_REFERS_TO_VALVE_ONLY <- factor(
  clean_cohort$DEVICE_FAILURE_REFERS_TO_VALVE_ONLY,
  levels = c(
    "No failure",
    "Probably iatrogenic",
    "Probably intrinsic",
    "Missing or unknown"
  )
)


summary(clean_cohort$DEVICE_FAILURE_REFERS_TO_VALVE_ONLY)



### Create the binary variable
clean_cohort$binary_known_device_failure_valve <- ifelse(
  is.na(clean_cohort$DEVICE_FAILURE_REFERS_TO_VALVE_ONLY) | 
    clean_cohort$DEVICE_FAILURE_REFERS_TO_VALVE_ONLY %in% c("No failure", "Missing or unknown"), 
  0, 
  1
)

sum(clean_cohort$binary_known_device_failure_valve)


##########################################################################################################################################################
########### STEP 25: Create new varaible valve_NOT_succesfully_deployed ##################################################################################
##########################################################################################################################################################

# Create the new variable
clean_cohort$valve_NOT_succesfully_deployed <- ifelse(
  clean_cohort$VALVE_SUCCESSFULLY_DEPLOYED == "No", "Yes",
  ifelse(clean_cohort$VALVE_SUCCESSFULLY_DEPLOYED == "Yes", "No",
         ifelse(clean_cohort$VALVE_SUCCESSFULLY_DEPLOYED == "Unknown" | is.na(clean_cohort$VALVE_SUCCESSFULLY_DEPLOYED), "Unknown or missing", NA)
  )
)

clean_cohort$valve_NOT_succesfully_deployed[is.na(clean_cohort$VALVE_SUCCESSFULLY_DEPLOYED)] <- "Unknown or missing"

# Convert the new variable to a factor
clean_cohort$valve_NOT_succesfully_deployed <- factor(
  clean_cohort$valve_NOT_succesfully_deployed,
  levels = c("Yes", "No", "Unknown or missing")
)


summary(clean_cohort$valve_NOT_succesfully_deployed)


### Create the binary variable
clean_cohort$binary_known_valve_NOT_succesfully_deployed <- ifelse(
  clean_cohort$valve_NOT_succesfully_deployed == "Yes", 
  1, 
  0
)


sum(clean_cohort$binary_known_valve_NOT_succesfully_deployed)



##########################################################################################################################################################
########### STEP 26: Create new binary variable FURTHER_VALVE_INTERVEN_B4_DISCHARG########################################################################
##########################################################################################################################################################
summary(clean_cohort$FURTHER_VALVE_INTERVEN_B4_DISCHARG)

clean_cohort$binary_further_valve_intervention_b4_discharge <- ifelse(
  clean_cohort$FURTHER_VALVE_INTERVEN_B4_DISCHARG %in% c("None", "Unknown") | 
    is.na(clean_cohort$FURTHER_VALVE_INTERVEN_B4_DISCHARG), 
  0, 
  1
)

summary(clean_cohort$binary_further_valve_intervention_b4_discharge)
sum(clean_cohort$binary_further_valve_intervention_b4_discharge)
class(clean_cohort$binary_further_valve_intervention_b4_discharge)






##########################################################################################################################################################
########### STEP 27: Create new binary variable conversion to full sternotomy ########################################################################
##########################################################################################################################################################

# Create the binary variable
clean_cohort$binary_known_sternotomy <- ifelse(
  is.na(clean_cohort$CONV_TO_FULL_STERNOTOMY_DURING_PRO) | 
    clean_cohort$CONV_TO_FULL_STERNOTOMY_DURING_PRO %in% c("No", "Unknown"), 
  0, 
  1
)

summary(clean_cohort$binary_known_sternotomy)
sum(clean_cohort$binary_known_sternotomy)


summary(clean_cohort$CONV_TO_FULL_STERNOTOMY_DURING_PRO)





##########################################################################################################################################################
########### STEP 27: Create new binary variable TAMPONADE_DURING_POST_PROCEDURE###########################################################################
##########################################################################################################################################################

# Create the binary variable
clean_cohort$binary_known_tamponade <- ifelse(
  is.na(clean_cohort$TAMPONADE_DURING_POST_PROCEDURE) | 
    clean_cohort$TAMPONADE_DURING_POST_PROCEDURE == "No", 
  0, 
  1
)

sum(clean_cohort$binary_known_tamponade)




##########################################################################################################################################################
########### STEP 28: Create new binary variable peri-procedural MI ######################################################################################
##########################################################################################################################################################


# Create the binary variable
clean_cohort$binary_known_peri_proc_MI_72hrs_after_TAVI <- ifelse(
  is.na(clean_cohort$PERI_PROCEDURAL_MI_72HRS_AFTER_PRO) | 
    clean_cohort$PERI_PROCEDURAL_MI_72HRS_AFTER_PRO %in% c("No", "Unknown"), 
  0, 
  1
)

sum(clean_cohort$binary_known_peri_proc_MI_72hrs_after_TAVI)



##########################################################################################################################################################
########### STEP 29: Create new binary variable CVA up to discharge ######################################################################################
##########################################################################################################################################################

# Create the binary variable
clean_cohort$binary_known_CVA_b4_dx <- ifelse(
  is.na(clean_cohort$CVA_UP_TO_DX) | 
    clean_cohort$CVA_UP_TO_DX %in% c("No", "Unknown"), 
  0, 
  1
)

sum(clean_cohort$binary_known_CVA_b4_dx)




##########################################################################################################################################################
########### STEP 30: Create new binary variable vascular access complications ############################################################################
##########################################################################################################################################################

# Create the binary variable
clean_cohort$binary_known_vasc_access_site_complics <- ifelse(
  is.na(clean_cohort$VASCULAR_ACCESS_SITE_COMPLICATIONS) | 
    clean_cohort$VASCULAR_ACCESS_SITE_COMPLICATIONS %in% c("None", "Unknown"), 
  0, 
  1
)

sum(clean_cohort$binary_known_vasc_access_site_complics)

summary(clean_cohort$VASCULAR_ACCESS_SITE_COMPLICATIONS)


sum(clean_cohort$BINARY_known_ACUTE_KIDNEY_INJURY)


##########################################################################################################################################################
########### STEP 30: Create new composite variable: TAVI_Procedural_complications     ####################################################################
##########################################################################################################################################################

#make composite variable "TAVI_Procedural_complications " b/c the volumes of the complication variables are too small individually

# Variables to include:
  # DEVICE_FAILURE_REFERS_TO_VALVE_ONLY (n=85) 
  # valve_NOT_successfully_deployed (n=623)
  # FURTHER_VALVE_INTERVEN_B4_DISCHARG (n=188)
  # conversion to full sternotomy (n=79)
  # tamponade (n=153)
  # peri-procedural MI (n=60)
  # CVA (n=391)
  # vascular access complics (n=1,268)
  # acute kidney injury (n=596)



# Create the composite variable, which is having any of these binary variables 
clean_cohort$TAVI_procedural_complications <- ifelse(
  rowSums(clean_cohort[, c(
    "binary_known_device_failure_valve",
    "binary_known_valve_NOT_succesfully_deployed",
    "binary_further_valve_intervention_b4_discharge",
    "binary_known_sternotomy",
    "binary_known_tamponade",
    "binary_known_peri_proc_MI_72hrs_after_TAVI",
    "binary_known_CVA_b4_dx",
    "binary_known_vasc_access_site_complics",
    "BINARY_known_ACUTE_KIDNEY_INJURY"
  )] == 1, na.rm = TRUE) > 0, 
  1, 
  0
)


sum(clean_cohort$TAVI_procedural_complications)
summary(clean_cohort$TAVI_procedural_complications)




##########################################################################################################################################################
########### STEP 31: Create new variable:age groups    ####################################################################
##########################################################################################################################################################



# Create Age Groups in 10-year Blocks
clean_cohort <- clean_cohort %>%
  mutate(age_groups = case_when(
    age_at_tavi < 70 ~ "<70",
    age_at_tavi >= 70 & age_at_tavi < 80 ~ "70-79",
    age_at_tavi >= 80 & age_at_tavi < 90 ~ "80-89",
    age_at_tavi >= 90 ~ "90+"
  ))

# Convert to Factor (Ordered for Analysis)
clean_cohort <- clean_cohort %>%
  mutate(age_groups = factor(age_groups, levels = c("<70", "70-79", "80-89", "90+")))

# Check Age Group Distribution
table(clean_cohort$age_groups)







##########################################################################################################################################################
########### STEP 32: Create new variable: days since TAVI    ####################################################################
##########################################################################################################################################################


#a.	create a new variable: tavi_cips_disdate_plus6m =  tavi_cips_disdate + 6 months
#b. create a new variable: tavi_cips_disdate_plus365days =  tavi_cips_disdate + 365 days

library(dplyr)
library(lubridate)

clean_cohort <- clean_cohort %>%
  # Step 1: Add 180 days to the TAVI discharge date
  mutate(
    tavi_cips_disdate_plus180days = tavi_cips_disdate + days(180),
    tavi_cips_disdate_plus365_days = tavi_cips_disdate + days(365),
    tavi_cips_disdate_plus90_days = tavi_cips_disdate + days(90)
  )




##########################################################################################################################################################
########### STEP 33: Create new variable: cardiac_rehab_follow_up_days    ####################################################################
##########################################################################################################################################################

#ACTIONS:
#a.	create a new variable: tavi_cips_disdate_plus6m =  tavi_cips_disdate + 6 months
#b.	create another new variable: cardiac_rehab_end_follow_up = min(tavi_cips_disdate_plus6m, date_of_death)
#c.	create final variable:  cardiac_rehab_follow_up_days = date_diff(cardiac_rehab_end_follow_up, tavi_cips_disdate)

library(dplyr)
library(lubridate)

clean_cohort <- clean_cohort %>%
  # Step 1: Add 180 days to the TAVI discharge date
  mutate(
    tavi_cips_disdate_plus180days = tavi_cips_disdate + days(180),
    
    # Step 2: Calculate the end of follow-up for cardiac rehab
    cardiac_rehab_end_follow_up = pmin(tavi_cips_disdate_plus180days, date_of_death, na.rm = TRUE),
    
    # Step 3: Calculate the follow-up days for cardiac rehab
    cardiac_rehab_follow_up_days = as.numeric(difftime(cardiac_rehab_end_follow_up, tavi_cips_disdate, units = "days"))
  )



summary(clean_cohort$cardiac_rehab_follow_up_days)


##########################################################################################################################################################
########### STEP 34: Create new variable: year of tavi discharge (categorical)    ####################################################################
##########################################################################################################################################################



# Extract the year from TAVI discharge date
clean_cohort <- clean_cohort %>%
  mutate(year_tavi_discharge = year(tavi_cips_disdate))  # Extract year


#convert to a factor to compare different years independently

clean_cohort <- clean_cohort %>%
  mutate(year_tavi_discharge = as.factor(year_tavi_discharge))

summary(clean_cohort$year_tavi_discharge)




#######################################################################################################################################################


