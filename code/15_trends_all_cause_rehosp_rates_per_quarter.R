#####################################################################################################################################################
#####  Incidence of outcome  Explore the rate of all cause rehosp by the total rehosp follow-up time    ########################
#####################################################################################################################################################


### all cause rehosp rate per quarter,
####aggregate the data by Discharge_TAVI_year_quarter and compute the sum of out_readm_all_cause_flag == 1 and rehosp censor time days for each quarter,then divide the two to get the rate (all cause rehosp rate over follwoup time in days)

library(dplyr)


# Define study start and end dates
study_end <- as.Date("2024-03-31")

library(lubridate)  # Ensure lubridate is loaded

clean_cohort <- clean_cohort %>%
  # Step 1: start of followup window
  mutate(
    tavi_cips_disdate_plus180days = tavi_cips_disdate + days(180),
    quarter_180days_post_tavi = paste0(year(tavi_cips_disdate_plus180days), "-Q", quarter(tavi_cips_disdate_plus180days)),  # Get year & quarter
    
    # Step 2: Calculate the end of follow-up for all cause rehosp
    all_cause_rehosp_end_follow_up = pmin(study_end, date_of_death, na.rm = TRUE),
    
    # Step 3: Calculate the follow-up days for all cause rehosp
    all_cause_rehosp_follow_up_days = as.numeric(difftime(all_cause_rehosp_end_follow_up, tavi_cips_disdate_plus180days, units = "days"))
  )





summary(clean_cohort$all_cause_rehosp_follow_up_days)



# Check where all_cause_rehosp_follow_up_days is negative (this is normal b/c start of follwoup window is 180 days post TAVI discharge, so need to filter out deaths within 180 days)

count_invalid_all_cause_rehosp_follow_up_days <- sum(clean_cohort$all_cause_rehosp_follow_up_days <= 0, na.rm = TRUE)
print(count_invalid_all_cause_rehosp_follow_up_days)

##note 1295 (rounded to nearest 5) patients have a negative or 0 all cause rehosp followup time, this equates to deaths within 180 days of tavi discharge

##check number deaths within 180 days of tavi discharge
library(dplyr)

# Calculate deaths within 180 days of TAVI discharge
deaths_within_180_days <- clean_cohort %>%
  filter(!is.na(date_of_death) & !is.na(tavi_cips_disdate)) %>%  # Ensure dates exist
  mutate(days_since_discharge = as.numeric(date_of_death - tavi_cips_disdate)) %>%  # Compute difference
  filter(days_since_discharge <= 180 & days_since_discharge > 0)  # Keep deaths within 180 days

# Count the number of deaths within 180 days
num_deaths_180_days <- nrow(deaths_within_180_days)

# Print result
print(paste("Number of deaths within 180 days of TAVI discharge:", num_deaths_180_days))

### 1295 (rounded to nearest 5) deaths within 180 days of TAVI discharge


# Grouping by rehab_indicator and summarizing count
count_invalid_all_cause_rehosp_follow_up_days <- clean_cohort %>%
  group_by(exp_cardiac_rehab_6m) %>%  # Group by rehab exposure
  summarise(
    count_invalid = sum(
      all_cause_rehosp_follow_up_days <= 0, na.rm = TRUE
    ) 
  )

# Print results
print(count_invalid_all_cause_rehosp_follow_up_days)
### 40 (rounded to nearest 5) deaths in rehab group and 1250 (rounded to nearest 5) deaths in no CR group, within 180 days of tavi discharge




##Set Negative Values to NA 
##If negative values indicate invalid data that cannot be fixed, exclude them by setting such cases to NA:

clean_cohort <- clean_cohort %>%
  mutate(
    all_cause_rehosp_follow_up_days = ifelse(all_cause_rehosp_follow_up_days < 0, NA, all_cause_rehosp_follow_up_days)
  )

summary(clean_cohort$all_cause_rehosp_follow_up_days)




###calculate the rehosp rate per quarter,
####aggregate the data by quarter_180days_post_tavi and compute the sum of out_readm_all_cause_date >= tavi_cips_disdate_plus180days == 1 and all_cause_rehosp_follow_up_days for each quarter,then divide the two to get the rate (rehosp rate over follwoup time in days)

library(dplyr)

# Step 1: Filter the data for valid follow-up days
clean_cohort <- clean_cohort %>%
  filter(all_cause_rehosp_follow_up_days > 0)

summary(clean_cohort$all_cause_rehosp_follow_up_days)


#select rehosps greater than 180days post tavi

clean_cohort <- clean_cohort |> 
  mutate(rehosp_all_cause_post180 = ifelse((out_readm_all_cause_date >= tavi_cips_disdate_plus180days), 1, 0))

# Replace NA values with 0
clean_cohort$rehosp_all_cause_post180[is.na(clean_cohort$rehosp_all_cause_post180)] <- 0



all_cause_rehosp_post180_summary <- clean_cohort |>  
  group_by(quarter_180days_post_tavi, binary_exposed_rehab) |>  # Group by rehab exposure
  summarise(
    total_patients = round(n_distinct(person_id) / 5) * 5,  # Unique patients in each group, rounded
    total_all_cause_rehosp = round(sum(rehosp_all_cause_post180) / 5) * 5,  # Total all_cause rehospitalizations, rounded
    avg_all_cause_rehosp_per_patient = total_all_cause_rehosp / total_patients  # Avg HF rehosp per patient
  )

print(all_cause_rehosp_post180_summary)




# Calculate all_cause rehosp rate per quarter
all_cause_rehosp_rate_per_quarter <- clean_cohort %>%
  group_by(quarter_180days_post_tavi) %>%   ##THIS is the quarter of TAVI_DISCHARGE plus 180days on the x-axis
  summarise(
    total_all_cause_rehosps = 5*round(sum(rehosp_all_cause_post180 == 1, na.rm = TRUE)/5),
    total_all_cause_rehosp_follow_up_days = sum(all_cause_rehosp_follow_up_days, na.rm = TRUE)
  ) %>%
  mutate(
    all_cause_rehosp_rate = total_all_cause_rehosps / total_all_cause_rehosp_follow_up_days
  )

# Step 3: View the resulting table
print(all_cause_rehosp_rate_per_quarter)




library(ggplot2)

# Convert rehosp rate to per 10,000 follow-up days
all_cause_rehosp_rate_per_quarter <- all_cause_rehosp_rate_per_quarter %>%
  mutate(all_cause_rehosp_rate_per_10000_days = all_cause_rehosp_rate * 10000)

print(all_cause_rehosp_rate_per_quarter)

# Save the rounded table as a CSV file for extraction
output_file_all_cause_rehosp_rate_per_quarter_total_cohort <- 
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_rehosp/All_cause_rehosp_rate_per_quarter_trends/all_cause_rehosp_rate_per_quarter_total_cohort.csv"
write.csv(all_cause_rehosp_rate_per_quarter, file = output_file_all_cause_rehosp_rate_per_quarter_total_cohort, row.names = FALSE)




ggplot(all_cause_rehosp_rate_per_quarter, aes(x = quarter_180days_post_tavi, y = all_cause_rehosp_rate_per_10000_days, group = 1)) +
  # Add blue shading for COVID period (March 2020 - December 2021)
  annotate("rect",
           xmin = which(all_cause_rehosp_rate_per_quarter$quarter_180days_post_tavi == "2020-Q1"),
           xmax = which(all_cause_rehosp_rate_per_quarter$quarter_180days_post_tavi == "2021-Q4"),
           ymin = 0, ymax = Inf, fill = "blue", alpha = 0.1) +
  # Add line for rehosp rate per quarter
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_text(aes(label = sprintf("%.2f", all_cause_rehosp_rate_per_10000_days)), vjust = -1, size = 3.5, color = "black") +
  scale_y_continuous(limits = c(0, 9)) + 
  labs(
    x = "Quarter (180 days post TAVI)",
    y = "All Cause rehospitilisation Rate per 10,000 Person Days",
    title = "All Cause Rehospitilisation Rate per Quarter",
    subtitle = "Number of all cause rehospitilisations per 10,000 Person Days of follow-up by quarter"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)
  )








###exposed vs not exposed to rehab


all_cause_rehosp_rate_per_quarter_per_group <- clean_cohort %>%
  group_by(quarter_180days_post_tavi, binary_exposed_rehab) %>%  # Group by rehab exposure
  summarise(
    total_all_cause_rehosps = 5*round(sum(rehosp_all_cause_post180 == 1, na.rm = TRUE)/5),
    total_all_cause_rehosp_follow_up_days = sum(all_cause_rehosp_follow_up_days, na.rm = TRUE)
  ) %>%
  mutate(
    all_cause_rehosp_rate = total_all_cause_rehosps / total_all_cause_rehosp_follow_up_days,
    all_cause_rehosp_rate_per_10000_days = all_cause_rehosp_rate * 10000
  )


print(all_cause_rehosp_rate_per_quarter_per_group)



# Save the rounded table as a CSV file for extraction
output_file_all_cause_rehosp_rate_per_quarter_per_group <- 
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_rehosp/All_cause_rehosp_rate_per_quarter_trends/all_cause_rehosp_rate_per_quarter_per_group.csv"
write.csv(all_cause_rehosp_rate_per_quarter_per_group, file = output_file_all_cause_rehosp_rate_per_quarter_per_group, row.names = FALSE)





library(ggplot2)


# Convert rehosp rate per 10,000 follow-up days
all_cause_rehosp_rate_per_quarter_per_group <- all_cause_rehosp_rate_per_quarter_per_group %>%
  mutate(all_cause_rehosp_rate_per_10000_days_per_group = all_cause_rehosp_rate * 10000)



# Plot with two lines for rehab exposure vs non-exposure
ggplot(all_cause_rehosp_rate_per_quarter_per_group, aes(x = quarter_180days_post_tavi, 
                                                 y = all_cause_rehosp_rate_per_10000_days_per_group, 
                                                 color = binary_exposed_rehab, 
                                                 group = binary_exposed_rehab)) + 
  
  
  # Add lines and points for each group (exposed vs unexposed)
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  # Add labels
  geom_text(aes(label = sprintf("%.2f", all_cause_rehosp_rate_per_10000_days_per_group)), 
            vjust = -1, size = 3.5, color = "black") +
  
  # Customize labels
  
  scale_y_continuous(limits = c(0, 9)) +
  labs(
    x = "Quarter (180 days post TAVI)",
    y = "All Cause Rehospitalization Rate per 10,000 Person Days",
    title = "All Cause Rehospitalization Rate per Quarter",
    subtitle = "Number of all cause rehospitalizations per 10,000 Person Days of follow-up by quarter",
    color = "Cardiac Rehabilitation Exposure"
  ) +
  
  # Improve theme
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)
  )












library(ggplot2)
library(dplyr)

# Define lockdown shading periods using "YYYY-QX" format
lockdown_shading <- data.frame(
  quarter = c("2020-Q2", "2020-Q4", "2021-Q1"),  # UK lockdowns
  ymin = 0, ymax = Inf  # Full Y-axis shading
)



# Plot with UK Lockdown Shading using geom_tile()
ggplot(all_cause_rehosp_rate_per_quarter_per_group, aes(x = quarter_180days_post_tavi, 
                                                 y = all_cause_rehosp_rate_per_10000_days_per_group, 
                                                 color = binary_exposed_rehab, 
                                                 group = binary_exposed_rehab)) + 
  # Add UK lockdown shading (blue)
  geom_tile(data = lockdown_shading, 
            aes(x = quarter, y = ymin, fill = "UK Lockdown"),
            width = 1, height = Inf, alpha = 0.2, inherit.aes = FALSE) +
  
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  geom_text(aes(label = sprintf("%.2f", all_cause_rehosp_rate_per_10000_days_per_group)),
            position = position_jitter(width = 0.3, height = 0),
            vjust = -1.2, size = 2.5, color = "black") +
  
  # Ensure X-axis shows quarter labels properly
  scale_x_discrete(name = "Quarter") +
  
  scale_y_continuous(limits = c(0, 9)) +
  
  scale_fill_manual(values = c("UK Lockdown" = "blue")) +  # Add shading legend
  labs(
    x = "Quarter",
    y = "All Cause Rehospitalization Rate per 10,000 Person Days",
    title = "All Cause Rehospitalization Rate per Quarter",
    subtitle = "Shaded in blue: UK Lockdowns (Q2 2020, Q4 2020, Q1 2021)",
    fill = "Shaded Period",  # Legend title for shading
    color = "Cardiac Rehabilitation Exposure"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Ensure readable quarter labels
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)
  )


#######################################################################################################################################################################################################


##Improve format with covid shading
# Add quarter_index for numeric x-axis
all_cause_rehosp_rate_per_quarter_per_group$quarter_index <- as.numeric(factor(all_cause_rehosp_rate_per_quarter_per_group$quarter_180days_post_tavi))

# Define precise lockdown shading windows based on index
lockdown_shading <- data.frame(
  xmin = c(7.9, 10.6, 11.0),   # Start of each lockdown
  xmax = c(8.6, 10.9, 11.8),  # End of each lockdown
  ymin = 0,
  ymax = Inf,
  lockdown = c("1st Lockdown", "2nd Lockdown", "3rd Lockdown")
)

# Plot 
all_cause_rehospitalisation_rate_per_quarter_CRvs_NoCR_figure5_manuscript <- ggplot(all_cause_rehosp_rate_per_quarter_per_group, aes(x = quarter_index, 
                                                                                                                      y = all_cause_rehosp_rate_per_10000_days_per_group, 
                                                                                                                      color = binary_exposed_rehab, 
                                                                                                                      group = binary_exposed_rehab)) +
  
  # Lockdown shading
  geom_rect(data = lockdown_shading,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = lockdown),
            inherit.aes = FALSE, alpha = 0.2) +
  
  # Line and points
  geom_line(size = 1.5) +
  geom_point(size = 2.5) +
  
  # Labels on points
  geom_text(aes(label = sprintf("%.2f", all_cause_rehosp_rate_per_10000_days_per_group)),
            position = position_jitter(width = 0.25, height = 0.1),
            vjust = -1.2, size = 2, color = "black") +
  
  # X-axis
  scale_x_continuous(
    breaks = all_cause_rehosp_rate_per_quarter_per_group$quarter_index,
    labels = all_cause_rehosp_rate_per_quarter_per_group$quarter_180days_post_tavi,
    name = "Quarter"
  ) +
  
  # Y-axis
  scale_y_continuous(limits = c(0,9.5)) +
  
  # Fill legend
  scale_fill_manual(values = c("1st Lockdown" = "blue", 
                               "2nd Lockdown" = "blue", 
                               "3rd Lockdown" = "blue")) +
  
  # Labels
  labs(
    y = "All-Cause Rehospitalisation Rate per 10,000 Person Days",
    title = "All-Cause Rehospitalisation Rate per Quarter",
    subtitle = "Shaded in blue: UK Lockdowns (Mar-May 2020, Nov-Dec 2020, Jan-Mar 2021)",
    fill = "Lockdown Period",
    color = "Rehabilitation Exposure"
  ) +
  
  # Theme
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 8),
    plot.subtitle = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.position = "right"
  )





ggsave("all_cause_rehospitalisation_rate_per_quarter_CRvs_NoCR_figure5_manuscript.pdf", width = 10, height = 6, units = "in")




