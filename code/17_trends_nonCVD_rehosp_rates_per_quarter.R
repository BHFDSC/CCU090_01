#####################################################################################################################################################
#####  Incidence of outcome  Explore the rate of Non CVD rehosp by the total rehosp follow-up time    ########################
#####################################################################################################################################################


### non cvd rehosp rate per quarter,
####aggregate the data by Discharge_TAVI_year_quarter and compute the sum of out_readm_non_cvd_flag == 1 and rehosp censor time days for each quarter,then divide the two to get the rate (non cvd rehosp rate over follwoup time in days)

library(dplyr)


# Define study start and end dates
study_end <- as.Date("2024-03-31")

library(lubridate)  # Ensure lubridate is loaded

clean_cohort <- clean_cohort %>%
  # Step 1: start of followup window
  mutate(
    tavi_cips_disdate_plus180days = tavi_cips_disdate + days(180),
    quarter_180days_post_tavi = paste0(year(tavi_cips_disdate_plus180days), "-Q", quarter(tavi_cips_disdate_plus180days)),  # Get year & quarter
    
    # Step 2: Calculate the end of follow-up for non cvd rehosp
    non_cvd_rehosp_end_follow_up = pmin(study_end, date_of_death, na.rm = TRUE),
    
    # Step 3: Calculate the follow-up days for non cvd rehosp
    non_cvd_rehosp_follow_up_days = as.numeric(difftime(non_cvd_rehosp_end_follow_up, tavi_cips_disdate_plus180days, units = "days"))
  )





summary(clean_cohort$non_cvd_rehosp_follow_up_days)



# Check where non_cvd_rehosp_follow_up_days is negative (this is normal b/c start of follwoup window is 180 days post TAVI discharge, so need to filter out deaths within 180 days)

count_invalid_non_cvd_rehosp_follow_up_days <- sum(clean_cohort$non_cvd_rehosp_follow_up_days <= 0, na.rm = TRUE)
print(count_invalid_non_cvd_rehosp_follow_up_days)

##note 1295 (rounded to nearest 5) patients have a negative or 0 rehosp followup time, this equates to deaths within 180 days of tavi discharge

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
count_invalid_non_cvd_rehosp_follow_up_days <- clean_cohort %>%
  group_by(exp_cardiac_rehab_6m) %>%  # Group by rehab exposure
  summarise(
    count_invalid = sum(
      non_cvd_rehosp_follow_up_days <= 0, na.rm = TRUE
    ) 
  )

# Print results
print(count_invalid_non_cvd_rehosp_follow_up_days)
### 40 (rounded to nearest 5) deaths in rehab group and 1250 (rounded to nearest 5) deaths in no CR group, within 180 days of tavi discharge




##Set Negative Values to NA 
##If negative values indicate invalid data that cannot be fixed, exclude them by setting such cases to NA:

clean_cohort <- clean_cohort %>%
  mutate(
    non_cvd_rehosp_follow_up_days = ifelse(non_cvd_rehosp_follow_up_days < 0, NA, non_cvd_rehosp_follow_up_days)
  )

summary(clean_cohort$non_cvd_rehosp_follow_up_days)




###calculate the rehosp rate per quarter,
####aggregate the data by quarter_180days_post_tavi and compute the sum of out_readm_non_cvd_date >= tavi_cips_disdate_plus180days == 1 and non_cvd_rehosp_follow_up_days for each quarter,then divide the two to get the rate (rehosp rate over follwoup time in days)

library(dplyr)

# Step 1: Filter the data for valid follow-up days
clean_cohort <- clean_cohort %>%
  filter(non_cvd_rehosp_follow_up_days > 0)

summary(clean_cohort$non_cvd_rehosp_follow_up_days)


#select rehosps greater than 180days post tavi

clean_cohort <- clean_cohort |> 
  mutate(rehosp_non_cvd_post180 = ifelse((out_readm_non_cvd_date >= tavi_cips_disdate_plus180days), 1, 0))

# Replace NA values with 0
clean_cohort$rehosp_non_cvd_post180[is.na(clean_cohort$rehosp_non_cvd_post180)] <- 0



non_cvd_rehosp_post180_summary <- clean_cohort |>  
  group_by(quarter_180days_post_tavi, binary_exposed_rehab) |>  # Group by rehab exposure
  summarise(
    total_patients = round(n_distinct(person_id) / 5) * 5,  # Unique patients in each group, rounded
    total_non_cvd_rehosp = round(sum(rehosp_non_cvd_post180) / 5) * 5,  # Total non_cvd rehospitalizations, rounded
    avg_non_cvd_rehosp_per_patient = total_non_cvd_rehosp / total_patients  # Avg non cvd rehosp per patient
  )

print(non_cvd_rehosp_post180_summary)




# Calculate non_cvd rehosp rate per quarter
non_cvd_rehosp_rate_per_quarter <- clean_cohort %>%
  group_by(quarter_180days_post_tavi) %>%   ##THIS is the quarter of TAVI_DISCHARGE plus 180days on the x-axis
  summarise(
    total_non_cvd_rehosps = 5*round(sum(rehosp_non_cvd_post180 == 1, na.rm = TRUE)/5),
    total_non_cvd_rehosp_follow_up_days = sum(non_cvd_rehosp_follow_up_days, na.rm = TRUE)
  ) %>%
  mutate(
    non_cvd_rehosp_rate = total_non_cvd_rehosps / total_non_cvd_rehosp_follow_up_days
  )

# Step 3: View the resulting table
print(non_cvd_rehosp_rate_per_quarter)




library(ggplot2)

# Convert rehosp rate to per 10,000 follow-up days
non_cvd_rehosp_rate_per_quarter <- non_cvd_rehosp_rate_per_quarter %>%
  mutate(non_cvd_rehosp_rate_per_10000_days = non_cvd_rehosp_rate * 10000)

print(non_cvd_rehosp_rate_per_quarter)

# Save the rounded table as a CSV file for extraction
output_file_non_cvd_rehosp_rate_per_quarter_total_cohort <- 
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/non cvd rehosp/non_cvd_rehosps_per_quarter_trends/non_cvd_rehosp_rate_per_quarter_total_cohort.csv"
write.csv(non_cvd_rehosp_rate_per_quarter, file = output_file_non_cvd_rehosp_rate_per_quarter_total_cohort, row.names = FALSE)




ggplot(non_cvd_rehosp_rate_per_quarter, aes(x = quarter_180days_post_tavi, y = non_cvd_rehosp_rate_per_10000_days, group = 1)) +
  # Add blue shading for COVID period (March 2020 - December 2021)
  annotate("rect",
           xmin = which(non_cvd_rehosp_rate_per_quarter$quarter_180days_post_tavi == "2020-Q1"),
           xmax = which(non_cvd_rehosp_rate_per_quarter$quarter_180days_post_tavi == "2021-Q4"),
           ymin = 0, ymax = Inf, fill = "blue", alpha = 0.1) +
  # Add line for rehosp rate per quarter
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_text(aes(label = sprintf("%.2f", non_cvd_rehosp_rate_per_10000_days)), vjust = -1, size = 3.5, color = "black") +
  scale_y_continuous(limits = c(0, 9)) + 
  labs(
    x = "Quarter (180 days post TAVI)",
    y = "Non-CVD rehospitilisation Rate per 10,000 Person Days",
    title = "Non-CVD Rehospitilisation Rate per Quarter",
    subtitle = "Number of non-CVD rehospitilisations per 10,000 person days of follow-up by quarter. COVID shading in blue"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)
  )








###exposed vs not exposed to rehab


non_cvd_rehosp_rate_per_quarter_per_group <- clean_cohort %>%
  group_by(quarter_180days_post_tavi, binary_exposed_rehab) %>%  # Group by rehab exposure
  summarise(
    total_non_cvd_rehosps = 5*round(sum(rehosp_non_cvd_post180 == 1, na.rm = TRUE)/5),
    total_non_cvd_rehosp_follow_up_days = sum(non_cvd_rehosp_follow_up_days, na.rm = TRUE)
  ) %>%
  mutate(
    non_cvd_rehosp_rate = total_non_cvd_rehosps / total_non_cvd_rehosp_follow_up_days,
    non_cvd_rehosp_rate_per_10000_days = non_cvd_rehosp_rate * 10000
  )


print(non_cvd_rehosp_rate_per_quarter_per_group)



# Save the rounded table as a CSV file for extraction
output_file_non_cvd_rehosp_rate_per_quarter_per_group <- 
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/non cvd rehosp/non_cvd_rehosps_per_quarter_trends/non_cvd_rehosp_rate_per_quarter_per_group.csv"
write.csv(non_cvd_rehosp_rate_per_quarter_per_group, file = output_file_non_cvd_rehosp_rate_per_quarter_per_group, row.names = FALSE)





library(ggplot2)


# Convert rehosp rate per 10,000 follow-up days
non_cvd_rehosp_rate_per_quarter_per_group <- non_cvd_rehosp_rate_per_quarter_per_group %>%
  mutate(non_cvd_rehosp_rate_per_10000_days_per_group = non_cvd_rehosp_rate * 10000)



# Plot with two lines for rehab exposure vs non-exposure
ggplot(non_cvd_rehosp_rate_per_quarter_per_group, aes(x = quarter_180days_post_tavi, 
                                                        y = non_cvd_rehosp_rate_per_10000_days_per_group, 
                                                        color = binary_exposed_rehab, 
                                                        group = binary_exposed_rehab)) + 
  
  
  # Add lines and points for each group (exposed vs unexposed)
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  # Add labels
  geom_text(aes(label = sprintf("%.2f", non_cvd_rehosp_rate_per_10000_days_per_group)), 
            vjust = -1, size = 3.5, color = "black") +
  
  # Customize labels
  
  scale_y_continuous(limits = c(0, 9)) +
  labs(
    x = "Quarter (180 days post TAVI)",
    y = "Non-CVD Rehospitalization Rate per 10,000 Person Days",
    title = "Non-CVD Rehospitalization Rate per Quarter",
    subtitle = "Number of non-CVD rehospitalizations per 10,000 Person Days of follow-up by quarter",
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
ggplot(non_cvd_rehosp_rate_per_quarter_per_group, aes(x = quarter_180days_post_tavi, 
                                                        y = non_cvd_rehosp_rate_per_10000_days_per_group, 
                                                        color = binary_exposed_rehab, 
                                                        group = binary_exposed_rehab)) + 
  # Add UK lockdown shading (blue)
  geom_tile(data = lockdown_shading, 
            aes(x = quarter, y = ymin, fill = "UK Lockdown"),
            width = 1, height = Inf, alpha = 0.2, inherit.aes = FALSE) +
  
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  geom_text(aes(label = sprintf("%.2f", non_cvd_rehosp_rate_per_10000_days_per_group)),
            position = position_jitter(width = 0.3, height = 0),
            vjust = -1.2, size = 2.5, color = "black") +
  
  # Ensure X-axis shows quarter labels properly
  scale_x_discrete(name = "Quarter") +
  
  scale_y_continuous(limits = c(0, 9)) +
  
  scale_fill_manual(values = c("UK Lockdown" = "blue")) +  # Add shading legend
  labs(
    x = "Quarter",
    y = "Non-CVD Rehospitalization Rate per 10,000 Person Days",
    title = "Non-CVD Rehospitalization Rate per Quarter",
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
non_cvd_rehosp_rate_per_quarter_per_group$quarter_index <- as.numeric(factor(non_cvd_rehosp_rate_per_quarter_per_group$quarter_180days_post_tavi))

# Define precise lockdown shading windows based on index
lockdown_shading <- data.frame(
  xmin = c(7.9, 10.6, 11.0),   # Start of each lockdown
  xmax = c(8.6, 10.9, 11.8),  # End of each lockdown
  ymin = 0,
  ymax = Inf,
  lockdown = c("1st Lockdown", "2nd Lockdown", "3rd Lockdown")
)

# Plot 
non_cvd_rehospitalisation_rate_per_quarter_CRvs_NoCR_figure5_manuscript <- ggplot(non_cvd_rehosp_rate_per_quarter_per_group, aes(x = quarter_index, 
                                                                                                                                     y = non_cvd_rehosp_rate_per_10000_days_per_group, 
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
  geom_text(aes(label = sprintf("%.2f", non_cvd_rehosp_rate_per_10000_days_per_group)),
            position = position_jitter(width = 0.225, height = 0),
            vjust = -1.2, size = 2, color = "black") +
  
  # X-axis
  scale_x_continuous(
    breaks = non_cvd_rehosp_rate_per_quarter_per_group$quarter_index,
    labels = non_cvd_rehosp_rate_per_quarter_per_group$quarter_180days_post_tavi,
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
    y = "Non-CVD Rehospitalisation Rate per 10,000 Person Days",
    title = "Non-Cardiovascular Rehospitalisation Rate per Quarter",
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




ggsave("non_cvd_rehospitalisation_rate_per_quarter_CRvs_NoCR_figure6_manuscript.pdf", width = 10, height = 6, units = "in")












