#####################################################################################################################################################
#####  Incidence of outcome  Explore the rate of all cause mortality by the total follow-up time    ########################
#####################################################################################################################################################


### all cause mortality rate per quarter,
####aggregate the data by Discharge_TAVI_year_quarter and compute the sum of date_of_death == 1 and mortality censor time days for each quarter,then divide the two to get the rate (all cause rehosp rate over follwoup time in days)

library(dplyr)


# Define study start and end dates
study_end <- as.Date("2024-03-31")

library(lubridate)  # Ensure lubridate is loaded

clean_cohort <- clean_cohort %>%
  # Step 1: start of followup window
  mutate(
    tavi_cips_disdate_plus180days = tavi_cips_disdate + days(180),
    quarter_180days_post_tavi = paste0(year(tavi_cips_disdate_plus180days), "-Q", quarter(tavi_cips_disdate_plus180days)),  # Get year & quarter
    
    # Step 2: Calculate the end of follow-up 
    all_cause_mortality_end_follow_up = pmin(study_end, date_of_death, na.rm = TRUE),
    
    # Step 3: Calculate the follow-up days
    all_cause_mortality_follow_up_days = as.numeric(difftime(all_cause_mortality_end_follow_up, tavi_cips_disdate_plus180days, units = "days")),
    
    # Define mortality indicator (1 for death, 0 for no death)
    mortality_indicator_all_cause = ifelse(!is.na(date_of_death), 1, 0)
  )



summary(clean_cohort$all_cause_mortality_follow_up_days)
sum(clean_cohort$mortality_indicator_all_cause)


# Check where all_cause_mortality_follow_up_days is negative (this is normal b/c start of follwoup window is 180 days post TAVI discharge, so need to filter out deaths within 180 days)

count_invalid_all_cause_mortality_follow_up_days <- sum(clean_cohort$all_cause_mortality_follow_up_days <= 0, na.rm = TRUE)
print(count_invalid_all_cause_mortality_follow_up_days)

##note 1295 (rounded to nearest 5) patients have a negative or 0 all cause mortality followup time, this equates to deaths within 180 days of tavi discharge

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
count_invalid_all_cause_mortality_follow_up_days <- clean_cohort %>%
  group_by(exp_cardiac_rehab_6m) %>%  # Group by rehab exposure
  summarise(
    count_invalid = sum(
      all_cause_mortality_follow_up_days <= 0, na.rm = TRUE
    ) 
  )

# Print results
print(count_invalid_all_cause_mortality_follow_up_days)
### 40 (rounded to nearest 5) deaths in rehab group and 1250 (rounded to nearest 5) deaths in no CR group, within 180 days of tavi discharge




##Set Negative Values to NA 
##If negative values indicate invalid data that cannot be fixed, exclude them by setting such cases to NA:

clean_cohort <- clean_cohort %>%
  mutate(
    all_cause_mortality_follow_up_days = ifelse(all_cause_mortality_follow_up_days < 0, NA, all_cause_mortality_follow_up_days)
  )

summary(clean_cohort$all_cause_mortality_follow_up_days)




###calculate the mortality rate per quarter,
####aggregate the data by quarter_180days_post_tavi and compute the sum of date_of_death >= tavi_cips_disdate_plus180days == 1 and all_cause_mortality_follow_up_days for each quarter,then divide the two to get the rate (mortality rate over follwoup time in days)

library(dplyr)

# Step 1: Filter the data for valid follow-up days
clean_cohort <- clean_cohort %>%
  filter(all_cause_mortality_follow_up_days > 0)

summary(clean_cohort$all_cause_mortality_follow_up_days)


#select deaths greater than 180days post tavi

clean_cohort <- clean_cohort |> 
  mutate(mortality_all_cause_post180 = ifelse((date_of_death >= tavi_cips_disdate_plus180days), 1, 0))

# Replace NA values with 0
clean_cohort$mortality_all_cause_post180[is.na(clean_cohort$mortality_all_cause_post180)] <- 0

sum(clean_cohort$mortality_all_cause_post180)
#NOTES: 7805 (rounded to nearest 5) deaths occurred post 180 days of tavi


all_cause_mortality_post180_summary <- clean_cohort |>  
  group_by(quarter_180days_post_tavi, binary_exposed_rehab) |>  # Group by rehab exposure
  summarise(
    total_patients = round(n_distinct(person_id) / 5) * 5,  # Unique patients in each group, rounded
    total_all_cause_mortality = round(sum(mortality_all_cause_post180) / 5) * 5,  # Total all_cause deaths, rounded
    avg_all_cause_mortality_per_patient = total_all_cause_mortality / total_patients  # Avg mortality rate per patient
  )

print(all_cause_mortality_post180_summary)




# Calculate all_cause mortality rate per quarter
all_cause_mortality_rate_per_quarter <- clean_cohort %>%
  group_by(quarter_180days_post_tavi) %>%   ##THIS is the quarter of TAVI_DISCHARGE plus 180days on the x-axis
  summarise(
    total_all_cause_mortality = 5*round(sum(mortality_all_cause_post180 == 1, na.rm = TRUE)/5),
    total_all_cause_mortality_follow_up_days = sum(all_cause_mortality_follow_up_days, na.rm = TRUE)
  ) %>%
  mutate(
    all_cause_mortality_rate = total_all_cause_mortality / total_all_cause_mortality_follow_up_days
  )

# Step 3: View the resulting table
print(all_cause_mortality_rate_per_quarter)




library(ggplot2)

# Convert mortality rate to per 10,000 follow-up days
all_cause_mortality_rate_per_quarter <- all_cause_mortality_rate_per_quarter %>%
  mutate(all_cause_mortality_rate_per_10000_days = all_cause_mortality_rate * 10000)

print(all_cause_mortality_rate_per_quarter)

# Save the rounded table as a CSV file for extraction
output_file_all_cause_mortality_rate_per_quarter_total_cohort <- 
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_mortality/Mortality rate per quarter trends/all_cause_mortality_rate_per_quarter_total_cohort.csv"
write.csv(all_cause_mortality_rate_per_quarter, file = output_file_all_cause_mortality_rate_per_quarter_total_cohort, row.names = FALSE)




ggplot(all_cause_mortality_rate_per_quarter, aes(x = quarter_180days_post_tavi, y = all_cause_mortality_rate_per_10000_days, group = 1)) +
  # Add blue shading for COVID period (March 2020 - December 2021)
  annotate("rect",
           xmin = which(all_cause_mortality_rate_per_quarter$quarter_180days_post_tavi == "2020-Q1"),
           xmax = which(all_cause_mortality_rate_per_quarter$quarter_180days_post_tavi == "2021-Q4"),
           ymin = 0, ymax = Inf, fill = "blue", alpha = 0.1) +
  # Add line for mortality rate per quarter
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_text(aes(label = sprintf("%.2f", all_cause_mortality_rate_per_10000_days)), vjust = -1, size = 3.5, color = "black") +
  scale_y_continuous(limits = c(0, 7)) + 
  labs(
    x = "Quarter (180 days post TAVI)",
    y = "All Cause Mortality Rate per 10,000 Person Days",
    title = "All Cause Mortality Rate per Quarter",
    subtitle = "Number of all cause deaths per 10,000 person days of follow-up by quarter. Shaded in blue: COVID period"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)
  )








###exposed vs not exposed to rehab


all_cause_mortality_rate_per_quarter_per_group <- clean_cohort %>%
  group_by(quarter_180days_post_tavi, binary_exposed_rehab) %>%  # Group by rehab exposure
  summarise(
    total_all_cause_mortality = 5*round(sum(mortality_all_cause_post180 == 1, na.rm = TRUE)/5),
    total_all_cause_mortality_follow_up_days = sum(all_cause_mortality_follow_up_days, na.rm = TRUE)
  ) %>%
  mutate(
    all_cause_mortality_rate = total_all_cause_mortality / total_all_cause_mortality_follow_up_days,
    all_cause_mortality_rate_per_10000_days = all_cause_mortality_rate * 10000,
    # NEW LINE: create display version of counts
    total_all_cause_mortality_display = ifelse(total_all_cause_mortality < 10, "<10", as.character(total_all_cause_mortality))
  )


print(all_cause_mortality_rate_per_quarter_per_group)




all_cause_mortality_rate_per_quarter_per_group %>%
  select(-total_all_cause_mortality) 

all_cause_mortality_rate_per_quarter_per_group <- all_cause_mortality_rate_per_quarter_per_group %>%
  select(
    quarter_180days_post_tavi,
    binary_exposed_rehab,
    total_all_cause_mortality_display,
    total_all_cause_mortality_follow_up_days,
    all_cause_mortality_rate,
    all_cause_mortality_rate_per_10000_days
  )

print(all_cause_mortality_rate_per_quarter_per_group)



# Save the rounded table as a CSV file for extraction
output_file_all_cause_mortality_rate_per_quarter_per_group <- 
  "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/all_cause_mortality/Mortality rate per quarter trends/all_cause_mortality_rate_per_quarter_per_group.csv"
write.csv(all_cause_mortality_rate_per_quarter_per_group, file = output_file_all_cause_mortality_rate_per_quarter_per_group, row.names = FALSE)





library(ggplot2)


# Convert mortality rate per 10,000 follow-up days
all_cause_mortality_rate_per_quarter_per_group <- all_cause_mortality_rate_per_quarter_per_group %>%
  mutate(all_cause_mortality_rate_per_10000_days_per_group = all_cause_mortality_rate * 10000)



# Plot with two lines for rehab exposure vs non-exposure
ggplot(all_cause_mortality_rate_per_quarter_per_group, aes(x = quarter_180days_post_tavi, 
                                                        y = all_cause_mortality_rate_per_10000_days_per_group, 
                                                        color = binary_exposed_rehab, 
                                                        group = binary_exposed_rehab)) + 
  
  
  # Add lines and points for each group (exposed vs unexposed)
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  # Add labels
  geom_text(aes(label = sprintf("%.2f", all_cause_mortality_rate_per_10000_days_per_group)), 
            vjust = -1, size = 3.5, color = "black") +
  
  # Customize labels
  
  scale_y_continuous(limits = c(0, 7)) +
  labs(
    x = "Quarter (180 days post TAVI)",
    y = "All Cause Mortality Rate per 10,000 Person Days",
    title = "All Cause Mortality Rate per Quarter",
    subtitle = "Number of all cause deaths per 10,000 person days of follow-up by quarter",
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
ggplot(all_cause_mortality_rate_per_quarter_per_group, aes(x = quarter_180days_post_tavi, 
                                                        y = all_cause_mortality_rate_per_10000_days_per_group, 
                                                        color = binary_exposed_rehab, 
                                                        group = binary_exposed_rehab)) + 
  # Add UK lockdown shading (blue)
  geom_tile(data = lockdown_shading, 
            aes(x = quarter, y = ymin, fill = "UK Lockdown"),
            width = 1, height = Inf, alpha = 0.2, inherit.aes = FALSE) +
  
  geom_line(size = 1) +
  geom_point(size = 2.2) +
  
  geom_text(aes(label = sprintf("%.2f", all_cause_mortality_rate_per_10000_days_per_group)),
            position = position_jitter(width = 0.4, height = 0.1),
            vjust = -1.2, size = 5, color = "black") +
  
  # Ensure X-axis shows quarter labels properly
  scale_x_discrete(name = "Quarter") +
  
  scale_y_continuous(limits = c(0, 7)) +
  
  scale_fill_manual(values = c("UK Lockdown" = "blue")) +  # Add shading legend
  labs(
    x = "Quarter",
    y = "All Cause Mortality Rate per 10,000 Person Days",
    title = "All Cause Mortality Rate per Quarter",
    subtitle = "Shaded in blue: UK Lockdowns (Q2 2020, Q4 2020, Q1 2021)",
    fill = "Shaded Period",  # Legend title for shading
    color = "Cardiac Rehabilitation Exposure"
  ) +
  
  theme_minimal(base_size = 21) +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),  # Ensure readable quarter labels
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )


# Step 5: Save PDF
ggsave("all_cause_mortality_plot_figure4_manuscript.pdf", plot = all_cause_mortality_plot_figure4_manuscript, width = 10, height = 6, units = "in")







################################################################################################################################################




##Improve format with covid shading
# Add quarter_index for numeric x-axis
all_cause_mortality_rate_per_quarter_per_group$quarter_index <- as.numeric(factor(all_cause_mortality_rate_per_quarter_per_group$quarter_180days_post_tavi))

# Define precise lockdown shading windows based on index
lockdown_shading <- data.frame(
  xmin = c(7.9, 10.6, 11.0),   # Start of each lockdown
  xmax = c(8.6, 10.9, 11.8),  # End of each lockdown
  ymin = 0,
  ymax = Inf,
  lockdown = c("1st Lockdown", "2nd Lockdown", "3rd Lockdown")
)

# Plot
all_cause_mortality_rate_per_quarter_CRvs_NoCR_figure4_manuscript <- ggplot(all_cause_mortality_rate_per_quarter_per_group, aes(x = quarter_index, 
                                                           y = all_cause_mortality_rate_per_10000_days_per_group, 
                                                           color = binary_exposed_rehab, 
                                                           group = binary_exposed_rehab)) +
  
  # Lockdown shading
  geom_rect(data = lockdown_shading,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = lockdown),
            inherit.aes = FALSE, alpha = 0.2) +
  
  # Line and points
  geom_line(size = 1.5) +
  geom_point(size = 3.2) +
  
  # Labels on points
  geom_text(aes(label = sprintf("%.2f", all_cause_mortality_rate_per_10000_days_per_group)),
            position = position_jitter(width = 0.1, height = 0.1),
            vjust = -1.2, size = 3, color = "black") +
  
  # X-axis
  scale_x_continuous(
    breaks = all_cause_mortality_rate_per_quarter_per_group$quarter_index,
    labels = all_cause_mortality_rate_per_quarter_per_group$quarter_180days_post_tavi,
    name = "Quarter"
  ) +
  
  # Y-axis
  scale_y_continuous(limits = c(0, 6)) +
  
  # Fill legend
  scale_fill_manual(values = c("1st Lockdown" = "blue", 
                               "2nd Lockdown" = "blue", 
                               "3rd Lockdown" = "blue")) +
  
  # Labels
  labs(
    y = "All-cause Mortality Rate per 10,000 Person Days",
    title = "All-cause Mortality Rate per Quarter",
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
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.position = "right"
  )




ggsave("all_cause_mortality_rate_per_quarter_CRvs_NoCR_figure4_manuscript.pdf", width = 11, height = 6, units = "in")















