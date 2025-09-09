#####################################################################################################################################################
#####  OBJECTIVE 1b:  Explore the rate of cardiac rehabilitation exposure by the total cardiac rehabilitation follow-up time ########################
#####################################################################################################################################################


##### Rehab rate defined as  = number who attended cardiac rehabilitation / total number of cardiac rehab follow-up time of TAVI patients############

##NOTE: BECAREFUL NOT TO OVERRIDE THE PREVIOUS CODE'S DATA FRAMES AND VECORS

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



### Check Data for Erroneous Dates
# Check where cardiac_rehab_follow_up_days is negative
clean_cohort %>%
  filter(cardiac_rehab_follow_up_days < 0) %>%
  select(tavi_cips_disdate, tavi_cips_disdate_plus180days, date_of_death, cardiac_rehab_end_follow_up, cardiac_rehab_follow_up_days) %>%
  print(n = Inf)


###NOTES: 0 patients died during the TAVI procedure - based on our data curation





###calculate the rehab rate per quarter,
####aggregate the data by Discharge_TAVI_year_quarter and compute the sum of binary_exposed_rehab == 1 and cardiac_rehab_follow_up_days for each quarter,then divide the two to get the rate (rehab rate over follwoup time in days)

library(dplyr)

# Step 1: Filter the data for valid follow-up days
clean_cohort <- clean_cohort %>%
  filter(cardiac_rehab_follow_up_days > 0)

# Step 2: Calculate rehab rate per quarter
cardiac_rehab_rate_per_quarter <- clean_cohort %>%
  group_by(Discharge_TAVI_year_quarter) %>%
  summarise(
    total_rehab_exposed = sum(exp_cardiac_rehab_6m == 1, na.rm = TRUE),
    total_follow_up_days = sum(cardiac_rehab_follow_up_days, na.rm = TRUE)
  ) %>%
  mutate(
    rehab_rate = total_rehab_exposed / total_follow_up_days
  )

# Step 3: View the resulting table
print(cardiac_rehab_rate_per_quarter)



#Calculate rehab rate per quarter, counts rounded to nearest 5
cardiac_rehab_rate_per_quarter_rounded <- clean_cohort %>%
  group_by(Discharge_TAVI_year_quarter) %>%
  summarise(
    total_rehab_exposed = round(sum(exp_cardiac_rehab_6m == 1, na.rm = TRUE) / 5)*5,
    total_follow_up_days = sum(cardiac_rehab_follow_up_days, na.rm = TRUE)
  ) %>%
  mutate(
    rehab_rate = total_rehab_exposed / total_follow_up_days
  )

# Step 3: View the resulting table
print(cardiac_rehab_rate_per_quarter_rounded)


# Save the rounded table as a CSV file for extraction
output_file_cardiac_rehab_rate_per_quarter_rounded <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/rates and vols/cardiac_rehab_rate_per_quarter_rounded.csv"
write.csv(cardiac_rehab_rate_per_quarter_rounded, file = output_file_cardiac_rehab_rate_per_quarter_rounded, row.names = FALSE)






#####visualize the rehab rate per quarter,  create a line or bar graph:
  
library(ggplot2)

ggplot(cardiac_rehab_rate_per_quarter_rounded, aes(x = Discharge_TAVI_year_quarter, y = rehab_rate)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Quarter",
    y = "Rehab Rate",
    title = "Rehab Rate per Quarter"
  ) +
  theme_minimal()




###IMPROVE VISUALISATION

library(ggplot2)

# Convert rehab rate to Rate of rehab per 1000 Person days
cardiac_rehab_rate_per_quarter_rounded <- cardiac_rehab_rate_per_quarter_rounded %>%
  mutate(rehab_rate_per_1000_days = rehab_rate * 1000)

ggplot(cardiac_rehab_rate_per_quarter_rounded, aes(x = Discharge_TAVI_year_quarter, y = rehab_rate_per_1000_days)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rehab_rate_per_1000_days)), vjust = -0.5, size = 3.5, color = "black") +
  labs(
    x = "Quarter",
    y = "Cardiac Rehab Rate per 1000 Person Days",
    title = "Rehab Rate per Quarter",
    subtitle = "Number of rehab exposures per 1000 Person Days of follow-up by quarter"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)
  )




###adjust the code to create a line graph instead of a bar graph to visualize the rehab rate per 1000 Person Days over quarters
library(ggplot2)

# Convert rehab rate to per 1000 follow-up days
cardiac_rehab_rate_per_quarter_rounded <- cardiac_rehab_rate_per_quarter_rounded %>%
  mutate(rehab_rate_per_1000_days = rehab_rate * 1000)

ggplot(cardiac_rehab_rate_per_quarter_rounded, aes(x = Discharge_TAVI_year_quarter, y = rehab_rate_per_1000_days, group = 1)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_text(aes(label = sprintf("%.2f", rehab_rate_per_1000_days)), vjust = -1, size = 3.5, color = "black") +
  labs(
    x = "Quarter",
    y = "Cardiac Rehab Rate per 1000 Person Days",
    title = "Rehab Rate per Quarter",
    subtitle = "Number of Cardiac rehab exposures per 1000 Person Days of follow-up by quarter"
  ) +
  scale_y_continuous(limits = c(0, 0.4), expand = c(0, 0)) +  # <- this sets y-min to 0
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)
  )





###adjust the code to create a line graph instead of a bar graph to visualize the rehab rate per 10,000 Person Days per quarter
library(ggplot2)

# Convert rehab rate to per 10,000 follow-up days
cardiac_rehab_rate_per_quarter_rounded <- cardiac_rehab_rate_per_quarter_rounded %>%
  mutate(rehab_rate_per_10000_days = rehab_rate * 10000)

ggplot(cardiac_rehab_rate_per_quarter_rounded, aes(x = Discharge_TAVI_year_quarter, y = rehab_rate_per_10000_days, group = 1)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_text(aes(label = sprintf("%.2f", rehab_rate_per_10000_days)), vjust = -1, size = 3.5, color = "black") +
  labs(
    x = "Quarter",
    y = "Cardiac Rehab Rate per 10,000 Person Days",
    title = "Rehab Rate per Quarter",
    subtitle = "Number of Cardiac rehab exposures per 10,000 Person Days of follow-up by quarter"
  ) +
  scale_y_continuous(limits = c(0, 4.0), expand = c(0, 0)) +  # <- this sets y-min to 0
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)
  )



### ADD RED SHADING DURING COVID ###

# Convert rehab rate to per 10,000 follow-up days
cardiac_rehab_rate_per_quarter_rounded <- cardiac_rehab_rate_per_quarter_rounded %>%
  mutate(rehab_rate_per_10000_days = rehab_rate * 10000)

ggplot(cardiac_rehab_rate_per_quarter_rounded, aes(x = Discharge_TAVI_year_quarter, y = rehab_rate_per_10000_days, group = 1)) +
  # Add red shading for COVID period (March 2020 - December 2021)
  annotate("rect",
           xmin = which(cardiac_rehab_rate_per_quarter$Discharge_TAVI_year_quarter == "2020-Q1"),
           xmax = which(cardiac_rehab_rate_per_quarter$Discharge_TAVI_year_quarter == "2021-Q4"),
           ymin = 0, ymax = Inf, fill = "red", alpha = 0.1) +
  # Add line for rehab rate per quarter
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 3) +
  # Add labels for values
  geom_text(aes(label = sprintf("%.2f", rehab_rate_per_10000_days)),
            vjust = -1, size = 3.5, color = 'black') +
  # Labels and title
  labs(
    x = "Quarter",
    y = "Cardiac Rehab Rate per 10,000 Person Days (post TAVI)",
    title = "Cardiac Rehabilitation Rate per Quarter",
    subtitle = "Number of Cardiac Rehab Exposures per 10,000 Person Days of Follow-up (post TAVI) by Quarter, with COVID period highlighted"
  ) +
  scale_y_continuous(limits = c(0, 4.0), expand = c(0, 0)) +  # <- this sets y-min to 0
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)
  )


###################################################################################################################################

#Update covid period to separate lockdowns

ggplot(cardiac_rehab_rate_per_quarter_rounded, aes(x = Discharge_TAVI_year_quarter, y = rehab_rate_per_10000_days, group = 1)) +
  
  # COVID period shading (March 2020 - December 2021)
  annotate("rect",
           xmin = which(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter == "2020-Q1"),
           xmax = which(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter == "2021-Q4"),
           ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") +
  
  # Lockdown period (highlighted separately)
  annotate("rect",
           xmin = which(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter == "2020-Q1"),
           xmax = which(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter == "2020-Q2"),
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "blue") +
  annotate("rect",
           xmin = which(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter == "2021-Q4"),
           xmax = which(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter == "2021-Q4") + 0.9,
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "blue") +
  
  # Rehab rate line and points
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 3) +
  
  # Text labels
  geom_text(aes(label = sprintf("%.2f", rehab_rate_per_10000_days)),
            vjust = -1, size = 3.5, color = "black") +
  
  # Axis labels and title
  labs(
    x = "Quarter",
    y = "Cardiac Rehab Rate per 10,000 Person Days (post TAVI)",
    title = "Cardiac Rehabilitation Rate per Quarter",
    subtitle = "Number of Cardiac Rehab Exposures per 10,000 Person Days of Follow-up (post TAVI) by Quarter,\nwith COVID and lockdown periods highlighted"
  ) +
  
  # Adjust y-axis scale
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  
  # Theme adjustments
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)
  )



###################################################################################################################################################
#Update covid period to separate lockdowns and better format

ggplot(cardiac_rehab_rate_per_quarter_rounded, aes(x = Discharge_TAVI_year_quarter, y = rehab_rate_per_10000_days, group = 1)) +
  
  # First lockdown (Mar 23-May 13, 2020)
  annotate("rect",
           xmin = which(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter == "2020-Q1"),
           xmax = which(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter == "2020-Q2"),
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "blue") +
  
  # Second lockdown (Nov 5-Dec 20, 2020)
  annotate("rect",
           xmin = which(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter == "2020-Q4"),
           xmax = which(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter == "2020-Q4") + 0.9,
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "blue") +
  
  # Third lockdown (Jan 6-Mar 29, 2021)
  annotate("rect",
           xmin = which(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter == "2021-Q1"),
           xmax = which(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter == "2021-Q1") + 0.9,
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "blue") +
  
  # Rehab rate line and points
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 3) +
  
  # Text labels
  geom_text(aes(label = sprintf("%.2f", rehab_rate_per_10000_days)),
            vjust = -1, size = 6.5, color = "black") +
  
  # Axis labels and title
  labs(
    x = "Quarter",
    y = "Cardiac Rehab Rate per 10,000 Person Days (post TAVI)",
    title = "Cardiac Rehabilitation Rate per Quarter",
    subtitle = "Number of Cardiac Rehab Exposures per 10,000 Person Days of Follow-up (post TAVI) by Quarter,\nwith COVID lockdown periods highlighted"
  ) +
  
  # Adjust y-axis scale
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  
  # Theme and font adjustments
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)
  )



###################################################################################################################################################
#Update covid period to separate lockdowns and even better format

#Add numeric quarter index
cardiac_rehab_rate_per_quarter_rounded$quarter_index <- seq_along(cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter)

#Create lockdown data
lockdowns <- data.frame(
  lockdown = c("COVID-19 lockdown period", "COVID-19 lockdown period", "COVID-19 lockdown period"),
  xmin = c(9.67, 12.33, 13.07),
  xmax = c(10.5, 12.9, 13.97)
)

#plot
cardiac_rehab_rate_per_quarter <- ggplot(cardiac_rehab_rate_per_quarter_rounded, aes(x = quarter_index, y = rehab_rate_per_10000_days)) +
  
  # Shaded lockdown periods with legend
  geom_rect(data = lockdowns,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = lockdown),
            inherit.aes = FALSE, alpha = 0.2) +
  
  # Rehab rate line and points
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 3) +
  
  # Text labels
  geom_text(aes(label = sprintf("%.2f", rehab_rate_per_10000_days)),
            vjust = -1, size = 3, color = "black") +
  
  # Set custom x-axis labels to quarter names
  scale_x_continuous(
    breaks = cardiac_rehab_rate_per_quarter_rounded$quarter_index,
    labels = cardiac_rehab_rate_per_quarter_rounded$Discharge_TAVI_year_quarter
  ) +
  
  # Y-axis
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  
  # Custom fill color for legend
  scale_fill_manual(values = c("COVID-19 lockdown period" = "blue")) +
  
  # Titles
  labs(
    x = "Quarter",
    y = "Cardiac Rehabilitation Rate per 10,000 Person Days (post TAVI)",
    title = "Cardiac Rehabilitation Rate per Quarter",
    subtitle = "Cardiac Rehab Exposures per 10,000 Person Days of Follow-up (post TAVI) by Quarter,\nwith COVID lockdown periods highlighted",
    fill = NULL
  ) +
  
  # Theme and text formatting
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.9, size = 11),
    plot.title = element_text(face = "bold", size = 8),
    plot.subtitle = element_text(size = 8),
    legend.position = "right"
  )


# Save as PDF
ggsave("cardiac_rehab_rate_per_quarter.pdf", plot = cardiac_rehab_rate_per_quarter, width = 10, height = 6, units = "in")



#####################################################################################################################################################
#####  Annual rates per 10,000 person days of followup                                                                       ########################
#####################################################################################################################################################

# Step 1: Convert quarter to year format
cardiac_rehab_rate_per_year <- clean_cohort %>%
  mutate(Discharge_TAVI_year = substr(as.character(Discharge_TAVI_year_quarter), 1, 4)) %>%  # Extract year
  group_by(Discharge_TAVI_year) %>%
  summarise(
    total_rehab_exposed = sum(exp_cardiac_rehab_6m == 1, na.rm = TRUE),
    total_follow_up_days = sum(cardiac_rehab_follow_up_days, na.rm = TRUE)
  ) %>%
  mutate(
    rehab_rate_per_10000_days = (total_rehab_exposed / total_follow_up_days) * 10000
  )

# Print results
print(cardiac_rehab_rate_per_year)




#Calculate rehab rate per year, counts rounded to nearest 5
cardiac_rehab_rate_per_year_rounded <- clean_cohort %>%
  mutate(Discharge_TAVI_year = substr(as.character(Discharge_TAVI_year_quarter), 1, 4)) %>%  # Extract year
  group_by(Discharge_TAVI_year) %>%
  summarise(
    total_rehab_exposed = round(sum(exp_cardiac_rehab_6m == 1, na.rm = TRUE) / 5)*5,
    total_follow_up_days = sum(cardiac_rehab_follow_up_days, na.rm = TRUE)
  ) %>%
  mutate(
    rehab_rate_per_10000_days = (total_rehab_exposed / total_follow_up_days) * 10000
  )

# Print results
print(cardiac_rehab_rate_per_year_rounded)


# Save the rounded table as a CSV file for extraction
output_file_cardiac_rehab_rate_per_year_rounded <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/rates and vols/cardiac_rehab_rate_per_year_rounded.csv"
write.csv(cardiac_rehab_rate_per_year_rounded, file = output_file_cardiac_rehab_rate_per_year_rounded, row.names = FALSE)





# Step 2: Visualize rehab rate per 10,000 follow-up days per year
library(ggplot2)

ggplot(cardiac_rehab_rate_per_year_rounded, aes(x = Discharge_TAVI_year, y = rehab_rate_per_10000_days, group = 1)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_text(aes(label = sprintf("%.2f", rehab_rate_per_10000_days)), vjust = -1, size = 3.5, color = "black") +
  labs(
    x = "Year",
    y = "Cardiac Rehab Rate per 10,000 Person-Days",
    title = "Rehab Rate per Year",
    subtitle = "Number of Cardiac Rehab Exposures per 10,000 Follow-Up Days"
  ) +
  scale_y_continuous(limits = c(0, 4.0), expand = c(0, 0)) +  # <- this sets y-min to 0
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)
  )





###################################################################################################################################





