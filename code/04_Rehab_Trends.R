##########################################################################################################################################################################
#####    OBJECTIVE 1a: Temporal trends in cardiac rehabilitation attendance post TAVI procedure before and after the Covid-19 pandemic, across England (2019 to 2023)  ####
##########################################################################################################################################################################


#Option to pull data or refresh from data cleaning
clean_cohort <- read.csv("clean_cohort.csv")
View(clean_cohort)
##################################################################################################################################################################
#####explore the volume of TAVI patients, who attended a cardiac rehabilitation program within 180 days post discharge from the TAVI procedure######################
####################################################################################################################################################################
install.packages(c("dplyr", "ggplot2", "lubridate", "tidyr", "colorspace"))

library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(colorspace)


# Step 1: Convert TAVI discharge date to Date format and create new variable for "TAVI discharge date year-month"  
clean_cohort <- clean_cohort %>%
  mutate(tavi_cips_disdate = as.Date(tavi_cips_disdate),
         Discharge_TAVI_year_month = floor_date(tavi_cips_disdate, "month"))

# Step 2: Aggregate data for all TAVI patients. Calculate the counts for all TAVI patients per month (based on TAVI discharge date)
all_tavi_cohort <- clean_cohort %>%
  group_by(Discharge_TAVI_year_month) %>%
  summarise(all_tavi_cohort = n(), .groups = "drop")



# Step 3: Aggregate data for those attending cardiac rehab. Calculate the counts for patients who attended cardiac rehab per month (based on TAVI discharge date)
cardiac_rehab_exposure <- clean_cohort %>%
  filter(exp_cardiac_rehab_6m == 1) %>%
  group_by(Discharge_TAVI_year_month) %>%
  summarise(cardiac_rehab_exposure = n(), .groups = "drop")




# Step 4: Merge the two datasets
temporal_data <- left_join(all_tavi_cohort, cardiac_rehab_exposure, by = "Discharge_TAVI_year_month") %>%
  replace_na(list(cardiac_rehab_exposure = 0))  # Replace NA with 0 for those not attending rehab



#Step 5: Create "Discharge_TAVI_year_quarter" column to format the graph by quarters. Using the discharge date column from "Discharge_TAVI__year_month" 
library(dplyr)

# Add discharge_TAVI_year_quarter column to clean cohort
clean_cohort <- clean_cohort %>%
  mutate(
    Discharge_TAVI_year_quarter = paste0(
      format(Discharge_TAVI_year_month, "%Y"),
      "-Q",
      (as.numeric(format(Discharge_TAVI_year_month, "%m")) - 1) %/% 3 + 1
    )
  )


# Add discharge_TAVI_year_quarter column to temporal_data
temporal_data <- temporal_data %>%
  mutate(
    Discharge_TAVI_year_quarter = paste0(
      format(Discharge_TAVI_year_month, "%Y"),
      "-Q",
      (as.numeric(format(Discharge_TAVI_year_month, "%m")) - 1) %/% 3 + 1
    )
  )




#Step 6: calculate the sum of values for each quarter from the temporal_data dataframe
# Summarize the data to calculate quarterly totals
quarterly_data_volumes <- temporal_data %>%
  group_by(Discharge_TAVI_year_quarter) %>%
  summarise(
    total_tavi_cohort = round(sum(all_tavi_cohort, na.rm = TRUE) / 5) * 5,
    total_cardiac_rehab_exposure = round(sum(cardiac_rehab_exposure, na.rm = TRUE) / 5) * 5
  )





# Create three summary tables (monthly, quarterly and yearly) and round the aggregated counts to the nearest 5 before exporting to CSV files:
# Load necessary libraries
library(dplyr)
library(lubridate)

# ----- 1) MONTHLY SUMMARY -----

monthly_summary_volumes <- temporal_data %>%
  group_by(Discharge_TAVI_year_month) %>%
  summarise(
    All_TAVI_Patients = sum(all_tavi_cohort, na.rm = TRUE),
    Cardiac_Rehab_Exposure = sum(cardiac_rehab_exposure, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # Round counts to nearest 5
  mutate(
    All_TAVI_Patients = round(All_TAVI_Patients / 5) * 5,
    Cardiac_Rehab_Exposure = round(Cardiac_Rehab_Exposure / 5) * 5,
    # Suppress values under 10
    All_TAVI_Patients = ifelse(All_TAVI_Patients < 10, NA, All_TAVI_Patients),
    Cardiac_Rehab_Exposure = ifelse(Cardiac_Rehab_Exposure < 10, NA, Cardiac_Rehab_Exposure)
  )


# Convert the monthly_summary_volumes' Discharge_TAVI_year_month to character
monthly_summary_volumes <- monthly_summary_volumes %>%
  mutate(Discharge_TAVI_year_month = as.character(Discharge_TAVI_year_month))

# Create a total row for monthly data
monthly_total_volumes <- monthly_summary_volumes %>%
  summarise(
    Discharge_TAVI_year_month = "Total",
    All_TAVI_Patients = round(sum(All_TAVI_Patients, na.rm = TRUE) / 5) * 5,
    Cardiac_Rehab_Exposure = round(sum(Cardiac_Rehab_Exposure, na.rm = TRUE) / 5) * 5
  )

# Create the monthly average row (excluding NA rows from suppression)
monthly_avg_volumes <- monthly_summary_volumes %>%
  summarise(
    Discharge_TAVI_year_month = "Monthly Average",
    All_TAVI_Patients = round(mean(All_TAVI_Patients, na.rm = TRUE) / 5) * 5,
    Cardiac_Rehab_Exposure = round(mean(Cardiac_Rehab_Exposure, na.rm = TRUE) / 5) * 5
  )

# Combine everything
monthly_summary_volumes <- bind_rows(monthly_summary_volumes, monthly_total_volumes, monthly_avg_volumes)
print(monthly_summary_volumes)

# Save the rounded table as a CSV file for extraction
output_file_monthly_summary_volumes <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/rates and vols/monthly_summary_volumes.csv"
write.csv(monthly_summary_volumes, file = output_file_monthly_summary_volumes, row.names = TRUE)





# ----- 2) QUARTERLY SUMMARY -----

quarterly_summary_volumes <- temporal_data %>%
  # Create a date column from the month string and a new Quarter variable
  mutate(
    date_col = as.Date(paste0(Discharge_TAVI_year_month, "-01")),
    Quarter = paste0(year(date_col), " Q", quarter(date_col))
  ) %>%
  group_by(Quarter) %>%
  summarise(
    All_TAVI_Patients = sum(all_tavi_cohort, na.rm = TRUE),
    Cardiac_Rehab_Exposure = sum(cardiac_rehab_exposure, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # Round counts to nearest 5
  mutate(
    All_TAVI_Patients = round(All_TAVI_Patients / 5) * 5,
    Cardiac_Rehab_Exposure = round(Cardiac_Rehab_Exposure / 5) * 5
  )


# Create a total row for quarterly data
quarterly_total_summary <- quarterly_summary_volumes %>%
  summarise(
    Quarter = "Total",
    All_TAVI_Patients = round(sum(All_TAVI_Patients, na.rm = TRUE) / 5) * 5,
    Cardiac_Rehab_Exposure = round(sum(Cardiac_Rehab_Exposure, na.rm = TRUE) / 5) * 5
  )

# ---- Add quarterly average row ----
quarterly_avg_summary <- quarterly_summary_volumes %>%
  summarise(
    Quarter = "Quarterly Average",
    All_TAVI_Patients = round(mean(All_TAVI_Patients, na.rm = TRUE) / 5) * 5,
    Cardiac_Rehab_Exposure = round(mean(Cardiac_Rehab_Exposure, na.rm = TRUE) / 5) * 5
  )

# Combine the summary with the total row
quarterly_summary_volumes <- bind_rows(quarterly_summary_volumes, quarterly_total_summary, quarterly_avg_summary)
print(quarterly_summary_volumes)

# Save the rounded table as a CSV file for extraction
output_file_quarterly_summary_volumes <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/rates and vols/quarterly_summary_volumes.csv"
write.csv(quarterly_summary_volumes, file = output_file_quarterly_summary_volumes, row.names = TRUE)




# ----- 3) YEARLY SUMMARY -----

yearly_summary_volumes <- temporal_data %>%
  # Create a date column and extract the year
  mutate(
    date_col = as.Date(paste0(Discharge_TAVI_year_month, "-01")),
    Year = as.character(year(date_col))
  ) %>%
  group_by(Year) %>%
  summarise(
    All_TAVI_Patients = sum(all_tavi_cohort, na.rm = TRUE),
    Cardiac_Rehab_Exposure = sum(cardiac_rehab_exposure, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # Round counts to nearest 5
  mutate(
    All_TAVI_Patients = round(All_TAVI_Patients / 5) * 5,
    Cardiac_Rehab_Exposure = round(Cardiac_Rehab_Exposure / 5) * 5
  )

# Create a total row for yearly data
yearly_total_volumes <- yearly_summary_volumes %>%
  summarise(
    Year = "Total",
    All_TAVI_Patients = round(sum(All_TAVI_Patients, na.rm = TRUE) / 5) * 5,
    Cardiac_Rehab_Exposure = round(sum(Cardiac_Rehab_Exposure, na.rm = TRUE) / 5) * 5
  )

# ---- Add yearly average row ----
yearly_avg_summary <- yearly_summary_volumes %>%
  summarise(
    Year = "Yearly Average",
    All_TAVI_Patients = round(mean(All_TAVI_Patients, na.rm = TRUE) / 5) * 5,
    Cardiac_Rehab_Exposure = round(mean(Cardiac_Rehab_Exposure, na.rm = TRUE) / 5) * 5
  )


# Combine the summary with the total row
yearly_summary_volumes <- bind_rows(yearly_summary_volumes, yearly_total_volumes, yearly_avg_summary)
print(yearly_summary_volumes)

# Save the rounded table as a CSV file for extraction
output_file_yearly_summary_volumes <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/rates and vols/yearly_summary_volumes.csv"
write.csv(yearly_summary_volumes, file = output_file_yearly_summary_volumes, row.names = TRUE)




### PLOT THE VOLUMES ###

# Create the bar graph - monthly volumes 
install.packages("farver")
library(farver)

library(ggplot2)

#  code for the bar plot
ggplot(temporal_data, aes(x = Discharge_TAVI_year_month)) +
  geom_bar(aes(y = all_tavi_cohort, fill = "All TAVI Patients"),
           stat = "identity", position = "dodge", color = "black") +
  geom_bar(aes(y = cardiac_rehab_exposure, fill = "Cardiac Rehab Exposure"),
           stat = "identity", position = "dodge", color = "black") +
  # Adding COVID period highlight
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-31"),
           ymin = 0, ymax = Inf, alpha = 0.2, fill = "red") +
  scale_fill_manual(values = c("All TAVI Patients" = "darkgreen", 
                               "Cardiac Rehab Exposure" = "darkblue")) +
  labs(x = "Time (Year-Month)",
       y = "Number of Patients",
       fill = "Group",
       title = "Monthly Temporal Trends in Cardiac Rehabilitation Exposure",
       subtitle = "COVID period (March 2020 - December 2021) highlighted in red") +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%b") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))



##show each month on x-axis
library(ggplot2)
library(dplyr)
library(lubridate)

# Example: Convert YYYY-MM to a Date object for the 1st of the month
temporal_data <- temporal_data %>%
  mutate(date_col = as.Date(paste0(Discharge_TAVI_year_month, "-01")))

ggplot(temporal_data, aes(x = date_col)) +
  geom_bar(aes(y = all_tavi_cohort, fill = "All TAVI Patients"),
           stat = "identity",
           position = "dodge",
           color = "black") +
  geom_bar(aes(y = cardiac_rehab_exposure, fill = "Cardiac Rehab Exposure"),
           stat = "identity",
           position = "dodge",
           color = "black") +
  # Highlight COVID period with a red rectangle
  annotate(
    "rect",
    xmin = as.Date("2020-03-01"),    # Start of COVID period
    xmax = as.Date("2021-12-31"),    # End of COVID period
    ymin = 0,
    ymax = Inf,
    alpha = 0.2,
    fill = "red"
  ) +
  scale_fill_manual(
    values = c(
      "All TAVI Patients" = "darkgreen",
      "Cardiac Rehab Exposure" = "darkblue"
    )
  ) +
  labs(
    x = "Month",
    y = "Number of Patients",
    fill = "Group",
    subtitle = "COVID period (March 2020 - December 2021) highlighted in red",
    title = "Monthly Trends in Cardiac Rehabilitation Exposure"
  ) +
  # Force a tick mark for each month and label only the month name
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b %Y",   # "%b" for abbreviated month (Jan, Feb, etc.)
    expand = c(0, 0)
  ) +
  # Optionally rotate labels to prevent overlap
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )




#updated example that shades the UK lockdown periods using actual lockdown dates, and labels each month on the x-axis with the abbreviated month name. In this example, three UK lockdown periods are defined:

#Lockdown 1: 23 March 2020 to 13 may 2020

#Lockdown 2: 5 November 2020 to 20 December 2020

#Lockdown 3: 6 January 2021 to 29 March 2021

# Convert the "YYYY-MM" column to a proper Date (assuming the 1st of the month)
temporal_data <- temporal_data %>%
  mutate(date_col = as.Date(paste0(Discharge_TAVI_year_month, "-01")))

temporal_data_rounded <- temporal_data %>%
  mutate(
    all_tavi_cohort = round(all_tavi_cohort / 5) * 5,
    cardiac_rehab_exposure = round(cardiac_rehab_exposure / 5) * 5,
    # Create label columns where values < 10 are shown as "<10"
    all_tavi_label = ifelse(all_tavi_cohort < 10, "<10", as.character(all_tavi_cohort)),
    rehab_label = ifelse(cardiac_rehab_exposure < 10, "<10", as.character(cardiac_rehab_exposure))
  )

# Define UK lockdown periods in a data frame
lockdowns <- data.frame(
  start = as.Date(c("2020-03-23", "2020-11-05", "2021-01-06")),
  end   = as.Date(c("2020-05-13", "2020-12-20", "2021-03-29"))
)

cardiac_rehab_monthly_volumes_plot <- ggplot(temporal_data_rounded, aes(x = date_col)) +
  # Add lockdown shading for each lockdown period
  geom_rect(data = lockdowns,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "red", alpha = 0.2, inherit.aes = FALSE) +
  # Bar for all TAVI patients
  geom_bar(aes(y = all_tavi_cohort, fill = "All TAVI Patients"),
           stat = "identity", position = "dodge", color = "black") +
  # Bar for cardiac rehab exposure
  geom_bar(aes(y = cardiac_rehab_exposure, fill = "Cardiac Rehab Exposure"),
           stat = "identity", position = "dodge", color = "black") +
  
  # Volume labels for All TAVI
  geom_text(aes(y = all_tavi_cohort, label = all_tavi_cohort),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 2.0, color = "black") +
  
  # Volume labels for Rehab exposure (slightly adjusted position)
  geom_text(aes(y = cardiac_rehab_exposure, label = rehab_label),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 1.9, color = "black") +
  
  scale_fill_manual(values = c(
    "All TAVI Patients" = "darkgrey",
    "Cardiac Rehab Exposure" = "darkblue"
  )) +
  labs(
    x = "Month",
    y = "Number of Patients",
    fill = "Group",
    subtitle = "COVID lockdowns highlighted in red",
    title = "Monthly Cardiac Rehabilitation Volumes post TAVI (January 2018 - March 2023)"
  ) +
  # Force a tick mark for each month and label with abbreviated month names
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b %y",  # abbreviated month names (e.g., Jan, Feb, Mar)
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 600),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(face = "bold", size = 8),
    plot.subtitle = element_text(size = 8),
    legend.position = "right",             
    legend.text = element_text(size = 8),  # smaller legend labels
    legend.title = element_text(size = 9)  # smaller legend title
  )


# Save as PDF
ggsave("cardiac_rehab_monthly_volumes_plot_figure1_manuscript.pdf", plot = cardiac_rehab_monthly_volumes_plot, width = 10, height = 6, units = "in")




# updated bar graph with time in quarters
# Create the bar graph
ggplot(quarterly_data_volumes, aes(x = Discharge_TAVI_year_quarter)) +
  # Add bars for all TAVI patients
  geom_bar(aes(y = total_tavi_cohort, fill = "All TAVI Patients"),
           stat = "identity", position = "dodge", color = "black",
           width = 0.7, alpha = 0.6) +
  # Add bars for cardiac rehab exposure
  geom_bar(aes(y = total_cardiac_rehab_exposure, fill = "Cardiac Rehab Exposure"),
           stat = "identity", position = "dodge", color = "black",
           width = 0.7, alpha = 0.9) +
  # Add labels above the bars
  geom_text(aes(y = total_tavi_cohort + 50, label = total_tavi_cohort),
            color = "darkgreen", size = 5, fontface = "bold", vjust = 0) +
  geom_text(aes(y = total_cardiac_rehab_exposure + 50, label = total_cardiac_rehab_exposure),
            color = "darkblue", size = 5, fontface = "bold", vjust = 0) +
  # Customize fill colors
  scale_fill_manual(values = c("All TAVI Patients" = "lightgreen",
                               "Cardiac Rehab Exposure" = "darkblue")) +
  labs(x = "Time (Year-Quarter)",
       y = "Number of Patients",
       title = "Quarterly Trends in Cardiac Rehabilitation Exposure",
       subtitle = "Summed values per quarter",
       fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))


##############################################################################################
###        UPDATED BAR GRAPH WITH IMPROVED FORMATTING SHOWING COVID period          ##########
##############################################################################################

# Summarize the data to calculate quarterly totals and round counts to the nearest 5
quarterly_data_volumes <- temporal_data %>%
  group_by(Discharge_TAVI_year_quarter) %>%
  summarise(
    total_tavi_cohort = round(sum(all_tavi_cohort, na.rm = TRUE) / 5) * 5,
    total_cardiac_rehab_exposure = round(sum(cardiac_rehab_exposure, na.rm = TRUE) / 5) * 5
  )

# Create the bar graph
ggplot(quarterly_data_volumes, aes(x = Discharge_TAVI_year_quarter)) +
  # Add COVID period highlight (March 2020 to December 2021)
  annotate("rect",
           xmin = which(quarterly_data_volumes$Discharge_TAVI_year_quarter == "2020-Q1"),
           xmax = which(quarterly_data_volumes$Discharge_TAVI_year_quarter == "2021-Q4"),
           ymin = 0, ymax = Inf, fill = "red", alpha = 0.1) +
  # Add bars for all TAVI patients
  geom_bar(aes(y = total_tavi_cohort, fill = "All TAVI Patients"),
           stat = "identity", position = "dodge", color = "black",
           width = 0.7, alpha = 0.6) +
  # Add bars for cardiac rehab exposure
  geom_bar(aes(y = total_cardiac_rehab_exposure, fill = "Cardiac Rehab Exposure"),
           stat = "identity", position = "dodge", color = "black",
           width = 0.7, alpha = 0.9) +
  # Add labels above the bars
  geom_text(aes(y = total_tavi_cohort + 50, label = total_tavi_cohort),
            color = "darkgreen", size = 5, fontface = "bold", vjust = 0) +
  geom_text(aes(y = total_cardiac_rehab_exposure + 50, label = total_cardiac_rehab_exposure),
            color = "darkblue", size = 5, fontface = "bold", vjust = 0) +
  # Customize fill colors
  scale_fill_manual(values = c("All TAVI Patients" = "lightgreen",
                               "Cardiac Rehab Exposure" = "darkblue")) +
  labs(x = "Time (Year-Quarter)",
       y = "Number of Patients",
       title = "Quarterly Trends in Cardiac Rehabilitation Exposure",
       subtitle = "Summed values per quarter with COVID period highlighted",
       fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))



###############################################################################

##############################################################################################
###        UPDATED BAR GRAPH WITH IMPROVED FORMATTING FOR CR GROUP          ##########
##############################################################################################

# Summarize the data to calculate quarterly totals and round counts to the nearest 5
quarterly_data_volumes <- temporal_data %>%
  group_by(Discharge_TAVI_year_quarter) %>%
  summarise(
    total_tavi_cohort = round(sum(all_tavi_cohort, na.rm = TRUE) / 5) * 5,
    total_cardiac_rehab_exposure = round(sum(cardiac_rehab_exposure, na.rm = TRUE) / 5) * 5
  )

library(ggplot2)

ggplot(quarterly_data_volumes, aes(x = Discharge_TAVI_year_quarter)) +
  # Add COVID period highlight
  annotate("rect", 
           xmin = which(quarterly_data_volumes$Discharge_TAVI_year_quarter == "2020-Q1"),
           xmax = which(quarterly_data_volumes$Discharge_TAVI_year_quarter == "2021-Q4"),
           ymin = 0, ymax = Inf, fill = "red", alpha = 0.1) +
  
  # Add bars for TAVI patients
  geom_bar(aes(y = total_tavi_cohort, fill = "All TAVI Patients"),
           stat = "identity", position = "dodge", color = "black", width = 0.7, alpha = 0.8) +
  
  # Add scaled bars for Cardiac Rehab Exposure
  geom_bar(aes(y = total_cardiac_rehab_exposure, fill = "Cardiac Rehab Exposure"),
           stat = "identity", position = "dodge", color = "black", width = 0.5, alpha = 0.9) +
  
  # Add labels above the bars
  geom_text(aes(y = total_tavi_cohort, label = total_tavi_cohort), 
            color = "darkgreen", size = 5, fontface = "bold", vjust = -0.5) +
  geom_text(aes(y = total_cardiac_rehab_exposure + 50, label = total_cardiac_rehab_exposure), 
            color = "#002060", size = 5, fontface = "bold", vjust = -0.5) +
  
  # Customize fill colors
  scale_fill_manual(values = c("All TAVI Patients" = "darkgreen", 
                               "Cardiac Rehab Exposure" = "#002060")) +
  
  # Customize axes and titles
  labs(x = "Time (Year-Quarter)",
       y = "Number of Patients",
       title = "Quarterly Trends in Cardiac Rehabilitation Exposure",
       subtitle = "Summed values per quarter with COVID period highlighted",
       fill = "Group") +
  
  # Improve theme for clarity
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))








##############################################################################################################################################


##############################################################################################
###        UPDATED DATA for ANNUAL VOLUMES                                          ##########
##############################################################################################

#Updated Code for Annual Aggregation



# Calculate the sum of values for each year
# Summarize the data to calculate annual totals
annual_data_volumes <- temporal_data %>%
  mutate(Discharge_TAVI_year = format(Discharge_TAVI_year_month, "%Y")) %>%
  group_by(Discharge_TAVI_year) %>%
  summarise(
    total_tavi_cohort = round(sum(all_tavi_cohort, na.rm = TRUE) / 5) * 5,
    total_cardiac_rehab_exposure = round(sum(cardiac_rehab_exposure, na.rm = TRUE) / 5) * 5
  )


# Updated bar graph with time in years
ggplot(annual_data_volumes, aes(x = Discharge_TAVI_year)) +
  # Add bars for all TAVI patients
  geom_bar(aes(y = total_tavi_cohort, fill = "All TAVI Patients"), 
           stat = "identity", position = "dodge", color = "black") +
  
  # Add bars for cardiac rehab exposure
  geom_bar(aes(y = total_cardiac_rehab_exposure, fill = "Cardiac Rehab Exposure"), 
           stat = "identity", position = "dodge", color = "black") +
  
  # Add labels above bars
  geom_text(aes(y = total_tavi_cohort + 100, label = total_tavi_cohort),
            color = "darkgreen", size = 5, fontface = "bold", vjust = 0) +
  geom_text(aes(y = total_cardiac_rehab_exposure + 100, label = total_cardiac_rehab_exposure),
            color = "#002060", size = 5, fontface = "bold", vjust = 0) +
  
  # Customize fill colors
  scale_fill_manual(values = c("All TAVI Patients" = "darkgreen",
                               "Cardiac Rehab Exposure" = "#002060")) +
  
  labs(x = "Year", 
       y = "Number of Patients", 
       title = "Annual Trends in Cardiac Rehabilitation Exposure",
       subtitle = "Summed values per year",
       fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))



##To calculate the annual growth rate of TAVI patients and cardiac rehab exposure from 2018 to 2022, we need to:

##Compute the year-over-year percentage change in the total number of TAVI patients and cardiac rehab exposures.

# Filter data for years 2018 to 2022
annual_growth_data <- annual_data_volumes %>%
  filter(Discharge_TAVI_year >= 2018 & Discharge_TAVI_year <= 2022) %>%
  arrange(Discharge_TAVI_year) %>%
  mutate(
    growth_tavi_cohort = (total_tavi_cohort / lag(total_tavi_cohort) - 1) * 100,
    growth_cardiac_rehab_exposure = (total_cardiac_rehab_exposure / lag(total_cardiac_rehab_exposure) - 1) * 100
  )

# Display the growth table
print(annual_growth_data)   
View(annual_growth_data)    


# Save the rounded table as a CSV file for extraction
output_file_year_on_year_growth_data <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/rates and vols/annual_growth_data.csv"
write.csv(annual_growth_data, file = output_file_year_on_year_growth_data, row.names = FALSE)



# Extract values for 2018 and 2022
tavi_2018 <- annual_data_volumes %>% filter(Discharge_TAVI_year == 2018) %>% pull(total_tavi_cohort)
tavi_2022 <- annual_data_volumes %>% filter(Discharge_TAVI_year == 2022) %>% pull(total_tavi_cohort)

rehab_2018 <- annual_data_volumes %>% filter(Discharge_TAVI_year == 2018) %>% pull(total_cardiac_rehab_exposure)
rehab_2022 <- annual_data_volumes %>% filter(Discharge_TAVI_year == 2022) %>% pull(total_cardiac_rehab_exposure)

# Compute percentage change
tavi_change <- ((tavi_2022 - tavi_2018) / tavi_2018) * 100
rehab_change <- ((rehab_2022 - rehab_2018) / rehab_2018) * 100

# Display the results
growth_2018_to_2022_vols <- data.frame(
  Metric = c("TAVI Patients", "Cardiac Rehab Exposure"),
  Change_Percentage = c(tavi_change, rehab_change)
)

print(growth_2018_to_2022_vols)

# Save the rounded table as a CSV file for extraction
output_file_growth_2018_to_2022_vols <- "D:/PhotonUser/My Files/Home Folder/Justin Braver/Outputs/rates and vols/growth_2018_to_2022_vols.csv"
write.csv(growth_2018_to_2022_vols, file = output_file_growth_2018_to_2022_vols, row.names = FALSE)





##############################################################################################################################################

##To calculate the quarterly growth rate of TAVI patients and cardiac rehab exposure from 2018 to 2023, we need to:

##Compute the quarter-over-quarter percentage change in the total number of TAVI patients and cardiac rehab exposures.

# Filter data for years 2018 to 2022
quartelry_growth_data <- quarterly_data_volumes %>%
  arrange(Discharge_TAVI_year_quarter) %>%
  mutate(
    growth_tavi_cohort = (total_tavi_cohort / lag(total_tavi_cohort) - 1) * 100,
    growth_cardiac_rehab_exposure = (total_cardiac_rehab_exposure / lag(total_cardiac_rehab_exposure) - 1) * 100
  )

# Display the growth table
print(quartelry_growth_data)   
View(quartelry_growth_data)    





##############################################################################################################################################

