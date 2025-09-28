


# script to connect RStudio to Databricks

rm(list = ls())



# load required packages
library(odbc)
library(DBI)
library(dbplyr)
library(dplyr)
library(lubridate)
library(pbs)

# create connection to databricks
con <- DBI::dbConnect(odbc::odbc(),
                      "Databricks",
                      Password = rstudioapi::askForPassword("Enter token:"))


## connecting to data tables

# collaboration database
dbc = "dsa_391419_j3w9t_collab"

# connect to CCU090 cohort table
cohort_tbl <- tbl(con, in_schema(dbc, "ccu090_01_out_combined"))

# read in the whole cohort data into dataframe
cohort_dat <- cohort_tbl %>%
  collect()

# or read in summarised data e.g. the number of records/patients per month, by cardiac rehab status
cohort_monthly_summary <- cohort_tbl %>%
  mutate(adm_ym = (year(`5_06_ADMISSION_DATE_FOR_PRO_1ST_HOSPITAL`)*100)+month(`5_06_ADMISSION_DATE_FOR_PRO_1ST_HOSPITAL`)) %>%
  group_by(adm_ym, exp_cardiac_rehab_6m) %>%
  summarise(recs = n(),
            pats = n_distinct(person_id)) %>%
  arrange(adm_ym) %>%
  collect()



